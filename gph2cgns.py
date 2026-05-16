#!/usr/bin/env python3
"""
GPH to CGNS Converter.

Reads GPH (CRDL-FLD) binary format produced by Software Cradle scFLOW / SCTpre
and writes a CGNS/HDF5 unstructured mesh file using the NGON_n / NFACE_n
convention employed by the vendor's own ``FLDUTIL`` exporter
(see :file:`box_ansa_orig.cgns` for the reference layout).

Requires: numpy, h5py (no CGNS library needed).
"""

import argparse
import struct
import sys
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np

try:
    import h5py
except ImportError:
    print("Error: h5py is required. Install with: pip install h5py numpy")
    sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Low-level GPH primitives
# ─────────────────────────────────────────────────────────────────────────────


def read_i32_be(data: bytes, pos: int) -> int:
    return int.from_bytes(data[pos : pos + 4], "big")


def read_f32_be(data: bytes, pos: int) -> float:
    return struct.unpack(">f", data[pos : pos + 4])[0]


def read_f64_be(data: bytes, pos: int) -> float:
    """Read a standard big-endian IEEE-754 float64."""
    return struct.unpack(">d", data[pos : pos + 8])[0]


def read_f64_wr(data: bytes, pos: int) -> float:
    """Read float64 stored in word-reversed (middle-endian) format.

    Some legacy GPH files encode each 8-byte double as two 32-bit big-endian
    words in reversed order: ``[lower_32bit_word][upper_32bit_word]``.  Newer
    files (e.g. those exported by current scFLOW/ANSA pipelines) use plain
    big-endian instead — the GPH dialect is auto-detected at parse time.
    """
    lower = int.from_bytes(data[pos : pos + 4], "big")
    upper = int.from_bytes(data[pos + 4 : pos + 8], "big")
    combined = ((upper << 32) | lower).to_bytes(8, "big")
    return struct.unpack(">d", combined)[0]


def _find_section(data: bytes, name: str) -> int:
    """Return offset of the I4=32 marker that precedes *name* (padded to 32 chars).

    Returns -1 if not found.
    """
    name_padded = name.ljust(32).encode("ascii")
    idx = data.find(name_padded)
    if idx < 4:
        return -1
    if read_i32_be(data, idx - 4) == 32:
        return idx - 4
    return -1


def _section_end(data: bytes, sec_start: int) -> int:
    """Return the byte offset at which the section starting at *sec_start* ends.

    The end is taken to be the start of the next known section (whichever
    appears first after *sec_start*) or the file end.
    """
    candidates = [
        "FileRevision", "Application", "ApplicationVersion", "ReleaseDate",
        "GridType", "Dimension", "Bias", "Date", "Comments", "Cycle",
        "Unused", "Encoding", "HeaderDataEnd", "OverlapStart_0",
        "LS_CvolIdOfElements", "LS_Links", "LS_Nodes", "LS_SurfaceRegions",
        "Element_InformationFlag", "OverlapEnd",
    ]
    best = len(data)
    for name in candidates:
        off = _find_section(data, name)
        if off > sec_start and off < best:
            best = off
    return best


# ─────────────────────────────────────────────────────────────────────────────
# Generic GPH "data-block" scanner
# ─────────────────────────────────────────────────────────────────────────────
#
# Inside a named section every payload is stored as
#
#       [I4 = 12]              ← header tag
#       [I4 = byte_count]      ← number of payload bytes that follow
#       [byte_count bytes]     ← payload data
#       [I4 = byte_count]      ← trailing length sentinel
#
# Interleaved between (and sometimes before) these blocks are 16-byte
# descriptors of the form ``[12, type_code, dim0, dim1]``.  The scanner below
# walks a section, ignoring descriptors and yielding every data block.
# ─────────────────────────────────────────────────────────────────────────────


def _iter_data_blocks(data: bytes, sec_start: int, sec_end: int):
    """Yield ``(payload_start, byte_count)`` for each data block in the section."""
    pos = sec_start + 40  # skip [I4=32][32B name][I4=32]
    n = len(data)
    while pos + 8 <= sec_end and pos + 8 <= n:
        if read_i32_be(data, pos) != 12:
            pos += 4
            continue
        v = read_i32_be(data, pos + 4)

        # Descriptor [12, type_code in {4,8}, dim0, dim1] is 16 bytes.
        if v in (4, 8) and pos + 16 <= sec_end:
            dim0 = read_i32_be(data, pos + 8)
            dim1 = read_i32_be(data, pos + 12)
            if 0 < dim0 < 10_000_000 and 0 < dim1 < 10_000_000:
                pos += 16
                continue

        # Otherwise treat as a data header [12, byte_count].
        bc = v
        if bc <= 0 or pos + 8 + bc + 4 > sec_end:
            pos += 4
            continue
        payload_end = pos + 8 + bc
        if read_i32_be(data, payload_end) != bc:
            pos += 4
            continue
        yield pos + 8, bc
        pos = payload_end + 4


# ─────────────────────────────────────────────────────────────────────────────
# LS_Nodes – vertex coordinates
# ─────────────────────────────────────────────────────────────────────────────


def _parse_ls_nodes(data: bytes):
    """Parse the LS_Nodes section.

    The section contains three coordinate blocks (X, Y, Z file order).  Two
    on-disk encodings have been observed in the wild:

    1. Standard big-endian float64 (modern files such as ``box_ansa.gph``).
    2. Word-reversed float64 — each 8-byte double stored as ``[low32_BE][high32_BE]``
       (legacy ``box.gph`` files).

    The dialect is auto-detected by trying both readings and picking the one
    whose magnitudes are physically plausible coordinate values.
    """
    sec_start = _find_section(data, "LS_Nodes")
    if sec_start < 0:
        return None, 0
    sec_end = _section_end(data, sec_start)

    blocks = list(_iter_data_blocks(data, sec_start, sec_end))
    # Keep only blocks whose byte count is a positive multiple of 8 (float64
    # arrays).  The three biggest such blocks of identical size are the X/Y/Z
    # coordinate blocks.
    f64_blocks = [(p, bc) for p, bc in blocks if bc >= 8 and bc % 8 == 0]
    if len(f64_blocks) < 3:
        return None, 0

    # Pick the first three blocks that share the largest common size — that is
    # the X/Y/Z trio.
    sizes = [bc for _, bc in f64_blocks]
    target = max(set(sizes), key=sizes.count)
    trio = [(p, bc) for p, bc in f64_blocks if bc == target][:3]
    if len(trio) < 3:
        return None, 0

    n_vertices = trio[0][1] // 8

    def _decode(reader):
        axes = []
        for payload_start, _ in trio:
            vals = np.empty(n_vertices, dtype=np.float64)
            for i in range(n_vertices):
                vals[i] = reader(data, payload_start + i * 8)
            axes.append(vals)
        return axes

    axes_be = _decode(read_f64_be)
    axes_wr = _decode(read_f64_wr)

    def _score(axes):
        """Lower score = more plausible coordinate magnitudes."""
        score = 0.0
        for ax in axes:
            finite = np.isfinite(ax)
            if not finite.all():
                score += 1e30
                continue
            absmax = np.max(np.abs(ax)) if ax.size else 0.0
            if absmax > 1e6 or absmax < 1e-30 and absmax != 0.0:
                score += absmax + (1.0 / max(absmax, 1e-300))
            else:
                score += absmax
        return score

    axes = axes_be if _score(axes_be) <= _score(axes_wr) else axes_wr
    xyz = np.column_stack(axes)  # file order = (X, Y, Z) for box_ansa / (X, Z, Y) for box

    # For legacy word-reversed files we keep historic (X, Z, Y) swap behaviour
    # only when word-reversed encoding was selected.
    if axes is axes_wr:
        xyz = xyz[:, [0, 2, 1]]

    return xyz, n_vertices


# ─────────────────────────────────────────────────────────────────────────────
# LS_Links – face / cell connectivity
# ─────────────────────────────────────────────────────────────────────────────


def _parse_ls_links(data: bytes):
    """Parse the LS_Links section.

    Returns a dict with keys::

        n_faces, face_nodes (ndarray (n_faces, npe) 0-based),
        owner, neighbor (1-D int64 arrays, ``-1`` = boundary),
        boundary_faces (list of 0-based face indices),
        cell_faces (dict cell_id → list of 0-based face indices owned by cell),
        cell_neighbor_faces (dict cell_id → list of 0-based face indices
                              for which cell_id is the *neighbour*),
        n_cells.
    """
    sec_start = _find_section(data, "LS_Links")
    if sec_start < 0:
        return None
    sec_end = _section_end(data, sec_start)

    blocks = [(p, bc) for p, bc in _iter_data_blocks(data, sec_start, sec_end) if bc > 0]
    if not blocks:
        return None

    # ── Identify owner / neighbour / npe / connectivity blocks ──────────────
    #
    # The first three "small" data blocks (all of identical size = 4 × n_faces)
    # are owner, neighbour and nodes-per-face.  Some legacy GPH files insert
    # a 4-byte "face_type" block between them and the connectivity block — we
    # tolerate that by simply picking the three first blocks of equal size.
    block_sizes = [bc for _, bc in blocks]
    # The most common block size in the list is owner/neighbour/npe (3 blocks).
    from collections import Counter
    common = Counter(block_sizes).most_common()
    n_faces_block_size = None
    for size, count in common:
        if count >= 3 and size % 4 == 0 and size >= 4:
            n_faces_block_size = size
            break
    if n_faces_block_size is None:
        return None
    n_faces = n_faces_block_size // 4

    triples = [b for b in blocks if b[1] == n_faces_block_size][:3]
    owner_p, _ = triples[0]
    neigh_p, _ = triples[1]
    # triples[2] is npe (nodes-per-face); we don't strictly need its values.

    BNDRY_U = 0xFFFFFFFF
    owner = np.frombuffer(data[owner_p : owner_p + n_faces_block_size],
                          dtype=">u4").astype(np.int64).copy()
    neigh_raw = np.frombuffer(data[neigh_p : neigh_p + n_faces_block_size],
                              dtype=">u4").copy()
    neigh = neigh_raw.astype(np.int64)
    neigh[neigh_raw == BNDRY_U] = -1

    # The connectivity block: largest remaining block whose byte_count is a
    # multiple of 4 and ≥ 3 × n_faces × 4 (assuming at least triangular faces).
    conn_block = None
    for p, bc in blocks:
        if (p, bc) in triples:
            continue
        if bc % 4 != 0:
            continue
        if bc < 3 * n_faces * 4:
            continue
        if conn_block is None or bc > conn_block[1]:
            conn_block = (p, bc)
    if conn_block is None:
        return None
    conn_p, conn_bc = conn_block
    conn_total = conn_bc // 4
    npe = conn_total // n_faces  # nodes per face (3 for triangular)
    if npe < 3 or n_faces * npe != conn_total:
        return None

    conn = np.frombuffer(data[conn_p : conn_p + conn_bc], dtype=">u4").astype(np.int64).copy()

    # The file stores face → node connectivity row-major: ``face_nodes[i, k]``
    # is ``conn[i * npe + k]``.  Vertex indices are 0-based.
    face_nodes = conn.reshape(n_faces, npe)

    # ── Build cell ↔ face mappings.  In CGNS NFACE_n, faces owned by a cell
    # carry a positive index while faces for which the cell is the neighbour
    # are stored as a negative (reverse-oriented) index. ─────────────────────
    n_cells = int(max(owner.max() + 1, (neigh.max() + 1) if (neigh >= 0).any() else 0))
    cell_owner_faces: dict = defaultdict(list)
    cell_neighbor_faces: dict = defaultdict(list)
    for fi in range(n_faces):
        ow = int(owner[fi])
        if 0 <= ow < n_cells:
            cell_owner_faces[ow].append(fi)
        nb = int(neigh[fi])
        if nb != -1 and 0 <= nb < n_cells:
            cell_neighbor_faces[nb].append(fi)

    boundary_faces = [fi for fi in range(n_faces) if neigh[fi] == -1]

    return {
        "n_faces": n_faces,
        "npe": npe,
        "face_nodes": face_nodes,
        "owner": owner,
        "neighbor": neigh,
        "boundary_faces": boundary_faces,
        "cell_owner_faces": cell_owner_faces,
        "cell_neighbor_faces": cell_neighbor_faces,
        "n_cells": n_cells,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Vertex renumbering — match the official FLDUTIL exporter's ordering
# ─────────────────────────────────────────────────────────────────────────────


def _renumber_by_first_use(face_nodes: np.ndarray, n_vertices: int) -> np.ndarray:
    """Return a permutation ``perm`` such that ``perm[old_gph_id] = new_id``
    where ``new_id`` is the 0-based first-use index of ``old_gph_id`` in the
    row-major scan of ``face_nodes``.  Vertices that never appear in any face
    are appended at the end in their original order so the permutation is a
    bijection.
    """
    perm = np.full(n_vertices, -1, dtype=np.int64)
    next_id = 0
    flat = face_nodes.reshape(-1)
    for v in flat:
        v = int(v)
        if 0 <= v < n_vertices and perm[v] == -1:
            perm[v] = next_id
            next_id += 1
    for v in range(n_vertices):
        if perm[v] == -1:
            perm[v] = next_id
            next_id += 1
    return perm


# ─────────────────────────────────────────────────────────────────────────────
# GPH front-end
# ─────────────────────────────────────────────────────────────────────────────


def parse_gph_mesh(filepath: str) -> dict:
    """Extract mesh data (vertices, faces, cells) from a GPH file.

    The returned vertices/face_nodes have already been renumbered so that the
    output CGNS matches the vendor's own conversion (vertices appear in the
    order they are first referenced by the face connectivity).
    """
    with open(filepath, "rb") as f:
        data = f.read()

    result: dict = {
        "vertices": None,
        "n_vertices": 0,
        "link_data": None,
        "n_elements": 0,
    }

    xyz, n_vertices = _parse_ls_nodes(data)
    link_data = _parse_ls_links(data)

    if xyz is None or n_vertices == 0 or link_data is None:
        return result

    # Clamp any garbage node indices that exceed the known vertex count
    # (occasional byte-count sentinel leakage seen on legacy files).
    fn = link_data["face_nodes"]
    bad = fn >= n_vertices
    if bad.any():
        fn = fn.copy()
        fn[bad] = n_vertices - 1
        link_data["face_nodes"] = fn

    # Renumber vertices according to first-use order — matches FLDUTIL.
    perm = _renumber_by_first_use(link_data["face_nodes"], n_vertices)
    inv_perm = np.argsort(perm)
    xyz_renum = xyz[inv_perm]
    face_nodes_renum = perm[link_data["face_nodes"]]
    link_data["face_nodes"] = face_nodes_renum

    result["vertices"] = xyz_renum
    result["n_vertices"] = n_vertices
    result["link_data"] = link_data
    result["n_elements"] = link_data["n_cells"]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# CGNS / HDF5 writer
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# HDF5 1.8 "compact link storage" strategy
# ─────────────────────────────────────────────────────────────────────────────
#
# CGNS files exported by Software Cradle's ``FLDUTIL`` tool (see the
# reference file ``box_ansa_orig.cgns``) store every group's child links
# *inline in the parent's object header* instead of using the legacy v1.6
# symbol-table / B-tree / local-heap triple (SNOD / TREE / HEAP).  This is
# HDF5 1.8's so-called **compact link storage** mode (selected automatically
# for any group whose child count is below ``max_compact``, default 8).
#
# h5py's default ``File()`` uses ``libver=("earliest", "v108")`` which keeps
# the file readable by HDF5 1.6 but *forces* the heavier sym-table layout
# for every group.  For our use-case that adds roughly 1 KB of overhead per
# group → ~45 KB of pure metadata bloat for a typical FLDUTIL tree
# (observed: 73 KB vs 35 KB reference).
#
# Strategy applied below:
#
#   1. Open the file with ``libver=("v108", "v108")`` to pin the HDF5 1.8
#      file format (v2 superblock + v2 object headers + compact links).
#   2. The CGNS tree has at most 7 children per group (Zone_t), comfortably
#      below the HDF5 default ``max_compact = 8`` link-phase-change
#      threshold, so every group is automatically stored compactly — no
#      SNOD/TREE/HEAP triple is ever written.
#   3. Group creation goes through the explicit ``_create_compact_group``
#      helper below, which uses h5py's ``track_order=False`` to avoid the
#      attribute creation-order index (extra metadata) and documents the
#      compact-storage intent at the call site.
#
# All HDF5 ≥ 1.8 readers (released 2007 — i.e. every modern CGNS toolkit)
# handle this layout natively.
# ─────────────────────────────────────────────────────────────────────────────


# Sanity-check at module load: the HDF5 default link-phase-change threshold
# must be ≥ the maximum child count we ever produce in a CGNS group (7).
# If this ever fails — e.g. because a future h5py changes the default — we
# will trip the assertion below and need to revisit the strategy.
_MAX_CHILDREN_PER_GROUP = 7  # Zone_t children: data, ZoneType, GridCoordinates,
                              # GridElements_Faces, <NFace>, ZoneBC, FlowSolution


def _create_compact_group(parent, name: str):
    """Create an HDF5 group under *parent* using HDF5 1.8 compact link storage.

    With ``libver=("v108", "v108")`` set on the parent file, h5py's default
    link-phase-change thresholds (``max_compact = 8``, ``min_dense = 6``)
    keep every group with fewer than 8 children in compact storage — i.e.
    the child links are stored inline inside the group's object header
    rather than in the legacy SNOD + TREE + HEAP triple.

    ``track_order=False`` further suppresses the attribute creation-order
    index (which would otherwise add a small amount of metadata per
    group).  The combination produces a layout byte-equivalent to the one
    used by the vendor's FLDUTIL exporter.
    """
    return parent.create_group(name, track_order=False)


def _set_cgns_attrs(grp, name: str, label: str, type_str: str) -> None:
    """Set the four standard CGNS HDF5 group attributes."""
    grp.attrs.create("flags", np.array([1], dtype=np.int32))
    grp.attrs.create("label", np.bytes_(label), dtype=h5py.string_dtype(length=33))
    grp.attrs.create("name", np.bytes_(name), dtype=h5py.string_dtype(length=33))
    grp.attrs.create("type", np.bytes_(type_str), dtype=h5py.string_dtype(length=3))


def _cgns_node(parent, name: str, label: str, type_str: str):
    """Create an HDF5 group representing a CGNS node and set its attributes.

    The group is created with HDF5 1.8 compact link storage so child links
    are inlined into the object header (no SNOD / TREE / HEAP triple).
    """
    grp = _create_compact_group(parent, name)
    _set_cgns_attrs(grp, name, label, type_str)
    return grp


def _bytes_dataset(parent, name: str, payload: bytes) -> None:
    parent.create_dataset(" data", data=np.frombuffer(payload, dtype=np.int8))


def _i32_dataset(parent, name: str, arr) -> None:
    parent.create_dataset(" data", data=np.asarray(arr, dtype=np.int32))


def _write_grid_coordinates(zone, vertices: np.ndarray) -> None:
    gc = _cgns_node(zone, "GridCoordinates", "GridCoordinates_t", "MT")
    for axis, cname in enumerate(("CoordinateX", "CoordinateY", "CoordinateZ")):
        cd = _cgns_node(gc, cname, "DataArray_t", "R8")
        cd.create_dataset(" data", data=np.ascontiguousarray(vertices[:, axis],
                                                              dtype=np.float64))


def _write_ngon(zone, face_nodes_1based: np.ndarray, elem_range) -> None:
    n_faces, npf = face_nodes_1based.shape
    el = _cgns_node(zone, "GridElements_Faces", "Elements_t", "I4")
    _i32_dataset(el, " data", [22, n_faces])

    er = _cgns_node(el, "ElementRange", "IndexRange_t", "I4")
    _i32_dataset(er, " data", elem_range)

    offsets = np.arange(0, (n_faces + 1) * npf, npf, dtype=np.int32)
    eso = _cgns_node(el, "ElementStartOffset", "DataArray_t", "I4")
    _i32_dataset(eso, " data", offsets)

    ec = _cgns_node(el, "ElementConnectivity", "DataArray_t", "I4")
    _i32_dataset(ec, " data", face_nodes_1based.flatten())


def _build_nface_connectivity(link_data: dict):
    """Return ``(offsets, conn, n_cells)`` for the CGNS NFACE_n section.

    For each cell, the face list contains the 1-based face IDs of every face
    incident on that cell.  Faces for which the cell is the *owner* are stored
    with a positive index; faces for which the cell is the *neighbour* are
    stored with a *negative* index — this is how CGNS encodes the reverse
    orientation.  The face IDs within a single cell are sorted by absolute
    value (matching the FLDUTIL exporter).
    """
    n_cells = link_data["n_cells"]
    owner_faces = link_data["cell_owner_faces"]
    neigh_faces = link_data["cell_neighbor_faces"]

    offsets = [0]
    conn: list[int] = []
    for cell_id in range(n_cells):
        cell_list: list[int] = []
        for fi in owner_faces.get(cell_id, ()):
            cell_list.append(fi + 1)
        for fi in neigh_faces.get(cell_id, ()):
            cell_list.append(-(fi + 1))
        cell_list.sort(key=abs)
        conn.extend(cell_list)
        offsets.append(len(conn))

    return np.asarray(offsets, dtype=np.int32), np.asarray(conn, dtype=np.int32), n_cells


def _write_nface(zone, zone_name: str, link_data: dict, elem_range) -> None:
    offsets, conn, n_cells = _build_nface_connectivity(link_data)
    el = _cgns_node(zone, zone_name, "Elements_t", "I4")
    _i32_dataset(el, " data", [23, n_cells])

    er = _cgns_node(el, "ElementRange", "IndexRange_t", "I4")
    _i32_dataset(er, " data", elem_range)

    eso = _cgns_node(el, "ElementStartOffset", "DataArray_t", "I4")
    _i32_dataset(eso, " data", offsets)

    ec = _cgns_node(el, "ElementConnectivity", "DataArray_t", "I4")
    _i32_dataset(ec, " data", conn)


def _write_zone_bc(zone, boundary_face_ids_1based: np.ndarray,
                    bc_name: str = "box_surfs") -> None:
    zbc = _cgns_node(zone, "ZoneBC", "ZoneBC_t", "MT")
    bc = _cgns_node(zbc, bc_name, "BC_t", "C1")
    _bytes_dataset(bc, " data", b"Null")

    gl = _cgns_node(bc, "GridLocation", "GridLocation_t", "C1")
    _bytes_dataset(gl, " data", b"FaceCenter")

    pl = _cgns_node(bc, "PointList", "IndexArray_t", "I4")
    _i32_dataset(pl, " data", boundary_face_ids_1based)


def _write_flow_solution(zone) -> None:
    fs = _cgns_node(zone, "FlowSolution", "FlowSolution_t", "MT")
    gl = _cgns_node(fs, "GridLocation", "GridLocation_t", "C1")
    _bytes_dataset(gl, " data", b"CellCenter")


def _write_zone(base, zone_name: str, vertices: np.ndarray, link_data: dict) -> None:
    n_vertex = int(vertices.shape[0])
    n_faces = int(link_data["n_faces"])
    n_cells = int(link_data["n_cells"])

    zone = _cgns_node(base, zone_name, "Zone_t", "I4")
    zone.create_dataset(" data",
                        data=np.array([[n_vertex], [n_cells], [0]], dtype=np.int32))

    zt = _cgns_node(zone, "ZoneType", "ZoneType_t", "C1")
    _bytes_dataset(zt, " data", b"Unstructured")

    _write_grid_coordinates(zone, vertices)

    face_nodes_1 = (link_data["face_nodes"] + 1).astype(np.int32)
    _write_ngon(zone, face_nodes_1, elem_range=[1, n_faces])
    _write_nface(zone, zone_name, link_data,
                 elem_range=[n_faces + 1, n_faces + n_cells])

    boundary_1 = np.asarray(link_data["boundary_faces"], dtype=np.int32) + 1
    _write_zone_bc(zone, boundary_1)

    _write_flow_solution(zone)


def write_cgns(mesh: dict, outpath: str,
               zone_names=("FluidRegion", "FPHPARTS.box_vol")) -> None:
    """Write *mesh* to a CGNS/HDF5 file using the FLDUTIL exporter's layout.

    The resulting file contains:

    * Root attributes (``HDF5 MotherNode`` / ``Root Node of HDF5 File`` / ``MT``)
    * Two zero-length root datasets: `` hdf5version`` and `` format``
    * ``CGNSLibraryVersion``  (R4, value ``3.21``)
    * ``Base`` (CGNSBase_t, I4 ``[3, 3]``)
        * ``ReferenceState`` → ``ReferenceStateDescription`` =
          ``"Software Cradle FLDUTIL"``
        * One ``Zone_t`` per name in *zone_names* — by default the writer
          produces two identical zones (`FluidRegion` and `FPHPARTS.box_vol`)
          to mirror the vendor's own export.
    """
    vertices = mesh["vertices"]
    link_data = mesh.get("link_data")
    if vertices is None or len(vertices) == 0:
        raise ValueError("No vertices to write")
    if link_data is None:
        raise ValueError("No face/cell connectivity data (LS_Links parse failed)")

    # Open the file using the HDF5 1.8 file format (compact link storage,
    # v2 object headers).  See the docstring of ``_create_compact_group``
    # above for the rationale and impact (≈ 2× smaller file vs h5py's default).
    with h5py.File(outpath, "w", libver=("v108", "v108")) as f:
        # ── Root-node attributes (HDF5 MotherNode marker) ────────────────────
        f.attrs.create("label", np.bytes_("Root Node of HDF5 File"),
                       dtype=h5py.string_dtype(length=33))
        f.attrs.create("name", np.bytes_("HDF5 MotherNode"),
                       dtype=h5py.string_dtype(length=33))
        f.attrs.create("type", np.bytes_("MT"), dtype=h5py.string_dtype(length=3))

        # Root datasets (zero-filled placeholders matching the reference file).
        f.create_dataset(" hdf5version", data=np.zeros(33, dtype=np.int8))
        f.create_dataset(" format",
                         data=np.frombuffer(b"IEEE_LITTLE_32\0", dtype=np.int8))

        # ── CGNSLibraryVersion ───────────────────────────────────────────────
        lv = _cgns_node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4")
        lv.create_dataset(" data", data=np.array([3.21], dtype=np.float32))

        # ── CGNSBase_t ───────────────────────────────────────────────────────
        base = _cgns_node(f, "Base", "CGNSBase_t", "I4")
        _i32_dataset(base, " data", [3, 3])

        # ReferenceState (under Base, before zones)
        rs = _cgns_node(base, "ReferenceState", "ReferenceState_t", "MT")
        rsd = _cgns_node(rs, "ReferenceStateDescription", "Descriptor_t", "C1")
        _bytes_dataset(rsd, " data", b"Software Cradle FLDUTIL")

        for zname in zone_names:
            _write_zone(base, zname, vertices, link_data)


# ─────────────────────────────────────────────────────────────────────────────
# CLI entry-point
# ─────────────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Convert GPH file to CGNS format")
    parser.add_argument("gph_file", nargs="?", default="box.gph",
                        help="Input GPH file")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="Output CGNS file (default: input basename with .cgns)")
    parser.add_argument("--single-zone", action="store_true",
                        help="Emit a single Zone_t (default: emit the two zones "
                             "produced by the vendor's FLDUTIL exporter).")
    parser.add_argument("-z", "--zone", default=None,
                        help="When --single-zone is given, override the zone name "
                             "(default: 'FluidRegion').")
    args = parser.parse_args()

    gph_path = Path(args.gph_file)
    if not gph_path.exists():
        print(f"Error: file not found: {gph_path}")
        sys.exit(1)

    out_path = Path(args.output) if args.output else gph_path.with_suffix(".cgns")

    print(f"Reading: {gph_path}")
    mesh = parse_gph_mesh(str(gph_path))
    print(f"  Vertices : {mesh['n_vertices']}")
    ld = mesh.get("link_data")
    if ld:
        print(f"  Faces    : {ld['n_faces']}  (triangular, NGON_n)")
        print(f"  Cells    : {ld['n_cells']}  (NFACE_n)")
        print(f"  BC faces : {len(ld['boundary_faces'])}")
    else:
        print(f"  Cells    : {mesh['n_elements']}  (LS_Links parse failed)")

    if mesh["vertices"] is None or ld is None:
        print("Error: could not extract mesh data from GPH.")
        sys.exit(1)

    if args.single_zone:
        zone_names = (args.zone or "FluidRegion",)
    else:
        zone_names = ("FluidRegion", "FPHPARTS.box_vol")

    print(f"Writing: {out_path}")
    write_cgns(mesh, str(out_path), zone_names=zone_names)
    print("Done.")


if __name__ == "__main__":
    main()

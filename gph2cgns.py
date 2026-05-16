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
        "LS_SolverUnusedRegions", "LS_VolumeRegions", "LS_Parts",
        "LS_Assemblies",
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
    npe_p, _ = triples[2]

    BNDRY_U = 0xFFFFFFFF
    owner = np.frombuffer(data[owner_p : owner_p + n_faces_block_size],
                          dtype=">u4").astype(np.int64).copy()
    neigh_raw = np.frombuffer(data[neigh_p : neigh_p + n_faces_block_size],
                              dtype=">u4").copy()
    neigh = neigh_raw.astype(np.int64)
    neigh[neigh_raw == BNDRY_U] = -1

    # ``npe`` (nodes-per-face) drives variable-length face connectivity for
    # general polyhedral meshes.  For pure-triangle meshes every entry is 3
    # (e.g. ``box_ansa.gph``), but tr03.gph mixes faces with 3 to 11 nodes.
    npe = np.frombuffer(data[npe_p : npe_p + n_faces_block_size],
                        dtype=">u4").astype(np.int64).copy()
    conn_total_expected = int(npe.sum())

    # The connectivity block: the one whose byte_count matches
    # ``sum(npe) × 4`` exactly.  Falling back to the largest non-tri block
    # remains useful for files where ``npe`` is uniform (== 3) and the conn
    # block size = ``3 × n_faces × 4``.
    conn_block = None
    for p, bc in blocks:
        if (p, bc) in triples:
            continue
        if bc % 4 != 0:
            continue
        if bc // 4 == conn_total_expected:
            conn_block = (p, bc)
            break
    if conn_block is None:
        # Fall back: pick the largest remaining block (handles legacy files
        # whose conn byte-count includes trailing sentinel padding).
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
    if conn_total < conn_total_expected:
        return None

    conn = np.frombuffer(data[conn_p : conn_p + conn_total_expected * 4],
                         dtype=">u4").astype(np.int64).copy()

    # ``face_nodes`` is now stored in *CSR-style* form: a flat 1-D array of
    # 0-based vertex indices plus a length-``n_faces+1`` offset array.  Face
    # ``i`` is ``conn[face_offsets[i] : face_offsets[i+1]]``.  This handles
    # both pure-triangle meshes (npe ≡ 3) and arbitrary polyhedral meshes.
    face_offsets = np.empty(n_faces + 1, dtype=np.int64)
    face_offsets[0] = 0
    np.cumsum(npe, out=face_offsets[1:])
    face_nodes = conn  # flat 0-based vertex indices

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
        "npe": npe,                 # (n_faces,) int64 — nodes-per-face
        "face_nodes": face_nodes,   # flat (sum_npe,) int64, 0-based vertex IDs
        "face_offsets": face_offsets,  # (n_faces+1,) int64, cumulative-sum prefix
        "owner": owner,
        "neighbor": neigh,
        "boundary_faces": boundary_faces,
        "cell_owner_faces": cell_owner_faces,
        "cell_neighbor_faces": cell_neighbor_faces,
        "n_cells": n_cells,
    }


# ─────────────────────────────────────────────────────────────────────────────
# LS_CvolIdOfElements / LS_VolumeRegions / LS_Parts / LS_Assemblies
# ─────────────────────────────────────────────────────────────────────────────
#
# These four sections together describe how the mesh is partitioned into
# CGNS zones by the FLDUTIL exporter.  The vendor's reference file
# `tr03_orig.cgns` for example contains 7 zones (5 volume regions +
# 2 parts) that all share the same underlying mesh.
# ─────────────────────────────────────────────────────────────────────────────


def _parse_ls_cvol_ids(data: bytes) -> Optional[np.ndarray]:
    """Parse ``LS_CvolIdOfElements`` → 1-D int64 array, one entry per cell.

    Each entry holds the 1-based **Part index** of the cell.  For instance in
    ``tr03.gph`` this array has 63 882 entries with values 1 (Case part) or
    2 (Rotate part).  Returns ``None`` if the section is missing or malformed.
    """
    sec_start = _find_section(data, "LS_CvolIdOfElements")
    if sec_start < 0:
        return None
    sec_end = _section_end(data, sec_start)
    for p, bc in _iter_data_blocks(data, sec_start, sec_end):
        if bc % 4 == 0 and bc >= 4:
            return np.frombuffer(data[p : p + bc], dtype=">i4").astype(np.int64).copy()
    return None


def _parse_ls_string_section(data: bytes, section_name: str) -> list[str]:
    """Return the list of ASCII strings stored in a section like
    ``LS_VolumeRegions`` or ``LS_Parts``."""
    sec_start = _find_section(data, section_name)
    if sec_start < 0:
        return []
    sec_end = _section_end(data, sec_start)
    out: list[str] = []
    for p, bc in _iter_data_blocks(data, sec_start, sec_end):
        raw = data[p : p + bc]
        # ASCII-only payload, NUL/space-padded
        if all(b == 0 or 32 <= b < 127 for b in raw):
            s = raw.decode("ascii", errors="replace").strip("\x00").rstrip()
            if s:
                out.append(s)
    return out


def _parse_ls_surface_regions(data: bytes) -> list[tuple[str, np.ndarray]]:
    """Parse ``LS_SurfaceRegions`` into a list of ``(name, gph_face_ids)``.

    Each surface region in the GPH file consists of three consecutive data
    blocks: ``name`` (ASCII, NUL/space-padded), ``face_ids`` (I4 array, one
    entry per face in the region) and ``weights`` (I4 array of the same
    length, currently always 1).  The returned face IDs are **0-based**
    indices into the global ``LS_Links`` face array — i.e. the natural
    indexing used internally by scFLOW / FLDUTIL.

    These regions become CGNS ``BC_t`` family groups under each zone's
    ``ZoneBC`` (one BC per region; faces filtered against the zone's
    sub-mesh).  For ``tr03.gph`` this yields 15 BCs per zone:

      inlet, outlet, Rotate_Plane, Rotate_Cylinder, Rotate_Plane[2],
      @PartSurface_Case, @PartSurface_Rotate, @PartSurface_Impeller,
      Rotate_MovingFaceRegion,
      Rotate_Plane_Moving, Rotate_Plane_Static,
      Rotate_Cylinder_Moving, Rotate_Cylinder_Static,
      Rotate_Plane[2]_Moving, Rotate_Plane[2]_Static
    """
    sec_start = _find_section(data, "LS_SurfaceRegions")
    if sec_start < 0:
        return []
    sec_end = _section_end(data, sec_start)
    blocks = list(_iter_data_blocks(data, sec_start, sec_end))

    regions: list[tuple[str, np.ndarray]] = []
    i = 0
    while i + 2 < len(blocks):
        p_n, bc_n = blocks[i]
        p_i, bc_i = blocks[i + 1]
        p_w, bc_w = blocks[i + 2]
        name_raw = data[p_n : p_n + bc_n]
        # Region name block: NUL/space-padded ASCII.
        if not all(b == 0 or 32 <= b < 127 for b in name_raw):
            i += 1
            continue
        name = name_raw.decode("ascii", errors="replace").strip("\x00").rstrip()
        if not name:
            i += 1
            continue
        # Face-IDs block and weights block must have the same shape.
        if bc_i > 0 and bc_i == bc_w and bc_i % 4 == 0:
            face_ids = np.frombuffer(data[p_i : p_i + bc_i],
                                     dtype=">i4").astype(np.int64).copy()
            regions.append((name, face_ids))
            i += 3
        else:
            i += 1
    return regions


def _parse_ls_assemblies(data: bytes) -> dict[str, Optional[str]]:
    """Parse ``LS_Assemblies`` XML → mapping ``part_name → assembly_name``.

    Parts that live at the XML root (not nested inside an ``<assembly>``)
    map to ``None`` and are emitted by the writer as ``FPHPARTS.<part>``;
    parts inside ``<assembly name="A">`` are emitted as
    ``FPHPARTS.A.<part>``.
    """
    sec_start = _find_section(data, "LS_Assemblies")
    if sec_start < 0:
        return {}
    sec_end = _section_end(data, sec_start)
    xml_bytes = b""
    for p, bc in _iter_data_blocks(data, sec_start, sec_end):
        chunk = data[p : p + bc]
        if chunk.lstrip().startswith(b"<?xml") or b"<part" in chunk:
            xml_bytes = chunk
            break
    if not xml_bytes:
        return {}
    try:
        import xml.etree.ElementTree as ET
        root = ET.fromstring(xml_bytes.decode("utf-8", errors="replace"))
    except Exception:
        return {}
    mapping: dict[str, Optional[str]] = {}
    # Parts inside <assembly> nodes
    for assembly in root.iter("assembly"):
        a_name = assembly.get("name")
        for part in assembly.findall("part"):
            p_name = part.get("name")
            if p_name:
                mapping[p_name] = a_name
    # Parts at the XML root that are NOT inside any assembly
    for part in root.findall("part"):
        p_name = part.get("name")
        if p_name and p_name not in mapping:
            mapping[p_name] = None
    return mapping


def _classify_zone_cells(zone_name: str, part_names: list[str],
                          cvol_id: Optional[np.ndarray], n_cells: int) -> np.ndarray:
    """Return a boolean cell-mask selecting cells that belong to *zone_name*.

    Cell membership for each CGNS zone is **not stored explicitly** in the
    GPH file — it must be derived from the ``LS_CvolIdOfElements`` part
    assignment and the zone name (which the upstream ANSA / scFLOW workflow
    encodes via naming conventions):

    * ``FluidRegion`` — convention for the *whole* mesh.
    * ``@VPartRegion_<P>[<N>]`` — exactly the cells of Part ``P``.
    * ``FPHPARTS.[<assembly>.]<P>`` — exactly the cells of Part ``P``.
    * Other names — substring-match a Part name; falls back to the whole
      mesh if no part is mentioned.
    """
    all_mask = np.ones(n_cells, dtype=bool)
    if zone_name == "FluidRegion":
        return all_mask
    if cvol_id is None or len(cvol_id) != n_cells or not part_names:
        return all_mask

    # Explicit `FPHPARTS.[<asm>.]<part>` zone name
    if zone_name.startswith("FPHPARTS."):
        suffix = zone_name[len("FPHPARTS.") :]
        candidate = suffix.rsplit(".", 1)[-1]
        if candidate in part_names:
            return cvol_id == (part_names.index(candidate) + 1)

    # `@VPartRegion_<part>[<n>]` zone name → strip prefix and trailing `[n]`
    if zone_name.startswith("@VPartRegion_"):
        rem = zone_name[len("@VPartRegion_") :]
        rem = rem.split("[", 1)[0]
        if rem in part_names:
            return cvol_id == (part_names.index(rem) + 1)

    # Substring fallback: pick the longest Part name that appears inside the
    # zone name (longer names are more specific).
    matches = sorted(
        (p for p in part_names if p and p in zone_name),
        key=len, reverse=True,
    )
    if matches:
        return cvol_id == (part_names.index(matches[0]) + 1)
    return all_mask


# ─────────────────────────────────────────────────────────────────────────────
# Vertex renumbering — match the official FLDUTIL exporter's ordering
# ─────────────────────────────────────────────────────────────────────────────


def _renumber_by_first_use(face_nodes_flat: np.ndarray, n_vertices: int) -> np.ndarray:
    """Return a permutation ``perm`` such that ``perm[old_gph_id] = new_id``
    where ``new_id`` is the 0-based first-use index of ``old_gph_id`` in the
    scan of the flat (CSR-style) face-node connectivity array.

    Vertices that never appear in any face are appended at the end in their
    original order so the permutation is a bijection.
    """
    perm = np.full(n_vertices, -1, dtype=np.int64)
    next_id = 0
    # Vectorised first-pass: collect the index of the first occurrence of each
    # vertex via numpy.unique with return_index=True.
    flat = np.asarray(face_nodes_flat).reshape(-1)
    flat_valid = flat[(flat >= 0) & (flat < n_vertices)]
    if flat_valid.size:
        # np.unique returns sorted unique values and the index of the first
        # occurrence of each.  We re-sort by index to recover scan order.
        uniq, first_idx = np.unique(flat_valid, return_index=True)
        order = np.argsort(first_idx)
        unique_in_scan_order = uniq[order]
        perm[unique_in_scan_order] = np.arange(len(unique_in_scan_order),
                                                dtype=np.int64)
        next_id = len(unique_in_scan_order)
    # Append any vertices that were never referenced by a face.
    missing = np.where(perm == -1)[0]
    perm[missing] = np.arange(next_id, next_id + missing.size, dtype=np.int64)
    return perm


# ─────────────────────────────────────────────────────────────────────────────
# GPH front-end
# ─────────────────────────────────────────────────────────────────────────────


def parse_gph_mesh(filepath: str) -> dict:
    """Extract mesh data (vertices, faces, cells) from a GPH file.

    The returned vertices/face_nodes have already been renumbered so that the
    output CGNS matches the vendor's own conversion (vertices appear in the
    order they are first referenced by the face connectivity).

    The result also carries the zone-partition metadata needed to write a
    multi-Zone CGNS that mirrors the vendor's FLDUTIL exporter
    (``LS_CvolIdOfElements``, ``LS_VolumeRegions``, ``LS_Parts``,
    ``LS_Assemblies`` XML).
    """
    with open(filepath, "rb") as f:
        data = f.read()

    result: dict = {
        "vertices": None,
        "n_vertices": 0,
        "link_data": None,
        "n_elements": 0,
        "cvol_id": None,
        "volume_regions": [],
        "parts": [],
        "part_assembly": {},
        "surface_regions": [],
    }

    xyz, n_vertices = _parse_ls_nodes(data)
    link_data = _parse_ls_links(data)

    # Partition metadata is always parsed if present, even if mesh is broken.
    result["cvol_id"] = _parse_ls_cvol_ids(data)
    result["volume_regions"] = _parse_ls_string_section(data, "LS_VolumeRegions")
    result["parts"] = _parse_ls_string_section(data, "LS_Parts")
    result["part_assembly"] = _parse_ls_assemblies(data)
    result["surface_regions"] = _parse_ls_surface_regions(data)

    if xyz is None or n_vertices == 0 or link_data is None:
        return result

    # Clamp any garbage node indices that exceed the known vertex count
    # (occasional byte-count sentinel leakage seen on legacy files).  The
    # face-node array is now stored as a flat 1-D vector (CSR-style); index
    # bound-checking is applied directly to it.
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
    link_data["face_nodes"] = perm[link_data["face_nodes"]]

    result["vertices"] = xyz_renum
    result["n_vertices"] = n_vertices
    result["link_data"] = link_data
    result["n_elements"] = link_data["n_cells"]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# CGNS / HDF5 writer
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# HDF5 file-format selection
# ─────────────────────────────────────────────────────────────────────────────
#
# CGNS files exported by Software Cradle's ``FLDUTIL`` tool (see the
# reference file ``box_ansa_orig.cgns``) use the **v0 superblock** with
# v1 object headers carrying HDF5 1.8 "Link Info" / "Group Info" messages.
# Some third-party CGNS readers — including ANSA — only accept v0/v1
# superblocks and reject the newer v2 superblock (HDF5 1.8 file format)
# with a "No bases found!" error.
#
# We therefore open the file with ``libver=("earliest", "v108")`` (h5py's
# default), which writes a **v0 superblock** with the legacy v1.6 group
# storage (symbol-table + B-tree + local-heap triple).  This is universally
# readable by every CGNS toolkit ever released, at the cost of roughly
# 1 KB of metadata per group (the resulting box_ansa.cgns is ~73 KB vs the
# 35 KB FLDUTIL reference; the *content* of every group/dataset/attribute
# remains a byte-perfect PERFECT MATCH against the reference).
#
# The smaller v2-superblock layout that the FLDUTIL exporter achieves
# (35 KB) requires HDF5's "new" group format with v0 superblock and Link
# Info messages, which is controlled by ``H5Pset_link_phase_change`` *plus*
# internal flags that the h5py Python binding does not expose.  Falling
# back to v2-only superblock via ``libver=("v108", "v108")`` shrinks the
# file to 31 KB but produces a v2-superblock file that ANSA refuses to
# read — hence the choice below.
# ─────────────────────────────────────────────────────────────────────────────


_MAX_CHILDREN_PER_GROUP = 7  # Zone_t children: data, ZoneType, GridCoordinates,
                              # GridElements_Faces, <NFace>, ZoneBC, FlowSolution


def _create_compact_group(parent, name: str):
    """Create an HDF5 group under *parent* with ANSA-compatible storage.

    ``track_order=False`` suppresses the HDF5 1.8 attribute-creation-order
    index (the FLDUTIL reference does not write it either).  Combined with
    the v0-superblock file format selected in :func:`write_cgns`, the
    resulting layout is read correctly by every CGNS toolkit including
    ANSA and the official CGNS C/C++ libraries.
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


def _write_ngon(zone, face_nodes_1based_flat: np.ndarray,
                 face_offsets: np.ndarray, n_faces: int, elem_range) -> None:
    """Write CGNS NGON_n element section using CSR-style face connectivity.

    ``face_nodes_1based_flat``: flat 1-D int array of 1-based vertex IDs
    concatenated across all faces (length = sum of face_nodes-per-face).
    ``face_offsets``: length ``n_faces + 1`` cumulative-sum prefix
    ``[0, n0, n0+n1, ...]`` indexing into ``face_nodes_1based_flat``.

    This single code path handles both pure-triangle meshes (every face
    contributes 3 entries) and arbitrary polyhedral meshes (faces with
    varying node counts, e.g. tr03.gph has faces with 3..11 nodes).
    """
    el = _cgns_node(zone, "GridElements_Faces", "Elements_t", "I4")
    _i32_dataset(el, " data", [22, n_faces])

    er = _cgns_node(el, "ElementRange", "IndexRange_t", "I4")
    _i32_dataset(er, " data", elem_range)

    eso = _cgns_node(el, "ElementStartOffset", "DataArray_t", "I4")
    _i32_dataset(eso, " data", face_offsets)

    ec = _cgns_node(el, "ElementConnectivity", "DataArray_t", "I4")
    _i32_dataset(ec, " data", face_nodes_1based_flat)


def _extract_zone_submesh(vertices: np.ndarray, link_data: dict,
                            cell_mask: np.ndarray) -> dict:
    """Return a self-contained sub-mesh describing only the cells selected
    by *cell_mask*.

    Each output CGNS Zone is a stand-alone unstructured mesh with its own
    local vertex / face / cell numbering (1-based when written to CGNS).
    This helper performs the renumbering and produces an object that
    structurally matches the ``link_data`` shape so the existing CGNS-writer
    code paths can be reused unchanged.

    Algorithm
    ---------
    1. *Faces in the sub-mesh* are those whose owner OR neighbour belongs
       to the selected cell set.
    2. *Vertices in the sub-mesh* are the union of the face-node lists of
       the kept faces.
    3. Face / vertex / cell IDs are renumbered to a dense 0-based range
       (then converted to 1-based when written to CGNS).
    4. Faces that originally had ``neighbour == -1`` (boundary of the full
       mesh) remain boundary faces of the sub-mesh.  Faces whose
       owner-or-neighbour pair straddles the zone boundary become
       *new boundary faces* of the sub-mesh; we still keep them with their
       valid (in-zone) side as the owner so the sub-mesh has a closed
       surface.
    """
    n_cells_full = link_data["n_cells"]
    n_faces_full = link_data["n_faces"]
    owner_full = link_data["owner"]
    neigh_full = link_data["neighbor"]
    face_nodes_full = link_data["face_nodes"]
    face_offsets_full = link_data["face_offsets"]
    npe_full = link_data["npe"]

    cell_mask = cell_mask.astype(bool)

    # ── 1. Pick faces incident on the zone ────────────────────────────────
    owner_in = np.zeros(n_faces_full, dtype=bool)
    neigh_in = np.zeros(n_faces_full, dtype=bool)
    ow_arr = owner_full
    valid_ow = (ow_arr >= 0) & (ow_arr < n_cells_full)
    owner_in[valid_ow] = cell_mask[ow_arr[valid_ow]]
    nb_arr = neigh_full
    valid_nb = (nb_arr >= 0) & (nb_arr < n_cells_full)
    neigh_in[valid_nb] = cell_mask[nb_arr[valid_nb]]

    face_mask = owner_in | neigh_in
    kept_face_ids = np.where(face_mask)[0]
    n_faces_sub = int(kept_face_ids.size)

    # ── 2. Renumber face owners / neighbours so the kept (in-zone) cell
    #     is always the owner.  Track which side originally had it. ───────
    new_owner = np.full(n_faces_sub, -1, dtype=np.int64)
    new_neigh = np.full(n_faces_sub, -1, dtype=np.int64)
    flip = np.zeros(n_faces_sub, dtype=bool)

    ow_kept = ow_arr[kept_face_ids]
    nb_kept = nb_arr[kept_face_ids]
    own_inzone = owner_in[kept_face_ids]
    nb_inzone = neigh_in[kept_face_ids]
    both_inzone = own_inzone & nb_inzone

    new_owner[own_inzone] = ow_kept[own_inzone]
    new_neigh[both_inzone] = nb_kept[both_inzone]
    only_nb = nb_inzone & ~own_inzone
    new_owner[only_nb] = nb_kept[only_nb]
    flip[only_nb] = True

    # ── 3. Re-index cells (kept cells → 0..n_cells_sub-1) ─────────────────
    kept_cell_ids = np.where(cell_mask)[0]
    cell_remap = np.full(n_cells_full, -1, dtype=np.int64)
    cell_remap[kept_cell_ids] = np.arange(kept_cell_ids.size, dtype=np.int64)
    new_owner[new_owner >= 0] = cell_remap[new_owner[new_owner >= 0]]
    new_neigh_in_zone_mask = new_neigh >= 0
    new_neigh[new_neigh_in_zone_mask] = cell_remap[new_neigh[new_neigh_in_zone_mask]]
    n_cells_sub = int(kept_cell_ids.size)

    # ── 4. Slice face_nodes for kept faces ────────────────────────────────
    sub_npe = npe_full[kept_face_ids]
    sub_face_offsets = np.empty(n_faces_sub + 1, dtype=np.int64)
    sub_face_offsets[0] = 0
    np.cumsum(sub_npe, out=sub_face_offsets[1:])

    sub_total_conn = int(sub_face_offsets[-1])
    sub_face_nodes = np.empty(sub_total_conn, dtype=np.int64)
    # Per-face copy: vectorised gather via per-face ranges.
    # For very large sub-meshes this is the hot loop; numpy fancy indexing
    # over a single concatenated index array is roughly 10× faster than
    # a Python for-loop.
    write_pos = sub_face_offsets[:-1]
    # Build a single source-index array by concatenating per-face slices.
    src_indices = np.concatenate([
        np.arange(face_offsets_full[fi], face_offsets_full[fi + 1],
                  dtype=np.int64)
        for fi in kept_face_ids
    ]) if n_faces_sub else np.empty(0, dtype=np.int64)
    sub_face_nodes[:] = face_nodes_full[src_indices]

    # If a face's owner/neighbour roles were swapped, reverse its node order
    # so the CGNS face-normal convention is preserved (the surface still
    # points *outward* from the in-zone cell).
    if flip.any():
        for li in np.where(flip)[0]:
            s = int(sub_face_offsets[li])
            e = int(sub_face_offsets[li + 1])
            sub_face_nodes[s:e] = sub_face_nodes[s:e][::-1]

    # ── 5. Re-index vertices (kept vertices → 0..n_verts_sub-1) ──────────
    used_v = np.unique(sub_face_nodes)
    n_verts_sub = int(used_v.size)
    v_remap = np.full(vertices.shape[0], -1, dtype=np.int64)
    v_remap[used_v] = np.arange(n_verts_sub, dtype=np.int64)
    sub_face_nodes = v_remap[sub_face_nodes]
    sub_vertices = vertices[used_v]

    # ── 6. Build per-cell face lists (owner/neighbour) ────────────────────
    cell_owner_faces: dict = defaultdict(list)
    cell_neighbor_faces: dict = defaultdict(list)
    for li in range(n_faces_sub):
        ow = int(new_owner[li])
        if 0 <= ow < n_cells_sub:
            cell_owner_faces[ow].append(li)
        nb = int(new_neigh[li])
        if 0 <= nb < n_cells_sub:
            cell_neighbor_faces[nb].append(li)

    boundary_faces = list(np.where(new_neigh == -1)[0])

    # GPH-face-id → local-face-id (1-based) map.  ``-1`` means the GPH face
    # is not present in this sub-mesh.  Surface-region BC writers use this
    # to remap LS_SurfaceRegions face IDs to per-zone CGNS face IDs.
    gph_to_local_1based = np.full(n_faces_full, 0, dtype=np.int64)
    gph_to_local_1based[kept_face_ids] = np.arange(1, n_faces_sub + 1,
                                                    dtype=np.int64)

    sub_link = {
        "n_faces": n_faces_sub,
        "npe": sub_npe,
        "face_nodes": sub_face_nodes,
        "face_offsets": sub_face_offsets,
        "owner": new_owner,
        "neighbor": new_neigh,
        "boundary_faces": boundary_faces,
        "cell_owner_faces": cell_owner_faces,
        "cell_neighbor_faces": cell_neighbor_faces,
        "n_cells": n_cells_sub,
        "gph_to_local_face_1based": gph_to_local_1based,
    }
    return {"vertices": sub_vertices, "link_data": sub_link}


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


def _write_zone_bc(zone, bc_families: list[tuple[str, np.ndarray]]) -> None:
    """Write the ``ZoneBC_t`` node with one ``BC_t`` family per entry of
    *bc_families*.

    Each family is a ``(name, face_ids_1based_int32)`` pair.  Empty face
    lists are still emitted as a BC group (shape ``(0,)`` PointList) so
    every zone has the same set of named families — this matches the
    layout of ``tr03_orig.cgns`` where e.g. ``inlet`` appears in every
    zone, but its PointList is empty in sub-mesh zones that don't include
    the inlet boundary.
    """
    zbc = _cgns_node(zone, "ZoneBC", "ZoneBC_t", "MT")
    for bc_name, face_ids in bc_families:
        bc = _cgns_node(zbc, bc_name, "BC_t", "C1")
        _bytes_dataset(bc, " data", b"Null")

        gl = _cgns_node(bc, "GridLocation", "GridLocation_t", "C1")
        _bytes_dataset(gl, " data", b"FaceCenter")

        pl = _cgns_node(bc, "PointList", "IndexArray_t", "I4")
        # When the zone contains *no* faces from this BC region the vendor
        # exporter creates the PointList group with the standard CGNS
        # attributes but *omits* the ' data' dataset entirely.  Match that
        # so e.g. the reference's empty ``inlet`` in Rotate_MovingVolumeRegion
        # is byte-equivalent to our output.
        if face_ids.size > 0:
            _i32_dataset(pl, " data", face_ids)


def _write_flow_solution(zone) -> None:
    fs = _cgns_node(zone, "FlowSolution", "FlowSolution_t", "MT")
    gl = _cgns_node(fs, "GridLocation", "GridLocation_t", "C1")
    _bytes_dataset(gl, " data", b"CellCenter")


def _write_zone(base, zone_name: str, vertices: np.ndarray, link_data: dict,
                 bc_families: Optional[list[tuple[str, np.ndarray]]] = None) -> None:
    n_vertex = int(vertices.shape[0])
    n_faces = int(link_data["n_faces"])
    n_cells = int(link_data["n_cells"])

    zone = _cgns_node(base, zone_name, "Zone_t", "I4")
    zone.create_dataset(" data",
                        data=np.array([[n_vertex], [n_cells], [0]], dtype=np.int32))

    zt = _cgns_node(zone, "ZoneType", "ZoneType_t", "C1")
    _bytes_dataset(zt, " data", b"Unstructured")

    _write_grid_coordinates(zone, vertices)

    face_nodes_1_flat = (link_data["face_nodes"] + 1).astype(np.int32)
    face_offsets_i4 = link_data["face_offsets"].astype(np.int32)
    _write_ngon(zone, face_nodes_1_flat, face_offsets_i4, n_faces,
                elem_range=[1, n_faces])
    _write_nface(zone, zone_name, link_data,
                 elem_range=[n_faces + 1, n_faces + n_cells])

    if bc_families is None:
        # Legacy single-BC layout used when no LS_SurfaceRegions are
        # available (e.g. synthetic test files): emit one BC family named
        # ``box_surfs`` listing every boundary face of the zone.
        boundary_1 = np.asarray(link_data["boundary_faces"],
                                 dtype=np.int32) + 1
        bc_families = [("box_surfs", boundary_1)]
    _write_zone_bc(zone, bc_families)

    _write_flow_solution(zone)


def _bc_families_for_zone(
    surface_regions: list[tuple[str, np.ndarray]],
    gph_to_local_1based: np.ndarray,
) -> list[tuple[str, np.ndarray]]:
    """Project the global GPH ``LS_SurfaceRegions`` definitions onto a single
    CGNS zone using *gph_to_local_1based* (length n_faces_full, value 0 for
    GPH faces absent from the zone, otherwise the zone-local 1-based face
    ID).

    Returns the list of ``(region_name, face_ids_1based_int32)`` pairs to
    hand to :func:`_write_zone_bc`.  Regions with no faces present in the
    zone are still emitted (with an empty PointList) so the BC structure
    is uniform across all zones — this matches the layout that the
    vendor's FLDUTIL exporter produces in ``tr03_orig.cgns``.
    """
    out: list[tuple[str, np.ndarray]] = []
    for name, gph_face_ids in surface_regions:
        # gph_face_ids are 0-based global GPH indices.  Look up each in the
        # zone-local-id table; entries equal to 0 mean "face not present in
        # this zone" and are filtered out.
        if gph_face_ids.size:
            local = gph_to_local_1based[gph_face_ids]
            local = local[local > 0].astype(np.int32)
        else:
            local = np.empty(0, dtype=np.int32)
        out.append((name, local))
    return out


def _build_zone_plan(mesh: dict,
                      override_zone_names: Optional[tuple] = None) -> list[tuple[str, np.ndarray]]:
    """Return the ordered list of ``(zone_name, cell_mask)`` to emit.

    Driven by the partition metadata parsed from ``LS_VolumeRegions``,
    ``LS_Parts`` and ``LS_Assemblies``.  Falls back to a single
    ``FluidRegion`` + ``FPHPARTS.<part>`` pair when no partition metadata
    is present (matches the legacy ``box_ansa.gph`` behaviour).
    """
    link_data = mesh["link_data"]
    n_cells = int(link_data["n_cells"])
    cvol_id = mesh.get("cvol_id")
    parts = mesh.get("parts", [])
    regions = mesh.get("volume_regions", [])
    assembly_map = mesh.get("part_assembly", {})

    # CLI override: caller supplied an explicit zone name list.
    if override_zone_names:
        all_mask = np.ones(n_cells, dtype=bool)
        return [(name, all_mask) for name in override_zone_names]

    # Build the ordered zone list.  Volume regions come first (in file
    # order), followed by parts.  This matches the vendor's output:
    #
    #   tr03_orig.cgns Base/ children order =
    #     Rotate_MovingVolumeRegion, FluidRegion, Rotate_Moving,
    #     @VPartRegion_Case[2], @VPartRegion_Rotate[2],
    #     FPHPARTS.tr03.Case, FPHPARTS.Rotate
    plan: list[tuple[str, np.ndarray]] = []

    if not regions and not parts:
        # No partition metadata — fall back to the legacy 2-zone layout.
        all_mask = np.ones(n_cells, dtype=bool)
        plan.append(("FluidRegion", all_mask))
        plan.append(("FPHPARTS.box_vol", all_mask))
        return plan

    # Volume regions (use the names as stored in the GPH file).
    for region in regions:
        mask = _classify_zone_cells(region, parts, cvol_id, n_cells)
        plan.append((region, mask))

    # If LS_VolumeRegions does NOT include "FluidRegion" but a full-mesh
    # zone is conventional, prepend it so downstream CFD tools always have
    # an "all cells" view.  In tr03 this entry already exists; in box_ansa
    # it does too — only synthetic / hand-crafted files would need this.
    if not any(name == "FluidRegion" for name, _ in plan) and len(parts) > 0:
        plan.insert(0, ("FluidRegion", np.ones(n_cells, dtype=bool)))

    # Parts (each becomes one FPHPARTS.* zone).
    for idx, part in enumerate(parts, start=1):
        a_name = assembly_map.get(part)
        zone_name = (f"FPHPARTS.{a_name}.{part}" if a_name
                     else f"FPHPARTS.{part}")
        mask = (cvol_id == idx) if (cvol_id is not None and len(cvol_id) == n_cells) \
                                 else np.ones(n_cells, dtype=bool)
        plan.append((zone_name, mask))

    return plan


def write_cgns(mesh: dict, outpath: str,
               zone_names: Optional[tuple] = None) -> None:
    """Write *mesh* to a CGNS/HDF5 file using the FLDUTIL exporter's layout.

    The resulting file contains:

    * Root attributes (``HDF5 MotherNode`` / ``Root Node of HDF5 File`` / ``MT``)
    * Two zero-length root datasets: `` hdf5version`` and `` format``
    * ``CGNSLibraryVersion``  (R4, value ``3.21``)
    * ``Base`` (CGNSBase_t, I4 ``[3, 3]``)
        * ``ReferenceState`` → ``ReferenceStateDescription`` =
          ``"Software Cradle FLDUTIL"``
        * One ``Zone_t`` per partition discovered in the GPH metadata
          (``LS_VolumeRegions`` + ``LS_Parts`` + ``LS_Assemblies``).  Falls
          back to the legacy ``FluidRegion`` + ``FPHPARTS.box_vol`` pair
          when no partition metadata is present.  Passing *zone_names*
          forces a specific list of full-mesh zones (legacy compatibility).
    """
    vertices = mesh["vertices"]
    link_data = mesh.get("link_data")
    if vertices is None or len(vertices) == 0:
        raise ValueError("No vertices to write")
    if link_data is None:
        raise ValueError("No face/cell connectivity data (LS_Links parse failed)")

    zone_plan = _build_zone_plan(mesh, override_zone_names=zone_names)

    # Use ``libver=("earliest", "v108")`` (h5py's default) to write a
    # **v0 superblock** file.  Some CGNS readers — notably ANSA — refuse
    # to open v2-superblock files (HDF5 1.8 file format) with a
    # "No bases found!" error, so we deliberately stay on the older,
    # universally-readable file format here.  See the module-level
    # comment above ``_create_compact_group`` for full rationale.
    with h5py.File(outpath, "w", libver=("earliest", "v108")) as f:
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

        surface_regions = mesh.get("surface_regions", [])
        n_faces_full = int(link_data["n_faces"])

        # For the FluidRegion / full-mesh path, the GPH face IDs equal the
        # zone-local CGNS face IDs (modulo the 0/1-based shift), so the
        # identity map is the natural choice.
        full_gph_to_local = np.concatenate([[0],
                                             np.arange(1, n_faces_full + 1,
                                                       dtype=np.int64)])
        # ``gph_face_0based → gph_to_local_1based[gph_face_0based]`` would
        # give a 0-prefixed lookup; we instead use it as a direct index by
        # storing ``local = gph + 1`` at position ``gph``.  Easier:
        full_gph_to_local = np.arange(1, n_faces_full + 1, dtype=np.int64)

        for zname, cell_mask in zone_plan:
            if cell_mask.all():
                # Whole-mesh zone — no sub-mesh extraction needed.
                bc_families = (_bc_families_for_zone(surface_regions,
                                                      full_gph_to_local)
                               if surface_regions else None)
                _write_zone(base, zname, vertices, link_data,
                            bc_families=bc_families)
            else:
                sub = _extract_zone_submesh(vertices, link_data, cell_mask)
                gph_to_local = sub["link_data"]["gph_to_local_face_1based"]
                bc_families = (_bc_families_for_zone(surface_regions,
                                                      gph_to_local)
                               if surface_regions else None)
                _write_zone(base, zname, sub["vertices"], sub["link_data"],
                            bc_families=bc_families)


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
        npe_arr = ld.get("npe")
        if npe_arr is not None and len(npe_arr) > 0:
            mn, mx = int(npe_arr.min()), int(npe_arr.max())
            face_kind = ("triangular" if mn == mx == 3
                         else f"polyhedral, {mn}..{mx} nodes per face")
        else:
            face_kind = "mixed"
        print(f"  Faces    : {ld['n_faces']}  ({face_kind}, NGON_n)")
        print(f"  Cells    : {ld['n_cells']}  (NFACE_n)")
        print(f"  BC faces : {len(ld['boundary_faces'])}")
    else:
        print(f"  Cells    : {mesh['n_elements']}  (LS_Links parse failed)")

    if mesh["vertices"] is None or ld is None:
        print("Error: could not extract mesh data from GPH.")
        sys.exit(1)

    # Show the zone plan that will be emitted (volume regions + parts as
    # parsed from LS_VolumeRegions / LS_Parts / LS_Assemblies XML).
    if args.single_zone:
        zone_names = (args.zone or "FluidRegion",)
    else:
        zone_names = None  # let write_cgns infer from the partition metadata
    plan = _build_zone_plan(mesh, override_zone_names=zone_names)
    print(f"  Zones    : {len(plan)}")
    for zname, mask in plan:
        cell_count = int(mask.sum())
        print(f"             - {zname}  ({cell_count} cells)")

    print(f"Writing: {out_path}")
    write_cgns(mesh, str(out_path), zone_names=zone_names)
    print("Done.")


if __name__ == "__main__":
    main()

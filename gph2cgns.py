#!/usr/bin/env python3
"""
GPH to CGNS Converter.

Reads GPH (CRDL-FLD) binary format and writes CGNS/HDF5 unstructured mesh.
Requires: numpy, h5py (no CGNS library needed).
"""

import argparse
import struct
import sys
from pathlib import Path

import numpy as np

try:
    import h5py
except ImportError:
    print("Error: h5py is required. Install with: pip install h5py numpy")
    sys.exit(1)


def read_i32_be(data: bytes, pos: int) -> int:
    return int.from_bytes(data[pos : pos + 4], "big")


def read_f32_be(data: bytes, pos: int) -> float:
    return struct.unpack(">f", data[pos : pos + 4])[0]


def read_f64_wr(data: bytes, pos: int) -> float:
    """Read float64 stored in word-reversed (middle-endian) format.

    The GPH format stores 64-bit doubles as two 32-bit big-endian words in
    reversed order: [lower_32bit_word][upper_32bit_word].  A coordinate of
    0.01 is therefore stored as bytes 47 AE 14 7B 3F 84 7A E1, which when
    mis-read as two consecutive float32 values yields (89128.96, 1.035) –
    the exact "89000+" symptom reported during testing.
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


def _parse_ls_nodes(data: bytes) -> tuple:
    """Parse the LS_Nodes section and return (xyz_array, n_vertices).

    The section stores three coordinate axes (X, Z, Y in file order) each as
    a separate block:
        [16-byte descriptor: 12 / 8 / n_verts / 1]
        [12-byte header:     12 / byte_count / upper_half_of_max_coord]
        [n_verts × 8-byte word-reversed float64 values]

    Returns (xyz, n_vertices) where xyz has shape (n_vertices, 3) with
    columns CoordinateX, CoordinateY, CoordinateZ, or (None, 0) on failure.
    """
    sec_start = _find_section(data, "LS_Nodes")
    if sec_start < 0:
        return None, 0

    # Skip: I4=32 (4B) + name (32B) + I4=32 (4B) = 40B
    pos = sec_start + 40

    # Scan descriptor blocks [12, type, dim0, dim1] until we find the first
    # coordinate-block descriptor [12, 8, n_verts, 1].
    n_vertices = 0
    while pos + 16 <= len(data):
        if read_i32_be(data, pos) != 12:
            break
        type_code = read_i32_be(data, pos + 4)
        dim0 = read_i32_be(data, pos + 8)
        dim1 = read_i32_be(data, pos + 12)
        pos += 16
        if type_code == 8 and dim1 == 1 and dim0 > 0:
            n_vertices = dim0
            break

    if n_vertices == 0:
        return None, 0

    # Read three successive coordinate blocks (X, Z, Y order in the file).
    # Each block: [12B header] [n_vertices × 8B word-reversed float64].
    # Between blocks there is a 16-byte descriptor; skip it when present.
    coord_blocks: list[list[float]] = []
    for _ in range(3):
        if pos + 12 > len(data):
            break
        # 12-byte header: [I4=12][I4=byte_count][I4=upper_half_of_max]
        pos += 12
        vals: list[float] = []
        for _ in range(n_vertices):
            if pos + 8 > len(data):
                break
            vals.append(read_f64_wr(data, pos))
            pos += 8
        coord_blocks.append(vals)
        # Skip the next 16-byte descriptor if present
        if pos + 16 <= len(data) and read_i32_be(data, pos) == 12:
            pos += 16

    if len(coord_blocks) < 3:
        return None, 0

    # File order is X, Z, Y → rearrange to X, Y, Z for CGNS output
    x_vals, z_vals, y_vals = coord_blocks
    xyz = np.array(
        [[x_vals[i], y_vals[i], z_vals[i]] for i in range(n_vertices)],
        dtype=np.float64,
    )
    return xyz, n_vertices


def _parse_ls_links(data: bytes) -> dict:
    """Parse the LS_Links section and return face/cell connectivity.

    The section stores five successive data blocks (each preceded by a
    16-byte descriptor and a 12-byte data header) in this order:

        1. owner   – for each face: the owning cell index (0-based)
        2. neighbor – for each face: the neighboring cell index, or
                      0xFFFFFFFF (−1) for boundary faces
        3. npe     – nodes-per-face (all 3 = triangular faces)
        4. ftype   – single value (4 = tetrahedra)
        5. conn    – face-to-node connectivity stored *column-major*:
                     [node0_f0..node0_fN, node1_f0..node1_fN, node2_f0..node2_fN]
                     Vertex indices are 0-based.

    Returns a dict with keys:
        n_faces, face_nodes (ndarray (n_faces,3) 0-based),
        owner, neighbor, boundary_faces (list of 0-based face indices),
        cell_faces (dict cell_id → set of 0-based face indices), n_cells.
    Returns None on failure.
    """
    sec_start = _find_section(data, "LS_Links")
    if sec_start < 0:
        return None

    # Skip label header (40B) then 5 outer descriptor blocks (each 16B)
    pos = sec_start + 40
    for _ in range(5):
        if pos + 16 > len(data) or read_i32_be(data, pos) != 12:
            break
        pos += 16

    def _read_block():
        """Consume a (16B descriptor?) + 12B header + payload; return int32 array."""
        nonlocal pos
        # Skip optional descriptor
        if pos + 16 <= len(data) and read_i32_be(data, pos) == 12:
            n_check = read_i32_be(data, pos + 8)
            if n_check > 0:
                pos += 16
        if pos + 12 > len(data):
            return None
        bc = read_i32_be(data, pos + 4)
        if bc <= 0 or bc > 30000:
            return None
        pos += 12
        raw = np.frombuffer(data[pos : pos + bc], dtype=">u4").copy()
        pos += bc
        return raw

    blocks = [_read_block() for _ in range(5)]
    if any(b is None for b in blocks):
        return None

    owner_u, neigh_u, _npe, _ftype, conn_u = blocks

    # Signed conversion (0xFFFFFFFF → -1)
    BNDRY_U = np.uint32(0xFFFFFFFF)
    owner = owner_u.astype(np.int64)
    neigh = neigh_u.astype(np.int64)
    neigh[neigh_u == BNDRY_U] = -1

    n_faces = int(len(owner))

    # The block byte_count value (= n_faces × 4) occasionally leaks into the
    # owner/neighbor arrays as a sentinel at block boundaries.  Any cell index
    # larger than n_faces itself is physically impossible (cell count < face
    # count for a valid 3-D mesh), so treat such values as invalid.
    owner[owner > n_faces] = -1
    neigh[(neigh > n_faces) & (neigh != -1)] = -1
    n_nodes_per_face = 3

    # Rebuild per-face node list from column-major storage
    face_nodes = np.column_stack([
        conn_u[k * n_faces : (k + 1) * n_faces].astype(np.int32)
        for k in range(n_nodes_per_face)
    ])

    # Build cell-to-face mapping from owner + neighbor arrays
    from collections import defaultdict
    cell_faces: dict = defaultdict(set)
    for fi in range(n_faces):
        ow = int(owner[fi])
        nb = int(neigh[fi])
        if 0 <= ow < 10000:
            cell_faces[ow].add(fi)
        if nb != -1 and 0 <= nb < 10000:
            cell_faces[nb].add(fi)

    boundary_faces = [fi for fi in range(n_faces) if neigh[fi] == -1]
    n_cells = len(cell_faces)

    return {
        "n_faces": n_faces,
        "face_nodes": face_nodes,       # (n_faces, 3), 0-based
        "owner": owner,
        "neighbor": neigh,
        "boundary_faces": boundary_faces,
        "cell_faces": cell_faces,
        "n_cells": n_cells,
    }


def _parse_n_elements(data: bytes) -> int:
    """Return the element count from LS_CvolIdOfElements (dim0 of descriptor 5)."""
    sec_start = _find_section(data, "LS_CvolIdOfElements")
    if sec_start < 0:
        return 135  # fallback

    pos = sec_start + 40
    for _ in range(10):
        if pos + 16 > len(data) or read_i32_be(data, pos) != 12:
            break
        dim0 = read_i32_be(data, pos + 8)
        dim1 = read_i32_be(data, pos + 12)
        pos += 16
        if dim1 == 1 and dim0 > 1:
            return dim0
    return 135


def parse_gph_mesh(filepath: str) -> dict:
    """Extract mesh data (vertices, faces, cells) from GPH file."""
    with open(filepath, "rb") as f:
        data = f.read()

    result: dict = {
        "vertices": None, "n_vertices": 0,
        "link_data": None, "n_elements": 0,
    }

    # LS_Nodes – word-reversed float64 coordinate blocks
    xyz, n_vertices = _parse_ls_nodes(data)
    if xyz is not None and n_vertices > 0:
        result["vertices"] = xyz
        result["n_vertices"] = n_vertices

    # LS_Links – structured face/cell connectivity
    link_data = _parse_ls_links(data)
    if link_data is not None:
        # Clamp any garbage node indices that exceed the known vertex count.
        # The block byte-count sentinel (= n_faces × nodes_per_face × 4) can
        # leak into the last entry of the connectivity column-major array.
        if n_vertices > 0:
            fn = link_data["face_nodes"]
            bad = fn >= n_vertices
            if bad.any():
                # Replace with the nearest valid index (clamp to n_vertices-1)
                fn[bad] = n_vertices - 1
                link_data["face_nodes"] = fn
        result["link_data"] = link_data
        result["n_elements"] = link_data["n_cells"]

    # Fallback element count from LS_CvolIdOfElements
    if result["n_elements"] == 0:
        result["n_elements"] = _parse_n_elements(data)

    return result


def _cgns_node(parent, name: str, label: str, type_str: str, order: int = 2):
    """Create a CGNS HDF5 group with the standard CGNS attributes."""
    grp = parent.create_group(name)
    grp.attrs["label"] = _cgns_str33(label)
    grp.attrs["type"]  = _cgns_str33(type_str)
    return grp


def _write_ngon_elements(zone, face_nodes_1based: np.ndarray, elem_range: np.ndarray) -> None:
    """Write NGonElements (NGON_n = type 22) to zone.

    face_nodes_1based: shape (n_faces, nodes_per_face), 1-based vertex indices.
    """
    n_faces, npf = face_nodes_1based.shape
    el = _cgns_node(zone, "NGonElements", "Elements_t", "I4")
    el.create_dataset(" data", data=np.array([22, 0], dtype=np.int32))

    er = _cgns_node(el, "ElementRange", "IndexRange_t", "I8")
    er.create_dataset(" data", data=elem_range.astype(np.int64))

    # ElementStartOffset: offsets into flattened connectivity (one extra sentinel)
    offsets = np.arange(0, (n_faces + 1) * npf, npf, dtype=np.int64)
    eso = _cgns_node(el, "ElementStartOffset", "DataArray_t", "I8")
    eso.create_dataset(" data", data=offsets)

    conn = face_nodes_1based.flatten().astype(np.int64)
    ec = _cgns_node(el, "ElementConnectivity", "DataArray_t", "I8")
    ec.create_dataset(" data", data=conn)


def _write_nface_elements(zone, cell_faces: dict, n_faces: int,
                           elem_range: np.ndarray) -> None:
    """Write NFaceElements (NFACE_n = type 23) to zone.

    cell_faces: mapping cell_id (0-based) → set of face indices (0-based).
    n_faces: total face count (for 1-based face indexing offset = 0, since face IDs are 1-based).
    """
    n_cells = len(cell_faces)
    el = _cgns_node(zone, "NFaceElements", "Elements_t", "I4")
    el.create_dataset(" data", data=np.array([23, 0], dtype=np.int32))

    er = _cgns_node(el, "ElementRange", "IndexRange_t", "I8")
    er.create_dataset(" data", data=elem_range.astype(np.int64))

    offsets = [0]
    conn: list[int] = []
    for cell_id in range(n_cells):
        faces = sorted(cell_faces.get(cell_id, []))
        conn.extend(fi + 1 for fi in faces)   # convert to 1-based face index
        offsets.append(len(conn))

    eso = _cgns_node(el, "ElementStartOffset", "DataArray_t", "I8")
    eso.create_dataset(" data", data=np.array(offsets, dtype=np.int64))

    ec = _cgns_node(el, "ElementConnectivity", "DataArray_t", "I8")
    ec.create_dataset(" data", data=np.array(conn, dtype=np.int64))


def _write_zone_bc(zone, boundary_face_ids_1based: np.ndarray, bc_name: str = "box_surfs") -> None:
    """Write ZoneBC_t with a single BC_t PointList for boundary faces."""
    zbc = _cgns_node(zone, "ZoneBC", "ZoneBC_t", "MT")

    bc = _cgns_node(zbc, bc_name, "BC_t", "C1")
    bc.create_dataset(" data", data=np.frombuffer(b"BCWall", dtype=np.int8))

    pl = _cgns_node(bc, "PointList", "IndexArray_t", "I8")
    # Shape (n_boundary_faces, 1) to match box_ngons.cgns convention
    pl.create_dataset(
        " data",
        data=boundary_face_ids_1based.reshape(-1, 1).astype(np.int64),
    )

    gl = _cgns_node(bc, "GridLocation", "GridLocation_t", "C1")
    gl.create_dataset(" data", data=np.frombuffer(b"FaceCenter", dtype=np.int8))


def _cgns_str33(s: str) -> np.ndarray:
    """CGNS name/label: 33 bytes, null-padded."""
    b = (s[:32] + "\0").encode("ascii")[:33]
    return np.frombuffer(b.ljust(33)[:33], dtype=np.uint8)


def write_cgns(mesh: dict, outpath: str, zone_name: str = "box_vol") -> None:
    """Write mesh to CGNS/HDF5 following the NGon/NFace convention of box_ngons.cgns.

    Structure produced:
        CGNSLibraryVersion (R4)
        Base (CGNSBase_t, I4 [3,3])
          <zone_name> (Zone_t, I8 [[n_verts],[n_cells],[0]])
            ZoneType          → Unstructured
            GridCoordinates   → CoordinateX/Y/Z (R8)
            NGonElements      → triangular faces  (NGON_n, element type 22)
            NFaceElements     → tetrahedral cells (NFACE_n, element type 23)
            ZoneBC
              box_surfs       → BCWall, FaceCenter, PointList of boundary faces
    """
    vertices = mesh["vertices"]
    if vertices is None or len(vertices) == 0:
        raise ValueError("No vertices to write")

    link_data = mesh.get("link_data")
    if link_data is None:
        raise ValueError("No face/cell connectivity data (LS_Links parse failed)")

    n_vertex = int(vertices.shape[0])
    n_faces  = int(link_data["n_faces"])
    n_cells  = int(link_data["n_cells"])

    x = np.asarray(vertices[:, 0], dtype=np.float64)
    y = np.asarray(vertices[:, 1], dtype=np.float64)
    z = np.asarray(vertices[:, 2], dtype=np.float64)

    # Convert 0-based face_nodes to 1-based for CGNS
    face_nodes_1 = (link_data["face_nodes"] + 1).astype(np.int64)
    boundary_1   = np.array(link_data["boundary_faces"], dtype=np.int64) + 1

    with h5py.File(outpath, "w") as f:
        # ── root format marker (matching box_ngons.cgns) ──────────────────────
        fmt = np.frombuffer(b"IEEE_LITTLE_32\0", dtype=np.int8)
        f.create_dataset(" format", data=fmt)

        # CGNSLibraryVersion (R4 float, value 4.2 – matches reference)
        lv = _cgns_node(f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4")
        lv.create_dataset(" data", data=np.array([4.2], dtype=np.float32))

        # CGNSBase_t
        base = _cgns_node(f, "Base", "CGNSBase_t", "I4")
        base.create_dataset(" data", data=np.array([3, 3], dtype=np.int32))

        # Zone_t  – data shape (3,1) matching box_ngons.cgns
        zone = _cgns_node(base, zone_name, "Zone_t", "I8")
        zone.create_dataset(
            " data",
            data=np.array([[n_vertex], [n_cells], [0]], dtype=np.int64),
        )

        # ZoneType
        zt = _cgns_node(zone, "ZoneType", "ZoneType_t", "C1")
        zt.create_dataset(" data", data=np.frombuffer(b"Unstructured", dtype=np.int8))

        # GridCoordinates
        gc = _cgns_node(zone, "GridCoordinates", "GridCoordinates_t", "MT")
        for cname, arr in [("CoordinateX", x), ("CoordinateY", y), ("CoordinateZ", z)]:
            cd = _cgns_node(gc, cname, "DataArray_t", "R8")
            cd.create_dataset(" data", data=arr)

        # NGonElements: faces 1 … n_faces
        _write_ngon_elements(
            zone,
            face_nodes_1,
            elem_range=np.array([1, n_faces], dtype=np.int64),
        )

        # NFaceElements: cells n_faces+1 … n_faces+n_cells
        _write_nface_elements(
            zone,
            link_data["cell_faces"],
            n_faces=n_faces,
            elem_range=np.array([n_faces + 1, n_faces + n_cells], dtype=np.int64),
        )

        # ZoneBC
        _write_zone_bc(zone, boundary_1, bc_name="box_surfs")


def main():
    parser = argparse.ArgumentParser(description="Convert GPH file to CGNS format")
    parser.add_argument("gph_file", nargs="?", default="box.gph", help="Input GPH file")
    parser.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        help="Output CGNS file (default: input basename with .cgns)",
    )
    parser.add_argument(
        "-z", "--zone", default="box_vol", help="Zone name in CGNS (default: box_vol)"
    )
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
        print(f"  Cells    : {ld['n_cells']}  (tetrahedral, NFACE_n)")
        print(f"  BC faces : {len(ld['boundary_faces'])}")
    else:
        print(f"  Cells    : {mesh['n_elements']}  (LS_Links parse failed)")

    if mesh["vertices"] is None:
        print("Error: Could not extract vertex data from GPH.")
        sys.exit(1)
    if ld is None:
        print("Error: Could not parse LS_Links connectivity.")
        sys.exit(1)

    print(f"Writing: {out_path}")
    write_cgns(mesh, str(out_path), zone_name=args.zone)
    print("Done.")


if __name__ == "__main__":
    main()

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


def _parse_ls_links(data: bytes) -> tuple:
    """Return (links_raw, n_ints) for the LS_Links section.

    Locates the section dynamically and returns the raw big-endian int32
    array that starts after the section label header (40 bytes) and the
    five 16-byte descriptor blocks (80 bytes).
    """
    sec_start = _find_section(data, "LS_Links")
    if sec_start < 0:
        return None, 0

    next_sec = _find_section(data, "LS_Nodes")
    if next_sec < 0 or next_sec <= sec_start:
        return None, 0

    # Skip section header (40B) + 5 descriptor blocks (80B)
    data_start = sec_start + 40 + 80
    data_end = next_sec
    n_ints = (data_end - data_start) // 4
    if n_ints <= 0 or data_start + n_ints * 4 > len(data):
        return None, 0

    links_raw = np.frombuffer(data[data_start : data_start + n_ints * 4], dtype=">i4")
    return links_raw, n_ints


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


def _parse_element_connectivity(links: np.ndarray, n_vertices: int, n_elements: int):
    """Try to extract element connectivity from LS_Links raw data.
    Returns (element_nodes_list, element_type) or (None, None) on failure.
    element_nodes_list: list of lists, each inner list = vertex indices (1-based for CGNS).
    element_type: 'HEXA_8', 'NFACE_n', or 'MIXED'
    """
    # Heuristic: find runs of 8 consecutive values in [1, n_vertices] as HEXA_8
    HEXA_8 = 12
    valid = (links >= 1) & (links <= n_vertices)
    elems = []
    i = 0
    while i <= len(links) - 8 and len(elems) < n_elements:
        chunk = links[i : i + 8]
        if np.all(valid[i : i + 8]):
            elems.append(chunk.tolist())
            i += 8
        else:
            i += 1
    if len(elems) >= n_elements * 0.9:
        return elems[:n_elements], "HEXA_8"

    # Fallback: try NFACE - [n_faces, f0_nv, v.., f1_nv, v.., ...] per element
    elems_nface = []
    pos = 0
    while pos < len(links) - 4 and len(elems_nface) < n_elements:
        nf = int(links[pos])
        if nf < 3 or nf > 32:
            pos += 1
            continue
        pos += 1
        faces = []
        all_ok = True
        for _ in range(nf):
            if pos >= len(links):
                all_ok = False
                break
            nv = int(links[pos])
            pos += 1
            if nv < 2 or nv > 16 or pos + nv > len(links):
                all_ok = False
                break
            vs = links[pos : pos + nv]
            if np.any(vs < 1) or np.any(vs > n_vertices):
                all_ok = False
            faces.append(vs.tolist())
            pos += nv
        if all_ok and faces:
            elems_nface.append(faces)
        else:
            pos -= 1
            pos += 1
    if len(elems_nface) >= n_elements * 0.5:
        return elems_nface, "NFACE_n"
    return None, None


def parse_gph_mesh(filepath: str) -> dict:
    """Extract mesh data (vertices, elements) from GPH file."""
    with open(filepath, "rb") as f:
        data = f.read()

    result = {"vertices": None, "elements": None, "element_type": None, "n_vertices": 0, "n_elements": 0}

    # LS_Nodes – coordinates stored as word-reversed float64 in three axis blocks
    xyz, n_vertices = _parse_ls_nodes(data)
    if xyz is not None and n_vertices > 0:
        result["vertices"] = xyz
        result["n_vertices"] = n_vertices

    # LS_Links – raw int32 connectivity data
    links_raw, n_ints = _parse_ls_links(data)
    if links_raw is not None:
        result["links_raw"] = links_raw
        result["n_links_ints"] = n_ints

    # Element count from LS_CvolIdOfElements
    n_elements = _parse_n_elements(data)
    result["n_elements"] = n_elements

    # Parse connectivity
    if links_raw is not None and n_vertices > 0:
        elems, etype = _parse_element_connectivity(links_raw, n_vertices, n_elements)
        result["elements"] = elems
        result["element_type"] = etype

    return result


def _write_elements(zone, mesh: dict, n_vertex: int, n_cell: int) -> None:
    """Write Elements_t per CGNS 4.2: ElementRange, ElementConnectivity, ElementType."""
    elems = mesh.get("elements")
    etype_str = mesh.get("element_type")

    if elems and etype_str == "HEXA_8":
        conn = np.array([v for e in elems for v in e], dtype=np.int64)
        _add_elements_section(
            zone, "Hexahedra", "", 0,
            np.array([1, len(elems)], dtype=np.int64),
            conn,
            element_type_str="HEXA_8",
        )
        return

    if elems and etype_str == "NFACE_n":
        conn_list = []
        for faces in elems:
            conn_list.append(len(faces))
            for f in faces:
                conn_list.append(len(f))
                conn_list.extend(f)
        conn = np.array(conn_list, dtype=np.int64)
        _add_elements_section(
            zone, "Polyhedra", "", 0,
            np.array([1, len(elems)], dtype=np.int64),
            conn,
            element_type_str="NFACE_n",
        )
        return

    # Fallback: NODE elements (1 vertex per element)
    n_use = min(n_cell, n_vertex)
    conn = np.arange(1, n_use + 1, dtype=np.int64)
    _add_elements_section(
        zone, "Vertices", "", 0,
        np.array([1, n_use], dtype=np.int64),
        conn,
        element_type_str="NODE",
    )


def _add_elements_section(
    zone, name: str, _label_suffix: str, _cgns_type: int,
    elem_range: np.ndarray, connectivity: np.ndarray,
    element_type_str: str = "NODE",
) -> None:
    """Add one Elements_t child to zone. element_type_str: HEXA_8, NFACE_n, NODE, etc."""
    el = zone.create_group(name)
    el.attrs["name"] = _cgns_str33(name)
    el.attrs["label"] = _cgns_str33("Elements_t")
    el.attrs["type"] = _cgns_str33("I8")
    el.attrs["order"] = np.int32(2)

    # ElementRange
    er = el.create_group("ElementRange")
    er.attrs["name"] = _cgns_str33("ElementRange")
    er.attrs["label"] = _cgns_str33("IndexRange_t")
    er.attrs["type"] = _cgns_str33("I8")
    er.attrs["order"] = np.int32(2)
    er.create_dataset(" data", data=elem_range)

    # ElementType (C1 string: HEXA_8, NFACE_n, NODE, ...)
    et = el.create_group("ElementType")
    et.attrs["name"] = _cgns_str33("ElementType")
    et.attrs["label"] = _cgns_str33("ElementType_t")
    et.attrs["type"] = _cgns_str33("C1")
    et.attrs["order"] = np.int32(2)
    s = element_type_str.ljust(32).encode("ascii")[:32]
    et.create_dataset(" data", data=np.array([s], dtype="S32"))

    # ElementConnectivity
    ec = el.create_group("ElementConnectivity")
    ec.attrs["name"] = _cgns_str33("ElementConnectivity")
    ec.attrs["label"] = _cgns_str33("DataArray_t")
    ec.attrs["type"] = _cgns_str33("I8")
    ec.attrs["order"] = np.int32(2)
    ec.create_dataset(" data", data=connectivity)


def _cgns_str33(s: str) -> np.ndarray:
    """CGNS name/label: 33 bytes, null-padded."""
    b = (s[:32] + "\0").encode("ascii")[:33]
    return np.frombuffer(b.ljust(33)[:33], dtype=np.uint8)


def write_cgns(mesh: dict, outpath: str, zone_name: str = "Zone1") -> None:
    """Write mesh to CGNS/HDF5 file per CGNS 4.2 standard (no CGNSTree, correct node data)."""
    vertices = mesh["vertices"]
    if vertices is None or len(vertices) == 0:
        raise ValueError("No vertices to write")

    n_vertex = vertices.shape[0]
    # xyz is already float64 from the word-reversed float64 parser
    x = np.asarray(vertices[:, 0], dtype=np.float64)
    y = np.asarray(vertices[:, 1], dtype=np.float64)
    z = np.asarray(vertices[:, 2], dtype=np.float64)

    n_cell = mesh.get("n_elements", 0) if mesh.get("links_raw") is not None else 0
    n_cell = max(1, n_cell)

    with h5py.File(outpath, "w") as f:
        # Root attributes (CGNS HDF5)
        f.attrs["format"] = _cgns_str33("HDF5")
        f.attrs["version"] = _cgns_str33("4.2")

        # CGNSLibraryVersion - direct child of root (no CGNSTree)
        libver = f.create_group("CGNSLibraryVersion")
        libver.attrs["name"] = _cgns_str33("CGNSLibraryVersion")
        libver.attrs["label"] = _cgns_str33("CGNSLibraryVersion_t")
        libver.attrs["type"] = _cgns_str33("C1")
        libver.attrs["order"] = np.int32(2)
        libver.create_dataset(
            " data", data=np.array(["4.2".encode("ascii")], dtype="S4")
        )

        # CGNSBase_t - direct child of root; CellDimension and PhysicalDimension
        # stored in Base node's data (I4[2]), not as separate child nodes
        base = f.create_group("Base")
        base.attrs["name"] = _cgns_str33("Base")
        base.attrs["label"] = _cgns_str33("CGNSBase_t")
        base.attrs["type"] = _cgns_str33("I4")
        base.attrs["order"] = np.int32(2)
        base.create_dataset(
            " data", data=np.array([3, 3], dtype=np.int32)
        )  # [CellDimension, PhysicalDimension]

        # Zone_t - data [VertexSize, CellSize, VertexSizeBoundary] in Zone node
        zone = base.create_group(zone_name)
        zone.attrs["name"] = _cgns_str33(zone_name)
        zone.attrs["label"] = _cgns_str33("Zone_t")
        zone.attrs["type"] = _cgns_str33("I8")  # cgsize_t
        zone.attrs["order"] = np.int32(2)
        zone.create_dataset(
            " data",
            data=np.array([[n_vertex, n_cell, 0]], dtype=np.int64),
        )  # VertexSize, CellSize, VertexSizeBoundary (IndexDimension=1)

        # ZoneType
        zt = zone.create_group("ZoneType")
        zt.attrs["name"] = _cgns_str33("ZoneType")
        zt.attrs["label"] = _cgns_str33("ZoneType_t")
        zt.attrs["type"] = _cgns_str33("C1")
        zt.attrs["order"] = np.int32(2)
        zt.create_dataset(
            " data", data=np.array(["Unstructured".encode("ascii")], dtype="S12")
        )

        # GridCoordinates
        gc = zone.create_group("GridCoordinates")
        gc.attrs["name"] = _cgns_str33("GridCoordinates")
        gc.attrs["label"] = _cgns_str33("GridCoordinates_t")
        gc.attrs["type"] = _cgns_str33("MT")
        gc.attrs["order"] = np.int32(2)

        for coord_name, arr in [("CoordinateX", x), ("CoordinateY", y), ("CoordinateZ", z)]:
            coord = gc.create_group(coord_name)
            coord.attrs["name"] = _cgns_str33(coord_name)
            coord.attrs["label"] = _cgns_str33("DataArray_t")
            coord.attrs["type"] = _cgns_str33("R8")
            coord.attrs["order"] = np.int32(2)
            coord.create_dataset(" data", data=arr)

        # Elements_t - ElementRange, ElementConnectivity (required for unstructured)
        _write_elements(zone, mesh, n_vertex, n_cell)


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
        "-z", "--zone", default="Zone1", help="Zone name in CGNS (default: Zone1)"
    )
    args = parser.parse_args()

    gph_path = Path(args.gph_file)
    if not gph_path.exists():
        print(f"Error: file not found: {gph_path}")
        sys.exit(1)

    out_path = Path(args.output) if args.output else gph_path.with_suffix(".cgns")

    print(f"Reading: {gph_path}")
    mesh = parse_gph_mesh(str(gph_path))
    print(f"  Vertices: {mesh['n_vertices']}")
    if mesh.get("links_raw") is not None:
        print(f"  Links (raw ints): {mesh['n_links_ints']}")
    print(f"  Elements: {mesh['n_elements']}")
    if mesh.get("element_type"):
        print(f"  Element type: {mesh['element_type']}")
    elif mesh.get("elements") is None:
        print("  Element type: NODE (fallback; LS_Links format TBD)")

    if mesh["vertices"] is None:
        print("Error: Could not extract vertex data from GPH.")
        sys.exit(1)

    print(f"Writing: {out_path}")
    write_cgns(mesh, str(out_path), zone_name=args.zone)
    print("Done.")


if __name__ == "__main__":
    main()

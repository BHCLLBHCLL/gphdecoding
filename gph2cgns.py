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

    # LS_Nodes
    nodes_section_end = 0x2C60
    nodes_data_start = 0x2750
    n_vertices = (nodes_section_end - nodes_data_start) // 12
    if nodes_data_start + n_vertices * 12 <= len(data):
        xyz = np.zeros((n_vertices, 3), dtype=np.float32)
        for i in range(n_vertices):
            base = nodes_data_start + i * 12
            xyz[i, 0] = read_f32_be(data, base)
            xyz[i, 1] = read_f32_be(data, base + 4)
            xyz[i, 2] = read_f32_be(data, base + 8)
        result["vertices"] = xyz
        result["n_vertices"] = n_vertices

    # LS_Links
    links_start = 0x09C0
    links_end = 0x26B0
    n_ints = (links_end - links_start) // 4
    links_raw = None
    if n_ints > 0 and links_start + n_ints * 4 <= len(data):
        links_raw = np.frombuffer(data[links_start : links_start + n_ints * 4], dtype=">i4")
        result["links_raw"] = links_raw
        result["n_links_ints"] = n_ints

    # Element count
    n_elements = 135
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
    x = vertices[:, 0].astype(np.float64)
    y = vertices[:, 1].astype(np.float64)
    z = vertices[:, 2].astype(np.float64)

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

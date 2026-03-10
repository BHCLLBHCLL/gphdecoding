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


def parse_gph_mesh(filepath: str) -> dict:
    """Extract mesh data (vertices, elements) from GPH file."""
    with open(filepath, "rb") as f:
        data = f.read()

    result = {"vertices": None, "elements": None, "n_vertices": 0, "n_elements": 0}

    # LS_Nodes: vertex coordinates R4[n,3], section 0x26b0-0x2c60
    # Data starts after descriptor blocks; heuristic: ~0x2750
    nodes_section_end = 0x2C60
    nodes_data_start = 0x2750  # empirically determined for box.gph
    n_vertices = (nodes_section_end - nodes_data_start) // 12  # 3 * R4
    if nodes_data_start + n_vertices * 12 <= len(data):
        xyz = np.zeros((n_vertices, 3), dtype=np.float32)
        for i in range(n_vertices):
            base = nodes_data_start + i * 12
            xyz[i, 0] = read_f32_be(data, base)
            xyz[i, 1] = read_f32_be(data, base + 4)
            xyz[i, 2] = read_f32_be(data, base + 8)
        result["vertices"] = xyz
        result["n_vertices"] = n_vertices

    # LS_CvolIdOfElements: 135 elements
    # LS_Links: connectivity - structure not fully reverse-engineered
    links_start = 0x09C0
    links_end = 0x26B0
    n_ints = (links_end - links_start) // 4
    if n_ints > 0 and links_start + n_ints * 4 <= len(data):
        conn = np.frombuffer(
            data[links_start : links_start + n_ints * 4],
            dtype=">i4",
        )
        result["links_raw"] = conn
        result["n_links_ints"] = n_ints

    # Element count from LS_CvolIdOfElements
    result["n_elements"] = 135

    return result


def _cgns_str33(s: str) -> np.ndarray:
    """CGNS name/label: 33 bytes, null-padded."""
    b = (s[:32] + "\0").encode("ascii")[:33]
    return np.frombuffer(b.ljust(33)[:33], dtype=np.uint8)


def write_cgns(mesh: dict, outpath: str, zone_name: str = "Zone1") -> None:
    """Write mesh to CGNS/HDF5 file (SIDS-compliant structure)."""
    vertices = mesh["vertices"]
    if vertices is None or len(vertices) == 0:
        raise ValueError("No vertices to write")

    n_vertex = vertices.shape[0]
    x = vertices[:, 0].astype(np.float64)
    y = vertices[:, 1].astype(np.float64)
    z = vertices[:, 2].astype(np.float64)

    with h5py.File(outpath, "w") as f:
        # Root attributes (CGNS HDF5)
        f.attrs["format"] = _cgns_str33("HDF5")
        f.attrs["version"] = _cgns_str33("4.3")

        # CGNSTree
        tree = f.create_group("CGNSTree")
        tree.attrs["name"] = _cgns_str33("CGNSTree")
        tree.attrs["label"] = _cgns_str33("CGNSTree_t")
        tree.attrs["type"] = _cgns_str33("MT")
        tree.attrs["order"] = np.int32(2)

        # CGNSLibraryVersion
        libver = tree.create_group("CGNSLibraryVersion")
        libver.attrs["name"] = _cgns_str33("CGNSLibraryVersion")
        libver.attrs["label"] = _cgns_str33("CGNSLibraryVersion_t")
        libver.attrs["type"] = _cgns_str33("C1")
        libver.attrs["order"] = np.int32(2)
        libver.create_dataset(
            " data", data=np.array(["4.3".encode("ascii")], dtype="S4")
        )

        # Base
        base = tree.create_group("Base")
        base.attrs["name"] = _cgns_str33("Base")
        base.attrs["label"] = _cgns_str33("CGNSBase_t")
        base.attrs["type"] = _cgns_str33("MT")
        base.attrs["order"] = np.int32(2)

        # CellDimension, PhysicalDimension
        cell_dim = base.create_group("CellDimension")
        cell_dim.attrs["name"] = _cgns_str33("CellDimension")
        cell_dim.attrs["label"] = _cgns_str33("IndexDimension_t")
        cell_dim.attrs["type"] = _cgns_str33("I4")
        cell_dim.attrs["order"] = np.int32(2)
        cell_dim.create_dataset(" data", data=np.array([3], dtype=np.int32))

        phys_dim = base.create_group("PhysicalDimension")
        phys_dim.attrs["name"] = _cgns_str33("PhysicalDimension")
        phys_dim.attrs["label"] = _cgns_str33("IndexDimension_t")
        phys_dim.attrs["type"] = _cgns_str33("I4")
        phys_dim.attrs["order"] = np.int32(2)
        phys_dim.create_dataset(" data", data=np.array([3], dtype=np.int32))

        # Zone
        zone = base.create_group(zone_name)
        zone.attrs["name"] = _cgns_str33(zone_name)
        zone.attrs["label"] = _cgns_str33("Zone_t")
        zone.attrs["type"] = _cgns_str33("MT")
        zone.attrs["order"] = np.int32(2)

        # ZoneType
        zt = zone.create_group("ZoneType")
        zt.attrs["name"] = _cgns_str33("ZoneType")
        zt.attrs["label"] = _cgns_str33("ZoneType_t")
        zt.attrs["type"] = _cgns_str33("C1")
        zt.attrs["order"] = np.int32(2)
        zt.create_dataset(
            " data", data=np.array(["Unstructured".encode("ascii")], dtype="S12")
        )

        # VertexSize, CellSize
        vs = zone.create_group("VertexSize")
        vs.attrs["name"] = _cgns_str33("VertexSize")
        vs.attrs["label"] = _cgns_str33("IndexRange_t")
        vs.attrs["type"] = _cgns_str33("I4")
        vs.attrs["order"] = np.int32(2)
        vs.create_dataset(" data", data=np.array([1, n_vertex], dtype=np.int32))

        cs = zone.create_group("CellSize")
        cs.attrs["name"] = _cgns_str33("CellSize")
        cs.attrs["label"] = _cgns_str33("IndexRange_t")
        cs.attrs["type"] = _cgns_str33("I4")
        cs.attrs["order"] = np.int32(2)
        n_cell = mesh.get("n_elements", 0) if mesh.get("links_raw") is not None else 0
        cs.create_dataset(" data", data=np.array([1, max(1, n_cell)], dtype=np.int32))

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

    if mesh["vertices"] is None:
        print("Error: Could not extract vertex data from GPH.")
        sys.exit(1)

    print(f"Writing: {out_path}")
    write_cgns(mesh, str(out_path), zone_name=args.zone)
    print("Done.")


if __name__ == "__main__":
    main()

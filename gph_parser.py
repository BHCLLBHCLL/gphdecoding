#!/usr/bin/env python3
"""
GPH (Geometry/Polyhedron) Binary Format Parser and Reverse-Engineered Format Describer.

Parses box.gph and outputs a structured format description.
Byte order: Big-Endian for all multi-byte integers and metadata.
"""

import struct
import sys
from pathlib import Path


def read_i32_be(data: bytes, pos: int) -> int:
    """Read 4-byte big-endian integer."""
    return int.from_bytes(data[pos : pos + 4], "big")


def read_f32_be(data: bytes, pos: int) -> float:
    """Read 4-byte big-endian float."""
    return struct.unpack(">f", data[pos : pos + 4])[0]


def parse_gph(filepath: str) -> dict:
    """Parse GPH file and return structured description."""
    with open(filepath, "rb") as f:
        data = f.read()

    result = {
        "file_size": len(data),
        "byte_order": "big-endian",
        "sections": [],
        "header": {},
        "data_arrays": {},
    }

    pos = 0

    # --- File header ---
    rec0_len = read_i32_be(data, 0)
    rec0_data = data[4 : 4 + rec0_len].decode("ascii", errors="replace")
    result["header"]["format_id"] = rec0_data  # "CRDL-FLD"
    pos = 4 + rec0_len

    # --- Initial dimensions (often 4,4,4) ---
    v1 = read_i32_be(data, pos)
    v2 = read_i32_be(data, pos + 4)
    v3 = read_i32_be(data, pos + 8)
    result["header"]["dims"] = [v1, v2, v3]
    pos += 12

    # --- Named fields: pattern [0x20][32-char name][0x20][payload] ---
    FIELD_NAMES = [
        "FileRevision",
        "Application",
        "ApplicationVersion",
        "ReleaseDate",
        "GridType",
        "Dimension",
        "Bias",
        "Date",
        "Comments",
        "Cycle",
        "Unused",
        "Encoding",
        "UTF-8",
        "HeaderDataEnd",
        "OverlapStart_0",
        "LS_CvolIdOfElements",
        "LS_Links",
        "LS_Nodes",
        "LS_SurfaceRegions",
        "Element_InformationFlag",
        "OverlapEnd",
    ]

    # Locate fields by scanning for 32-byte blocks with known names
    def find_field(name: str):
        name_padded = name.ljust(32).encode("ascii")
        idx = data.find(name_padded)
        return idx - 4 if idx >= 4 and read_i32_be(data, idx - 4) == 32 else None

    # Parse scalar/metadata fields
    scalar_values = {}
    for name in [
        "FileRevision",
        "Application",
        "ApplicationVersion",
        "ReleaseDate",
        "GridType",
        "Dimension",
        "Bias",
        "Date",
        "Cycle",
        "Encoding",
    ]:
        name_pos = find_field(name)
        if name_pos and name_pos > 0:
            # Value typically 36 bytes after name start (skip 0x20 + name + 0x20 + 12-byte descriptor)
            val_pos = name_pos + 4 + 32 + 4 + 12 + 8  # heuristic
            if val_pos + 4 <= len(data):
                scalar_values[name] = read_i32_be(data, val_pos)

    # Extract string values by heuristics
    str_fields = {}
    app_pos = find_field("Application")
    if app_pos:
        # "SCTpre" appears at 0xb8 in earlier dump
        s = data[0xB4 : 0xB4 + 8].decode("ascii", errors="replace").strip()
        if s:
            str_fields["Application"] = s

    result["header"]["scalar_metadata"] = scalar_values
    result["header"]["string_fields"] = str_fields

    # --- Section layout (from reverse engineering) ---
    sections = [
        (0x0000, 0x001c, "file_header", "CRDL-FLD identifier + dims"),
        (0x001c, 0x0078, "FileRevision", "I4 scalar, e.g. 2025"),
        (0x0078, 0x00d8, "Application", "C1[8] e.g. SCTpre"),
        (0x00d8, 0x0134, "ApplicationVersion", "I4"),
        (0x0134, 0x0190, "ReleaseDate", "string/date"),
        (0x0190, 0x01ec, "GridType", "string"),
        (0x01ec, 0x0248, "Dimension", "I4"),
        (0x0248, 0x02a4, "Bias", "I4"),
        (0x02a4, 0x0300, "Date", "string"),
        (0x0300, 0x03a8, "Comments", "string"),
        (0x03a8, 0x04e0, "Cycle", "I4 + Unit"),
        (0x04e0, 0x0560, "Unused", "reserved"),
        (0x0560, 0x05d8, "Encoding", "e.g. UTF-8"),
        (0x05d8, 0x0600, "HeaderDataEnd", "marker"),
        (0x0600, 0x0628, "OverlapStart_0", "overlap marker"),
        (0x0628, 0x08dc, "LS_CvolIdOfElements", "I4[135] control volume IDs"),
        (0x08dc, 0x26b0, "LS_Links", "I4[] link/connectivity pairs"),
        (0x26b0, 0x2c60, "LS_Nodes", "R8[n,3] vertex coords (word-reversed float64, 3 axis blocks)"),
        (0x2c60, 0x42b0, "LS_SurfaceRegions", "surface region data"),
        (0x42b0, 0x45a4, "Element_InformationFlag", "element flags"),
        (0x45a4, len(data), "OverlapEnd", "trailer"),
    ]

    for start, end, name, desc in sections:
        result["sections"].append(
            {
                "offset_hex": f"0x{start:04X}",
                "end_hex": f"0x{end:04X}",
                "size": end - start,
                "name": name,
                "description": desc,
            }
        )

    # --- Data array extraction for validation ---
    # LS_CvolIdOfElements: n x I4, located dynamically
    # Structure: 40B label-header + 5×16B descriptors + 12B data-header + n×4B data
    cv_sec = b"LS_CvolIdOfElements" + b" " * 13  # 32 chars total
    cv_idx = data.find(cv_sec)
    if cv_idx >= 4 and read_i32_be(data, cv_idx - 4) == 32:
        cv_pos = cv_idx - 4 + 40  # skip section header
        # Skip 5 descriptors (each 16B) to find the one with dim0 = n_elements
        n_cv = 0
        for _ in range(10):
            if cv_pos + 16 > len(data) or read_i32_be(data, cv_pos) != 12:
                break
            d0 = read_i32_be(data, cv_pos + 8)
            d1 = read_i32_be(data, cv_pos + 12)
            cv_pos += 16
            if d1 == 1 and d0 > 1:
                n_cv = d0
                break
        if n_cv > 0:
            cv_pos += 12  # skip 12B data header
            arr = [read_i32_be(data, cv_pos + i * 4) for i in range(min(n_cv, 10))]
            result["data_arrays"]["LS_CvolIdOfElements_sample"] = arr
            result["data_arrays"]["LS_CvolIdOfElements_count"] = n_cv

    # LS_Nodes: vertex coordinates stored as three separate float64 axis blocks
    # Each block: [16B descriptor: 12/8/n/1] + [12B header: 12/bytes/upper_max]
    #             + [n × 8B word-reversed float64]
    # Word-reversed means each float64 is stored as [lower_32bit_BE][upper_32bit_BE].
    # e.g. 0.01 = 0x3F847AE147AE147B is stored as bytes 47 AE 14 7B 3F 84 7A E1.
    nodes_label = b"LS_Nodes" + b" " * 24
    nodes_idx = data.find(nodes_label)
    if nodes_idx >= 4:
        pos = nodes_idx + 32 + 4 + 4  # skip name(32) + second marker(4) + first desc(4 already counted)
        # Actually skip past label header: 4(marker)+32(name)+4(marker) = 40B from section start
        pos = nodes_idx - 4 + 40  # section_start + 40
        # Scan to find first [12, 8, n, 1] descriptor
        n_nodes = 0
        while pos + 16 <= len(data):
            if int.from_bytes(data[pos:pos+4], "big") != 12:
                break
            tc = int.from_bytes(data[pos+4:pos+8], "big")
            d0 = int.from_bytes(data[pos+8:pos+12], "big")
            d1 = int.from_bytes(data[pos+12:pos+16], "big")
            pos += 16
            if tc == 8 and d1 == 1 and d0 > 0:
                n_nodes = d0
                break
        if n_nodes > 0:
            sample = []
            axis_vals = []
            for axis in range(min(3, 3)):
                if pos + 12 > len(data):
                    break
                pos += 12  # skip 12B header
                vals = []
                for i in range(n_nodes):
                    if pos + 8 > len(data):
                        break
                    lower = int.from_bytes(data[pos:pos+4], "big")
                    upper = int.from_bytes(data[pos+4:pos+8], "big")
                    vals.append(struct.unpack(">d", ((upper << 32) | lower).to_bytes(8, "big"))[0])
                    pos += 8
                axis_vals.append(vals)
                if pos + 16 <= len(data) and int.from_bytes(data[pos:pos+4], "big") == 12:
                    pos += 16
            if len(axis_vals) == 3:
                for i in range(min(3, n_nodes)):
                    sample.append((axis_vals[0][i], axis_vals[2][i], axis_vals[1][i]))
                result["data_arrays"]["LS_Nodes_sample"] = sample
                result["data_arrays"]["LS_Nodes_count"] = n_nodes

    return result


def format_description() -> str:
    """Return the GPH format specification as text."""
    return """
================================================================================
                    GPH Binary Format Description (Reverse-Engineered)
================================================================================

1. OVERVIEW
-----------
  GPH appears to be a geometry/polyhedron mesh file format, likely from SCTpre
  or a CGNS-related CFD tool. The magic identifier is "CRDL-FLD" (8 bytes).
  Byte order: BIG-ENDIAN for all integers and floats.

2. FILE STRUCTURE
-----------------

  +------------------+----------------------------------------------------------+
  | Offset           | Content                                                  |
  +------------------+----------------------------------------------------------+
  | 0x0000 - 0x001C  | File header                                              |
  | 0x001C - 0x05D8  | Metadata section (named key-value fields)                 |
  | 0x05D8 - 0x0600  | HeaderDataEnd marker                                     |
  | 0x0600 - 0x0628  | OverlapStart_0                                           |
  | 0x0628 - 0x08DC  | LS_CvolIdOfElements (control volume IDs, I4 array)        |
  | 0x08DC - 0x26B0  | LS_Links (element connectivity, I4 array)                 |
  | 0x26B0 - 0x2C60  | LS_Nodes (vertex coordinates, R4[ n*3 ])                 |
  | 0x2C60 - 0x42B0  | LS_SurfaceRegions                                        |
  | 0x42B0 - 0x45A4  | Element_InformationFlag                                  |
  | 0x45A4 - EOF     | OverlapEnd trailer                                       |
  +------------------+----------------------------------------------------------+

3. RECORD FORMAT
----------------

  Each named field follows:

    [length: I4] = 0x20 (32)
    [name: C1[32]] = 32-byte ASCII name, space-padded
    [length: I4] = 0x20 (32)  -- optional, section marker
    [descriptor: variable]
      - 0x0C 0x04 [dim0] [dim1]  : dimensions, data type 04 = I4
      - 0x0C 0x08 [dim0] [dim1]  : data type 08 = R4 or I8
    [value / array data]

  Data type codes (inferred):
    - 04 : I4 (32-bit signed integer)
    - 08 : R4 (32-bit float) or I8 (64-bit int) depending on context

4. METADATA FIELDS (Header)
---------------------------

  | Field               | Type      | Example / Notes                    |
  |---------------------|-----------|-------------------------------------|
  | FileRevision        | I4        | 2025                                |
  | Application         | C1[8]     | "SCTpre"                            |
  | ApplicationVersion  | I4        | 1                                   |
  | ReleaseDate         | string    | date                                |
  | GridType            | string    | grid type                           |
  | Dimension           | I4        | spatial dimension (e.g. 3)          |
  | Bias                | I4        |                                    |
  | Date                | string    |                                    |
  | Comments            | string    |                                    |
  | Cycle               | I4        | simulation cycle                    |
  | Encoding            | C1        | "UTF-8"                             |

5. DATA ARRAYS
--------------

  LS_CvolIdOfElements  : I4[n]    - control volume ID per element (n from descriptor)
  LS_Links             : I4[]     - element connectivity (raw big-endian int32)
  LS_Nodes             : R8[n,3]  - vertex XYZ as three word-reversed float64 blocks
  LS_SurfaceRegions    : variable - surface boundary region definitions
  Element_InformationFlag : flags per element

  LS_Nodes detail:
    Each spatial axis (X, Z, Y in file order) is stored as one block:
      [16B descriptor: I4=12 / I4=8 / I4=n_verts / I4=1]
      [12B header:     I4=12 / I4=byte_count / I4=upper_half_of_max_coord]
      [n_verts × 8B]  one word-reversed float64 per vertex

    Word-reversed float64: each 8-byte double is stored as
      [lower_32bit_word_BE][upper_32bit_word_BE]
    Example: 0.01 (= 0x3F847AE147AE147B) is stored as bytes
      47 AE 14 7B  3F 84 7A E1
    When mis-read as two float32 values this yields (89128.96, 1.035),
    which is the root cause of the "89000+" coordinate bug.

6. ALIGNMENT & PADDING
----------------------
  - Records appear 4-byte aligned.
  - Names are fixed 32 bytes.
  - Block markers use 0x20 (32) as length prefix.

7. REFERENCES
-------------
  - Similar to ADF (Advanced Data Format) used in CGNS.
  - 32-char labels match ADF node label convention.
  - CRDL-FLD may mean "Card/Record Field" or vendor-specific.

================================================================================
"""


def main():
    args = sys.argv[1:]
    path = Path(args[0]) if args and not args[0].startswith("-") else Path(__file__).parent / "box.gph"

    if not path.exists():
        print(f"Error: file not found: {path}")
        sys.exit(1)

    print("=== GPH Format Parser ===\n")
    parsed = parse_gph(str(path))
    print(f"File: {path}")
    print(f"Size: {parsed['file_size']} bytes")
    print(f"Format ID: {parsed['header']['format_id']}")
    print(f"Header dims: {parsed['header']['dims']}")
    print()
    print("Sections:")
    for s in parsed["sections"]:
        print(f"  {s['offset_hex']}-{s['end_hex']} ({s['size']:5} B) {s['name']}: {s['description']}")
    print()
    if parsed["data_arrays"]:
        print("Data samples:")
        for k, v in parsed["data_arrays"].items():
            print(f"  {k}: {v}")
    print()
    print(format_description())


if __name__ == "__main__":
    main()

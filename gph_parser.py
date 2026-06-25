#!/usr/bin/env python3
"""
GPH (Geometry/Polyhedron) Binary Format Parser and Reverse-Engineered Format Describer.

Parses any .gph file and outputs a structured format description.
Byte order: Big-Endian for all multi-byte integers and metadata.
"""

import sys
from pathlib import Path

import numpy as np

from gph_model import (
    find_section,
    open_gph_buffer,
    parse_ls_assemblies_summary,
    parse_ls_cvol_ids,
    parse_ls_links_summary,
    parse_ls_nodes_vertices,
    parse_ls_parts,
    parse_ls_string_list,
    parse_ls_surface_regions_summary,
    format_part_cvol_spec,
    read_i32_be,
)


def parse_gph(filepath: str) -> dict:
    """Parse GPH file and return structured description."""
    with open_gph_buffer(filepath) as data:
        return _parse_gph_buffer(data, filepath)


def _parse_gph_buffer(data, filepath: str) -> dict:
    result = {
        "file_size": len(data),
        "byte_order": "big-endian",
        "sections": [],
        "header": {},
        "data_arrays": {},
        "mesh_summary": {},
        "partition": {},
    }

    # --- File header ---
    rec0_len = read_i32_be(data, 0)
    rec0_data = data[4 : 4 + rec0_len].decode("ascii", errors="replace")
    result["header"]["format_id"] = rec0_data
    v1 = read_i32_be(data, rec0_len + 4)
    v2 = read_i32_be(data, rec0_len + 8)
    v3 = read_i32_be(data, rec0_len + 12)
    result["header"]["dims"] = [v1, v2, v3]

    # --- Dynamic section layout ---
    section_names = [
        "FileRevision", "Application", "ApplicationVersion", "ReleaseDate",
        "GridType", "Dimension", "Bias", "Date", "Comments", "Cycle",
        "Unused", "Encoding", "HeaderDataEnd", "OverlapStart_0",
        "LS_CvolIdOfElements", "LS_Links", "LS_Nodes", "LS_SurfaceRegions",
        "LS_SolverUnusedRegions", "LS_VolumeRegions", "LS_Parts",
        "LS_Assemblies", "Element_InformationFlag", "OverlapEnd",
    ]
    found = []
    for name in section_names:
        off = find_section(data, name)
        if off >= 0:
            found.append((off, name))
    found.sort(key=lambda x: x[0])

    first_off = found[0][0] if found else len(data)
    layout = [(0, first_off, "file_header", "CRDL-FLD identifier + dims")]
    for i, (off, name) in enumerate(found):
        end = found[i + 1][0] if i + 1 < len(found) else len(data)
        layout.append((off, end, name, ""))

    for start, end, name, desc in layout:
        result["sections"].append({
            "offset_hex": f"0x{start:04X}",
            "end_hex": f"0x{end:04X}",
            "size": end - start,
            "name": name,
            "description": desc or _section_blurb(name),
        })

    # --- Mesh / partition samples ---
    sample, dialect, n_verts = parse_ls_nodes_vertices(data)
    if n_verts:
        result["data_arrays"]["LS_Nodes_count"] = n_verts
        result["data_arrays"]["LS_Nodes_dialect"] = dialect
        result["data_arrays"]["LS_Nodes_sample"] = sample

    links = parse_ls_links_summary(data)
    if links:
        result["mesh_summary"] = links

    # Parse cvol_id first — LS_Parts mapping uses its unique set as the authority.
    cvol = parse_ls_cvol_ids(data)

    parts = parse_ls_parts(data, cvol_id=cvol)
    if parts:
        result["partition"]["LS_Parts"] = [
            f"{name} (cvol={format_part_cvol_spec(cv)})" for name, cv in parts
        ]

    regions = parse_ls_string_list(data, "LS_VolumeRegions")
    if regions:
        result["partition"]["LS_VolumeRegions"] = regions

    surf = parse_ls_surface_regions_summary(data)
    if surf:
        result["partition"]["LS_SurfaceRegions"] = surf

    asm = parse_ls_assemblies_summary(data)
    if asm.get("part_paths") or asm.get("root_empty_prefix"):
        result["partition"]["LS_Assemblies"] = asm

    if cvol is not None:
        n_cv = len(cvol)
        sample = cvol[:10].tolist()
        result["data_arrays"]["LS_CvolIdOfElements_count"] = n_cv
        result["data_arrays"]["LS_CvolIdOfElements_sample"] = sample
        unique = sorted({int(x) for x in np.unique(cvol[:min(n_cv, 1_000_000)])})
        result["data_arrays"]["LS_CvolIdOfElements_unique"] = unique[:20]

    return result


def _section_blurb(name: str) -> str:
    blurbs = {
        "LS_CvolIdOfElements": "I4[n_cells] per-cell cvol_id (opaque Part id, not list index)",
        "LS_Links": "face topology: owner, neighbor, npe, conn (CSR; may split >1 GiB, multi-chunk)",
        "LS_Nodes": "R4/R8[n,3] vertex coords (three axis blocks; float32 or float64 BE)",
        "LS_SurfaceRegions": "named BC regions -> global face id lists",
        "LS_VolumeRegions": "volume region names (-> CGNS zones)",
        "LS_Parts": "part names + cvol spec (single id or membership list; see format_description §10)",
        "LS_Assemblies": "XML assembly tree for zone naming",
        "LS_SolverUnusedRegions": "solver-internal region names",
        "Element_InformationFlag": "per-element flags",
    }
    return blurbs.get(name, "")


def format_description() -> str:
    """Return the GPH format specification as text."""
    return """
================================================================================
                    GPH Binary Format Description (Reverse-Engineered)
================================================================================

1. OVERVIEW
-----------
  GPH is a geometry/polyhedron mesh file from Software Cradle scFLOW / SCTpre
  (and ANSA export pipelines).  Magic identifier: "CRDL-FLD" (8 bytes).
  Byte order: BIG-ENDIAN for all integers and floats.

2. FILE STRUCTURE (logical order)
---------------------------------
  Sections are located dynamically by scanning for 32-byte ASCII labels
  preceded by [I4=32].  Typical order:

    file_header, metadata fields, HeaderDataEnd, OverlapStart_0,
    LS_CvolIdOfElements, LS_Links, LS_Nodes,
    LS_SurfaceRegions, LS_SolverUnusedRegions, LS_VolumeRegions,
    LS_Parts, LS_Assemblies, Element_InformationFlag, OverlapEnd

  Offsets differ per file (box.gph, box_ansa.gph, tr03.gph, laptop_*.gph).

3. GENERIC DATA-BLOCK FORMAT
----------------------------
  Inside each named section, payloads use:

    [I4=12][I4=byte_count][payload][I4=byte_count]   (8-byte header + trailer)

  Interleaved 16-byte descriptors: [12, type_code, dim0, dim1]
    type 4 = I4, type 8 = float64 (R8) in coordinate blocks.

4. LS_Nodes - vertex coordinates
--------------------------------
  Three equal-sized axis blocks (X, Y, Z file order).  Type from 16-byte
  descriptors: type 4 = float32 (R4), type 8 = float64 (R8).  FPH meshes
  (e.g. tests/tr03_9.fph) use float32; ANSA exports use float64.

  gph_model.parse_ls_nodes_xyz() auto-selects among:
    - big-endian float32 (descriptor type 4)
    - standard big-endian float64
    - word-reversed float64 (legacy; disk order X,Z,Y -> permuted to X,Y,Z)

  Scoring uses _COORD_MIN_ABSMAX = 1e-4 m so float32 payload misread as
  float64 (~1e-13 denormals) loses to the correct decode.  Vertex count
  comes from descriptors (dim0 where dim0 > 1), not payload byte_count // 8.

5. LS_Links - face / cell topology
----------------------------------
  At least three equal I4 blocks of length n_faces:
    owner, neighbor (0xFFFFFFFF = boundary), npe (nodes per face).
  Optional 4-byte face_type block (legacy SCTpre).
  conn block: sum(npe) vertex indices, 0-based, CSR layout:
    face i nodes = conn[face_offsets[i]:face_offsets[i+1]].

  Pure-triangle meshes: all npe=3.  Mixed/polyhedral (e.g. tr03.gph):
  npe varies (3..11+).  Voxel meshes (laptop_simplified_voxel_v4.gph) use
  npe in {4,5,6,7} with ~10^7 faces and ~10^7 cells.  Denser meshes
  (laptop_simplified_denser_v2_gph.gph, ~5.9 GiB) reach ~10^8 faces.
  laptop_simplified_voxel_v6.gph (~4.9 GiB) has ~10^8 faces and ~5×10^7 cells.

  Large conn arrays (>~1 GiB): when sum(npe)*4 exceeds a single payload
  limit (~1073741824 bytes), scFLOW splits conn into multiple segments:

    [12, bc1][conn part 1][bc1]              <- standard block (bc1 often 1 GiB)
    [I4=1073741824][conn part 2 raw...]      <- bare byte_count + payload (no [I4=12])
    ... further 1 GiB bare chunks as needed ...
    [I4=1073741824][I4=bcN][conn final...]   <- last chunk (3+ segments): repeat
                                               1 GiB marker + actual payload bytes
                                               + data; misreading bcN as a vertex
                                               index causes face warping at the
                                               2nd GiB boundary (see DEV_SUMMARY
                                               §11.13).  Two-segment files use
                                               single-header [I4=bcN][payload].

  Conn block selection: if no block matches sum(npe) exactly, pick the
  largest non-triple I4 block with byte_count >= 12 (not 3*n_faces*4, which
  fails when the first conn segment is capped at 1 GiB on polyhedral meshes).

  gph2cgns / gph_model concatenate all segments before CSR indexing.
  parse_ls_links_summary reports conn_got, conn_chunks and conn_complete.

  Examples:
    laptop_simplified_voxel_v4.gph  - 2 conn chunks (~1.44 GiB total)
    laptop_simplified_denser_v2_gph.gph - 3 conn chunks (~2.05 GiB total)
    laptop_simplified_voxel_v6.gph - 2 conn chunks (~1.91 GiB total)

  gph2cgns imports _read_conn_continuations from gph_model (returns
  got, pos, n_continuations — callers must unpack all three).

  Files >512 MiB: gph2cgns, gph_parser and gphviewer memory-map the file
  instead of loading it entirely into RAM.

9. LS_CvolIdOfElements - per-cell partition id
----------------------------------------------
  I4[n_cells]: one cvol_id per volume cell — the opaque Part label stored in
  LS_Parts (NOT the 1-based index in the Part list).
  Example: laptop_simplified_voxel_less.gph uses cvol_ids {1, 9, 11}.
  Multi-region laptop models may use dozens of distinct cvol_ids (geometry
  sub-blocks); a single Part (e.g. air_domain) may own many of them via an
  explicit membership list (see §10).

10. LS_Parts - part definitions
-------------------------------
  Repeated records: 255-byte ASCII name block + post-name byte region.

  **Simple Part** (one cvol_id): post-name chain is typically [1, cvol_id]
  — leading 1 is a marker; the trailing value is the opaque scFLOW id.
  Examples: outlet11 -> [1,2]; rotation1 -> [1,7].

  **Composite Part** (many cvol_ids): background fluid parts such as
  air_domain in laptop_simplified_more_regions.gph use:
    [12,4,N,4]   N = length of the following list (NOT a cvol_id)
    I4[N]        explicit list of cvol_ids belonging to this Part
  Zone cells = all entries in LS_CvolIdOfElements whose value is in that set.

  Resolution (gph_model.parse_ls_parts -> PartCvolSpec = int | frozenset):
  1. Parse LS_CvolIdOfElements first; unique values form authoritative set S.
  2. For each Part, if [12,4,N,4] + I4[N] membership is present, return
     frozenset(members).
  3. Else pick the last [12,4,X,4] chain value in S (simple Part rule).
  4. gph2cgns / gphviewer use part_cvol_cell_mask() for zone cell masks.

  Helpers: format_part_cvol_spec(), part_cvol_cell_mask().
  Test: python tests/test_volume_zone_cells.py -v tests/laptop_simplified_more_regions.gph

11. LS_VolumeRegions / LS_Assemblies / LS_SurfaceRegions
--------------------------------------------------------
  LS_VolumeRegions: ASCII names -> one CGNS Zone_t each (FluidRegion, etc.).
  LS_Assemblies: XML tree; drives FPHPARTS.* / dotted zone names.
  LS_SurfaceRegions: triplets (name, face_ids I4[], weights I4[]) -> ZoneBC
  families (one BC_t per region; empty PointList when zone has no faces).

12. CGNS export (gph2cgns.py)
----------------------------
  Writes FLDUTIL-compatible CGNS/HDF5: NGON_n faces, NFACE_n cells (signed
  face ids), multi-zone layout from partition metadata, per-zone BC families.
  HDF5 superblock v0 (libver earliest) for ANSA compatibility - see
  GPH_FORMAT_SPEC.md section 9 and DEV_SUMMARY.md.

================================================================================
"""


def _default_sample_gph() -> Path:
    """Prefer tests/box_ansa.gph, then legacy repo-root names."""
    root = Path(__file__).resolve().parent
    for name in ("tests/box_ansa.gph", "box_ansa.gph", "box.gph"):
        p = root / name
        if p.is_file():
            return p
    return root / "tests" / "box_ansa.gph"


def main():
    args = sys.argv[1:]
    path = Path(args[0]) if args and not args[0].startswith("-") else _default_sample_gph()

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
        desc = f": {s['description']}" if s["description"] else ""
        print(f"  {s['offset_hex']}-{s['end_hex']} ({s['size']:5} B) {s['name']}{desc}")
    print()
    if parsed.get("mesh_summary"):
        print("Mesh topology:")
        for k, v in parsed["mesh_summary"].items():
            print(f"  {k}: {v}")
        print()
    if parsed.get("partition"):
        print("Partition metadata:")
        for k, v in parsed["partition"].items():
            print(f"  {k}: {v}")
        print()
    if parsed["data_arrays"]:
        print("Data samples:")
        for k, v in parsed["data_arrays"].items():
            print(f"  {k}: {v}")
    print()
    print(format_description())


if __name__ == "__main__":
    main()

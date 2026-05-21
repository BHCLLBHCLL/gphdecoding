#!/usr/bin/env python3
"""
GPH (Geometry/Polyhedron) Binary Format Parser and Reverse-Engineered Format Describer.

Parses any .gph file and outputs a structured format description.
Byte order: Big-Endian for all multi-byte integers and metadata.
"""

import sys
from pathlib import Path

from gph_model import (
    find_section,
    iter_data_blocks,
    open_gph_buffer,
    parse_ls_assemblies_summary,
    parse_ls_links_summary,
    parse_ls_nodes_vertices,
    parse_ls_parts,
    parse_ls_string_list,
    parse_ls_surface_regions_summary,
    read_i32_be,
    section_end,
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

    parts = parse_ls_parts(data)
    if parts:
        result["partition"]["LS_Parts"] = parts

    regions = parse_ls_string_list(data, "LS_VolumeRegions")
    if regions:
        result["partition"]["LS_VolumeRegions"] = regions

    surf = parse_ls_surface_regions_summary(data)
    if surf:
        result["partition"]["LS_SurfaceRegions"] = surf

    asm = parse_ls_assemblies_summary(data)
    if asm.get("part_paths") or asm.get("root_empty_prefix"):
        result["partition"]["LS_Assemblies"] = asm

    cv_sec = find_section(data, "LS_CvolIdOfElements")
    if cv_sec >= 0:
        sec_end = section_end(data, cv_sec)
        for p, bc in iter_data_blocks(data, cv_sec, sec_end):
            if bc % 4 == 0 and bc >= 4:
                n_cv = bc // 4
                sample = [read_i32_be(data, p + i * 4) for i in range(min(10, n_cv))]
                result["data_arrays"]["LS_CvolIdOfElements_count"] = n_cv
                result["data_arrays"]["LS_CvolIdOfElements_sample"] = sample
                import numpy as np
                cv_arr = np.frombuffer(data, dtype=">i4", count=n_cv, offset=p)
                result["data_arrays"]["LS_CvolIdOfElements_unique"] = sorted(
                    {int(x) for x in np.unique(cv_arr[:min(n_cv, 1_000_000)])}
                )[:20]
                break

    return result


def _section_blurb(name: str) -> str:
    blurbs = {
        "LS_CvolIdOfElements": "I4[n_cells] per-cell cvol_id (opaque Part id, not list index)",
        "LS_Links": "face topology: owner, neighbor, npe, conn (CSR; may split >1 GiB)",
        "LS_Nodes": "R8[n,3] vertex coords (three float64 axis blocks)",
        "LS_SurfaceRegions": "named BC regions -> global face id lists",
        "LS_VolumeRegions": "volume region names (-> CGNS zones)",
        "LS_Parts": "part names + cvol_id descriptors",
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
  Three equal-sized float64 axis blocks (X, Y, Z file order).
  Dialect auto-detection: standard big-endian float64 vs word-reversed
  (legacy); word-reversed files use X,Z,Y on disk -> permuted to X,Y,Z.

5. LS_Links - face / cell topology
----------------------------------
  At least three equal I4 blocks of length n_faces:
    owner, neighbor (0xFFFFFFFF = boundary), npe (nodes per face).
  Optional 4-byte face_type block (legacy SCTpre).
  conn block: sum(npe) vertex indices, 0-based, CSR layout:
    face i nodes = conn[face_offsets[i]:face_offsets[i+1]].

  Pure-triangle meshes: all npe=3.  Mixed/polyhedral (e.g. tr03.gph):
  npe varies (3..11+).  Voxel meshes (laptop_simplified_voxel_v4.gph) use
  npe in {4,5,6,7} with ~10^7 faces and ~10^7 cells.

  Large conn arrays (>~1 GiB): when sum(npe)*4 exceeds a single payload
  limit (~1073741824 bytes), the conn data is split:
    [12, bc1][conn part 1][bc1]  then  [I4=bc2][conn part 2 raw...]
  The continuation uses a bare I4 byte_count (no [I4=12] header tag).
  gph2cgns / gph_model concatenate both parts before CSR indexing.

  Files >512 MiB: gph2cgns, gph_parser and gphviewer memory-map the file
  instead of loading it entirely into RAM.

9. LS_CvolIdOfElements - per-cell partition id
----------------------------------------------
  I4[n_cells]: each cell's cvol_id matches the opaque id stored in the
  corresponding LS_Parts descriptor - NOT the 1-based index in LS_Parts.
  Example: laptop_simplified_voxel_less.gph uses cvol_ids {1, 9, 11}.

10. LS_Parts - part definitions
-------------------------------
  Repeated records: 255-byte ASCII name block + trailing descriptors.
  The cvol_id is the d0 field of the last [12, 4, cvol_id, 4] descriptor
  before the next part's name block.

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

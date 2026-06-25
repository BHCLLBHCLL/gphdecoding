# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

| Purpose | Command |
|---------|---------|
| Parse/inspect GPH format | `python3 gph_parser.py [file.gph]` |
| Convert GPH → CGNS/HDF5 | `python3 gph2cgns.py input.gph -o output.cgns` |
| Convert FPH → CGNS (with fields) | `python3 fph2cgns.py input.fph -o output.cgns` |
| Single-zone export | `python3 gph2cgns.py input.gph --single-zone -z ZoneName` |
| GPH Viewer GUI | `python3 gphviewer.py [file.gph]` |
| Headless viewer (CI/server) | `QT_QPA_PLATFORM=offscreen python3 gphviewer.py file.gph` |
| Coord dialect regression | `python3 tests/test_coord_score.py -v` |
| Lint | `~/.local/bin/ruff check *.py` |
| Syntax check | `python3 -m py_compile <file>.py` |

**Partition test:** `python3 tests/test_volume_zone_cells.py` compares volume-zone cell counts against `tests/*_orig.cgns`. Use `-v` to inspect LS_Parts chains and cvol_id mapping.

For full CGNS output validation, run `gph2cgns.py` on `tests/box_ansa.gph` and compare with `tests/box_ansa_orig.cgns`.

No build step — each script is a standalone entry-point run directly with `python3`.

## Architecture

This is a **flat Python CLI/GUI toolkit** for reverse-engineering and converting `.gph` files (proprietary CFD mesh format from Software Cradle scFLOW / SCTpre, magic `CRDL-FLD`, all integers/floats big-endian) into CGNS/HDF5. No web server, database, or services.

### File roles

- **`gph_model.py`** — Shared library: data types (`GphNode`, `GphDocument`), binary section scanners, heuristics for vertex-coordinate dialect detection, LS_Links/LS_Parts/LS_Assemblies parsers, mesh-preview builder. **Single source of truth** for GPH parsing logic, including `parse_ls_parts()` / cvol_id resolution.
- **`gph2cgns.py`** — Mesh-only converter: imports `parse_ls_nodes_xyz`, `_read_conn_continuations`, and `parse_ls_parts` from `gph_model`. Reads GPH via mmap for files >512 MiB. Writes HDF5 v0 superblock CGNS (libver earliest/v108) for ANSA compatibility.
- **`fph2cgns.py`** — FPH converter (mesh + FlowSolution field data): same mesh parsing as `gph2cgns` via `parse_ls_nodes_xyz`.
- **`gph_parser.py`** — CLI inspector: imports all parsing functions from `gph_model`, prints structured section layout and format description.
- **`gphviewer.py`** — PyQt GUI browser/editor: imports parsing functions from `gph_model`, builds tree view + hex dump + data table + 3D mesh preview via `MeshPreviewWidget`.

### Dependencies

`numpy`, `h5py` (converter only), `PyQt6` (viewer; falls back to `PyQt5` on Python < 3.12). System library `libegl1` required for PyQt6 rendering in headless environments.

## Key technical details

### GPH section discovery

Section offsets vary per file. All tools locate sections dynamically by scanning the buffer for 32-byte ASCII labels preceded by `[I4=32]`. The list of known section names is duplicated in `gph_model._SECTION_BOUNDARY_NAMES`, `gph2cgns._section_end`, and `gph_parser._parse_gph_buffer`. All three must stay in sync.

### Vertex coordinate dialect auto-detection

`LS_Nodes` contains three equal-sized axis blocks (X, Y, Z file order). Encodings:

1. **Big-endian float32** (FPH solver-result files, e.g. `tests/tr03_9.fph` — descriptor type 4)
2. **Standard big-endian float64** (modern ANSA/scFLOW GPH exports)
3. **Word-reversed float64** (legacy — each 8-byte double stored as `[low32_BE][high32_BE]`, disk axes X, Z, Y)

**Canonical parser**: `gph_model.parse_ls_nodes_xyz()` (used by `gph2cgns`, `fph2cgns`, `gph_parser`, `gphviewer`). It scores candidate decodings with `_score_coord_axes` using `_COORD_MIN_ABSMAX = 1e-4` m so float32 misread as float64 (~1e-13 denormals) loses. Vertex count comes from type descriptors (`dim0`), not `byte_count // 8`.

### LS_Links conn splitting (files > ~1 GiB connectivity)

When `sum(npe)*4` exceeds ~1 GiB, the connectivity array is split into multiple segments. The first segment uses the standard `[12, bc][payload][bc]` block format; continuation chunks use bare `[I4=byte_count][payload]` (no `[I4=12]` header). `_read_conn_continuations` in `gph_model` handles concatenation. The converter imports this directly.

### cvol_id resolution in LS_Parts

The `cvol_id` for each part is NOT the 1-based index. `parse_ls_parts` returns `PartCvolSpec = int | frozenset[int]`:

- **Simple part**: post-name chain `[1, cvol_id]` → single `int` (e.g. `outlet11` → 2).
- **Composite part**: `[12,4,N,4]` + `I4[N]` membership list → `frozenset` (e.g. `air_domain` in `laptop_simplified_more_regions.gph` owns 66 ids, 151375 cells). **Do not** treat `N` as a cvol_id.

**Canonical implementation**: `gph_model.parse_ls_parts(data, cvol_id=…)` (imported by `gph2cgns` as `_parse_ls_parts_with_cvol_ids`). Zone masks use `part_cvol_cell_mask(cvol_id, spec)`.

1. Parse `LS_CvolIdOfElements` first; unique values form set `S`.
2. Detect composite layout (`_parse_part_cvol_membership`) before single-id chain rule.
3. Simple parts: last chain value in `S`.
4. Always pass the per-cell cvol array when available.

Test composite cases: `python tests/test_volume_zone_cells.py -v tests/laptop_simplified_more_regions.gph`

### CGNS export layout (matching FLDUTIL)

- **NGON_n** for faces, **NFACE_n** for cells (signed face IDs: positive = owner, negative = neighbor)
- Multi-zone layout inferred from `LS_VolumeRegions` + `LS_Parts` + `LS_Assemblies` XML
- Surface regions map to `ZoneBC_t` with one `BC_t` per region; empty zones get BC nodes with no `data` dataset (matching the vendor exporter)
- HDF5 v0 superblock with `track_order=False` groups for ANSA compatibility
- Vertex renumbering by first-use order in face connectivity (matches FLDUTIL)
- Zone naming: paths with depth ≥ 2 use the dotted path verbatim; shallower paths get an `FPHPARTS.` prefix

### Memory-mapping

Files larger than 512 MiB are memory-mapped (`mmap`) by all three tools rather than loaded into RAM. This makes multi-gigabyte meshes (~8 GiB `laptop_simplified_denser.gph`) inspectable. Large-file mmap mode is **read-only** — the viewer disables Save/Save As for mmap'd files.

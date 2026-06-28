# AGENTS.md

## Cursor Cloud specific instructions

This is a flat Python CLI/GUI toolkit for converting `.gph` (proprietary CFD mesh format from SCTpre) to CGNS/HDF5. No web server, database, or long-running service is needed.

### Services

| Tool | Command | Notes |
|------|---------|-------|
| GPH Parser | `python3 gph_parser.py tests/box_ansa.gph` | Outputs section layout and format description |
| GPH→CGNS Converter | `python3 gph2cgns.py tests/box_ansa.gph -o output.cgns` | Mesh-only CGNS; supports float32/float64 LS_Nodes |
| FPH→CGNS Converter | `python3 fph2cgns.py tests/tr03_9.fph -o output.cgns` | FlowSolution default R4; `--flow-f64`, `--skip-fluid-region`, `--clip-flow 1` |
| GPH Viewer (GUI) | `QT_QPA_PLATFORM=offscreen python3 gphviewer.py tests/box_ansa.gph` | PyQt6 GUI; use `QT_QPA_PLATFORM=offscreen` in headless environments |
| Zone cell test | `python3 tests/test_volume_zone_cells.py -v` | Compare zone plan vs `*_orig.cgns`; `-v` shows simple/composite Part cvol specs |
| Coord dialect test | `python3 tests/test_coord_score.py -v` | LS_Nodes float32/float64 scoring; includes `tests/tr03_9.fph` |

### Lint / Test / Build

- **Lint**: `~/.local/bin/ruff check *.py` (minor pre-existing warnings exist in unused imports/variables)
- **Syntax check**: `python3 -m py_compile <file>.py`
- **Partition / zone test**: `python3 tests/test_volume_zone_cells.py` — compares each `tests/*.gph` against `{stem}_orig.cgns` (volume-zone cell counts). Add `-v` for LS_Parts descriptor chains, cvol_id histogram, and zone-selection notes.
- **Full conversion check**: run `gph2cgns.py` on a sample and compare with reference CGNS (e.g. `tests/box_ansa_orig.cgns`).
- **No build step** — scripts are run directly with `python3`.

### Dependencies

Installed via `pip install numpy h5py PyQt6`. The `requirements.txt` lists PyQt5 for Python < 3.12; on Python 3.12+ use PyQt6 instead. System library `libegl1` is required for PyQt6 rendering.

### Gotchas

- The viewer (`gphviewer.py`) auto-detects PyQt6 first, then falls back to PyQt5. On Python 3.12+ only PyQt6 works.
- `ruff` is installed to `~/.local/bin/` which may not be on PATH; use full path or add it to PATH.
- The sample meshes for regression live under **`tests/`** (`box_ansa.gph`, `laptop_simplified_voxel_less.gph`, matching `*_orig.cgns`). Larger samples (`tr03.gph`, `laptop_simplified*.gph`) remain in the repo root.
- **LS_Nodes parsing** is centralized in `gph_model.parse_ls_nodes_xyz()` (float32 / float64 BE / word-reversed). Do not duplicate dialect logic in converters or the viewer.
- **FPH conversion** uses `fph2cgns.py` only. FlowSolution defaults to **R4 (float32)**; `--flow-f64` for R8. `--skip-fluid-region` drops the entire `FluidRegion` Zone_t from output (not just field data). Zone filtering: `_filter_zone_plan()` in `fph2cgns.py`.
- **LS_Parts cvol_id mapping** lives in `gph_model.parse_ls_parts(data, cvol_id=…)` → `PartCvolSpec` (`int` or `frozenset` for composite parts like `air_domain`). Zone masks: `part_cvol_cell_mask()`. `gph2cgns` imports this — do not duplicate. Multi-region laptop files (`laptop_simplified_more_regions.gph`) use `[12,4,N,4]+I4[N]` membership lists; **N is a count, not a cvol_id**.

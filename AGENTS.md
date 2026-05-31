# AGENTS.md

## Cursor Cloud specific instructions

This is a flat Python CLI/GUI toolkit for converting `.gph` (proprietary CFD mesh format from SCTpre) to CGNS/HDF5. No web server, database, or long-running service is needed.

### Services

| Tool | Command | Notes |
|------|---------|-------|
| GPH Parser | `python3 gph_parser.py tests/box_ansa.gph` | Outputs section layout and format description |
| GPH→CGNS Converter | `python3 gph2cgns.py tests/box_ansa.gph -o output.cgns` | Core functionality; requires numpy + h5py |
| GPH Viewer (GUI) | `QT_QPA_PLATFORM=offscreen python3 gphviewer.py tests/box_ansa.gph` | PyQt6 GUI; use `QT_QPA_PLATFORM=offscreen` in headless environments |
| Zone cell test | `python3 tests/test_volume_zone_cells.py -v` | Compare GPH zone plan vs `*_orig.cgns`; `-v` prints LS_Parts / cvol_id details |

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
- **LS_Parts cvol_id mapping** lives in `gph_model.parse_ls_parts(data, cvol_id=…)`; always pass the `LS_CvolIdOfElements` array when available. `gph2cgns` imports this function — do not reintroduce a duplicate copy.

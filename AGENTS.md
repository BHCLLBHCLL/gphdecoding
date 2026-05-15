# AGENTS.md

## Cursor Cloud specific instructions

This is a flat Python CLI/GUI toolkit for converting `.gph` (proprietary CFD mesh format from SCTpre) to CGNS/HDF5. No web server, database, or long-running service is needed.

### Services

| Tool | Command | Notes |
|------|---------|-------|
| GPH Parser | `python3 gph_parser.py box.gph` | Outputs section layout and format description |
| GPH→CGNS Converter | `python3 gph2cgns.py box.gph -o output.cgns` | Core functionality; requires numpy + h5py |
| GPH Viewer (GUI) | `QT_QPA_PLATFORM=offscreen python3 gphviewer.py box.gph` | PyQt6 GUI; use `QT_QPA_PLATFORM=offscreen` in headless environments |

### Lint / Test / Build

- **Lint**: `~/.local/bin/ruff check *.py` (minor pre-existing warnings exist in unused imports/variables)
- **Syntax check**: `python3 -m py_compile <file>.py`
- **No automated test suite exists.** Validate by running the converter and comparing output with the reference file `box_ngons.cgns`.
- **No build step** — scripts are run directly with `python3`.

### Dependencies

Installed via `pip install numpy h5py PyQt6`. The `requirements.txt` lists PyQt5 for Python < 3.12; on Python 3.12+ use PyQt6 instead. System library `libegl1` is required for PyQt6 rendering.

### Gotchas

- The viewer (`gphviewer.py`) auto-detects PyQt6 first, then falls back to PyQt5. On Python 3.12+ only PyQt6 works.
- `ruff` is installed to `~/.local/bin/` which may not be on PATH; use full path or add it to PATH.
- The sample data file `box.gph` is committed in the repo root for testing.

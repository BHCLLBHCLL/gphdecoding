# gphdecoding

Reverse-engineer and convert **GPH** (Software Cradle scFLOW / SCTpre CFD mesh) to **CGNS/HDF5**.

## Files

| File | Role |
|------|------|
| **gph_model.py** | Shared parsing library (`parse_ls_parts`, `parse_ls_cvol_ids`, mesh topology, tree model) |
| **gph_parser.py** | CLI inspector — section layout, format description, partition summary |
| **gph2cgns.py** | GPH → CGNS/HDF5 converter (multi-zone, polyhedral NGON_n / NFACE_n) |
| **gphviewer.py** | PyQt GUI browser/editor (tree + hex dump + 3D preview) |
| **GPH_FORMAT_SPEC.md** | Reverse-engineered GPH binary format specification (Chinese) |
| **DEV_SUMMARY.md** | Development log and design decisions |
| **tests/test_volume_zone_cells.py** | Partition regression test (GPH zone plan vs `*_orig.cgns`) |

Sample meshes: `tests/box_ansa.gph`, `tests/laptop_simplified_voxel_less.gph` (with `*_orig.cgns` references); larger files `tr03.gph`, `laptop_simplified*.gph` in repo root.

## Usage

### Parse GPH format
```bash
python gph_parser.py [file.gph]
```
Outputs section layout, data samples, and full format description.

### Convert GPH to CGNS
```bash
python gph2cgns.py tests/box_ansa.gph -o box.cgns
python gph2cgns.py input.gph -o output.cgns -z ZoneName   # single-zone export
```
Requires: `pip install numpy h5py`

### Test volume-zone cell counts
```bash
python tests/test_volume_zone_cells.py          # all tests/*.gph vs *_orig.cgns
python tests/test_volume_zone_cells.py -v       # + LS_Parts chains, cvol_id details
python tests/test_volume_zone_cells.py -v tests/box_ansa.gph
```
Compares each GPH file's zone plan against a matching `{stem}_orig.cgns` reference when present.

### GPH Viewer (GUI)
```bash
python gphviewer.py [file.gph]
```
Tree + hex/table view + mesh preview. Requires PyQt6 (Python 3.12+) or PyQt5.

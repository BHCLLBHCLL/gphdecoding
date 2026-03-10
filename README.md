# gphdecoding

Figure out the data structure of the GPH file.

## Files

- **gph_parser.py** - Python script to parse `box.gph` and output format description
- **gph2cgns.py** - Convert GPH to CGNS/HDF5 mesh (vertices only; elements TBD)
- **gphviewer.py** - PyQt GUI to browse/edit GPH (like HDFView)
- **GPH_FORMAT_SPEC.md** - Reverse-engineered GPH binary format specification (Chinese)

## Usage

### Parse GPH format
```bash
python gph_parser.py [box.gph]
```
Outputs section layout, data samples, and full format description.

### Convert GPH to CGNS
```bash
python gph2cgns.py box.gph -o box.cgns
# or
python gph2cgns.py input.gph -o output.cgns -z Zone1
```
Requires: `pip install numpy h5py`

### GPH Viewer (GUI)
```bash
python gphviewer.py [box.gph]
```
Opens a PyQt window to browse the GPH structure (tree + hex/table view) and edit scalars or vertex coordinates, then save. Requires: `pip install PyQt5`

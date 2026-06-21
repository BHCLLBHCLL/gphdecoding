import h5py

for p in ["tests/tr03_9_converted2.cgns", "tests/box_ansa_orig.cgns", "tr03_orig.cgns"]:
    with h5py.File(p, "r") as f:
        fmt = bytes(f[" format"][...]).decode("ascii", errors="replace").rstrip("\x00")
        ver = bytes(f[" hdf5version"][...]).decode("ascii", errors="replace").rstrip("\x00")
        print(f"{p}:")
        print(f"  format = {fmt!r}")
        print(f"  hdf5version = {ver!r}")
        # Check raw bytes
        print(f"  format raw = {list(f[' format'][...])}")
        print(f"  hdf5version raw = {list(f[' hdf5version'][...])}")

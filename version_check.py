"""Check CGNS/HDF5 version info in reference and converted files."""
import h5py
import numpy as np


def check_file(path, out):
    out.append(f"\n=== {path} ===")
    with h5py.File(path, "r") as f:
        # Root attributes
        out.append("Root attributes:")
        for k in sorted(f.attrs.keys()):
            v = f.attrs[k]
            if isinstance(v, bytes):
                v = v.decode("ascii", errors="replace")
            if isinstance(v, np.ndarray):
                if v.dtype == np.int8:
                    v = bytes(v).decode("ascii", errors="replace").rstrip("\x00")
            out.append(f"  {k} = {v!r}")

        # CGNSLibraryVersion
        if "CGNSLibraryVersion" in f:
            v = f["CGNSLibraryVersion"]
            if isinstance(v, h5py.Group):
                out.append(f"CGNSLibraryVersion: Group, keys={list(v.keys())}")
                for k in v.keys():
                    item = v[k]
                    if isinstance(item, h5py.Dataset):
                        arr = item[...]
                        out.append(f"  {k} = {arr!r} (dtype={arr.dtype})")
                attrs = dict(v.attrs)
                for k in sorted(attrs.keys()):
                    av = attrs[k]
                    if isinstance(av, bytes):
                        av = av.decode("ascii", errors="replace")
                    out.append(f"  attr {k} = {av!r}")
            else:
                arr = v[...]
                out.append(f"CGNSLibraryVersion: {arr!r} (dtype={arr.dtype})")
        else:
            out.append("CGNSLibraryVersion: NOT FOUND")

        # HDF5 file version (from libver)
        out.append(f"HDF5 libver (read): {f.libver}")
        out.append(f"HDF5 driver: {f.driver}")

        # Base attributes
        if "Base" in f:
            base = f["Base"]
            out.append("Base attributes:")
            for k in sorted(base.attrs.keys()):
                v = base.attrs[k]
                if isinstance(v, bytes):
                    v = v.decode("ascii", errors="replace")
                out.append(f"  {k} = {v!r}")


def main():
    out = []
    out.append("Checking CGNS/HDF5 version info:")

    import os
    for p in ["tests/tr03_9_converted2.cgns", "tests/box_ansa_orig.cgns",
              "tr03_orig.cgns", "tests/box_ansa_fph_test.cgns"]:
        if os.path.exists(p):
            check_file(p, out)
        else:
            out.append(f"\n=== {p} === NOT FOUND")

    # Check HDF5 superblock version using h5py low-level API
    out.append("\n=== HDF5 superblock info ===")
    for p in ["tests/tr03_9_converted2.cgns", "tests/box_ansa_orig.cgns", "tr03_orig.cgns"]:
        if os.path.exists(p):
            with open(p, "rb") as fp:
                header = fp.read(16)
                out.append(f"{p}: first 16 bytes = {header.hex()}")
                # HDF5 signature is at offset 0: \x89HDF\r\n\x1a\n
                if header[:8] == b'\x89HDF\r\n\x1a\n':
                    out.append(f"  Valid HDF5 signature")
                    # Superblock version is at offset 8
                    out.append(f"  Superblock version: {header[8]}")
                else:
                    out.append(f"  NOT a valid HDF5 file!")

    with open("version_check.txt", "w", encoding="utf-8") as fp:
        fp.write("\n".join(out))
    print("Wrote version_check.txt")


if __name__ == "__main__":
    main()

"""Check dataset storage properties."""
import h5py


def check_datasets(path, out):
    out.append(f"\n=== {path} ===")
    with h5py.File(path, "r") as f:
        base = f["Base"]
        for zname in sorted(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel != "Zone_t":
                continue
            zone = base[zname]
            out.append(f"Zone: {zname}")

            # Check GridCoordinates
            if "GridCoordinates" in zone:
                gc = zone["GridCoordinates"]
                for cname in ["CoordinateX", "CoordinateY", "CoordinateZ"]:
                    if cname in gc:
                        ds = gc[cname][" data"]
                        out.append(f"  {cname}: shape={ds.shape} dtype={ds.dtype} "
                                  f"chunks={ds.chunks} compression={ds.compression} "
                                  f"maxshape={ds.maxshape}")

            # Check Elements
            for cname in zone.keys():
                citem = zone[cname]
                if not isinstance(citem, h5py.Group):
                    continue
                cattrs = dict(citem.attrs)
                clabel = cattrs.get("label", b"")
                if isinstance(clabel, bytes):
                    clabel = clabel.decode("ascii", errors="replace")
                if clabel == "Elements_t":
                    for subname in citem.keys():
                        sitem = citem[subname]
                        if isinstance(sitem, h5py.Dataset):
                            out.append(f"  {cname}/{subname}: shape={sitem.shape} "
                                      f"dtype={sitem.dtype} chunks={sitem.chunks} "
                                      f"compression={sitem.compression}")

            # Check FlowSolution
            if "FlowSolution" in zone:
                fs = zone["FlowSolution"]
                for vname in fs.keys():
                    vitem = fs[vname]
                    if isinstance(vitem, h5py.Group):
                        for subname in vitem.keys():
                            sitem = vitem[subname]
                            if isinstance(sitem, h5py.Dataset):
                                out.append(f"  FlowSolution/{vname}/{subname}: "
                                          f"shape={sitem.shape} dtype={sitem.dtype} "
                                          f"chunks={sitem.chunks} compression={sitem.compression}")
            break  # only first zone


def main():
    out = []
    import os
    for p in ["tests/tr03_9_converted2.cgns", "tests/box_ansa_orig.cgns",
              "tr03_orig.cgns", "tests/box_ansa_fph_test.cgns"]:
        if os.path.exists(p):
            check_datasets(p, out)

    with open("dataset_storage.txt", "w", encoding="utf-8") as fp:
        fp.write("\n".join(out))
    print("Wrote dataset_storage.txt")


if __name__ == "__main__":
    main()

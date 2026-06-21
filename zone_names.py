"""Compare zone names and check for special characters."""
import h5py


def list_zones(path, out):
    out.append(f"\n=== {path} ===")
    with h5py.File(path, "r") as f:
        if "Base" not in f:
            out.append("  No Base found")
            return
        base = f["Base"]
        for zname in sorted(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel == "Zone_t":
                zdata = base[zname][" data"][...]
                out.append(f"  Zone: '{zname}' (repr={repr(zname)}) "
                           f"vertices={zdata[0][0]} cells={zdata[1][0]}")
                # Check for special characters
                special = [c for c in zname if not (c.isalnum() or c in "._-")]
                if special:
                    out.append(f"    SPECIAL CHARS: {special}")


def main():
    out = []
    import os
    for p in ["tests/tr03_9_converted2.cgns", "tests/box_ansa_orig.cgns",
              "tr03_orig.cgns", "tests/box_ansa_fph_test.cgns"]:
        if os.path.exists(p):
            list_zones(p, out)
        else:
            out.append(f"\n=== {p} === NOT FOUND")

    with open("zone_names.txt", "w", encoding="utf-8") as fp:
        fp.write("\n".join(out))
    print("Wrote zone_names.txt")


if __name__ == "__main__":
    main()

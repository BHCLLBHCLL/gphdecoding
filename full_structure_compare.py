"""Detailed structural comparison between reference and converted CGNS."""
import h5py
import numpy as np


def dump_full_tree(g, path="", depth=0, max_depth=10, out=None):
    if out is None:
        out = []
    if depth > max_depth:
        return out
    for name in sorted(g.keys()):
        item = g[name]
        indent = "  " * depth
        if isinstance(item, h5py.Group):
            attrs = dict(item.attrs)
            label = attrs.get("label", b"")
            if isinstance(label, bytes):
                label = label.decode("ascii", errors="replace")
            attr_strs = []
            for k in sorted(attrs.keys()):
                if k == "label":
                    continue
                v = attrs[k]
                if isinstance(v, bytes):
                    v = v.decode("ascii", errors="replace")
                attr_strs.append(f"{k}={v}")
            attr_info = f" [{', '.join(attr_strs)}]" if attr_strs else ""
            out.append(f"{indent}{name}/ (label={label}{attr_info})")
            dump_full_tree(item, f"{path}/{name}", depth + 1, max_depth, out)
        elif isinstance(item, h5py.Dataset):
            attrs = dict(item.attrs)
            arr = item[...]
            # For large arrays, just show shape/dtype/range
            if arr.size > 20:
                if np.issubdtype(arr.dtype, np.number):
                    info = f"shape={arr.shape} dtype={arr.dtype} min={arr.min()} max={arr.max()}"
                else:
                    info = f"shape={arr.shape} dtype={arr.dtype}"
            else:
                info = f"value={arr.tolist()} dtype={arr.dtype}"
            out.append(f"{indent}{name} = ({info})")
    return out


def main():
    out = []

    # Dump first zone of reference in full detail
    out.append("=" * 80)
    out.append("REFERENCE: tests/tr03_9.cgns - first zone (full detail)")
    out.append("=" * 80)
    import os
    if os.path.exists("tests/tr03_9.cgns"):
        with h5py.File("tests/tr03_9.cgns", "r") as f:
            base = f["Base"]
            for zname in sorted(base.keys()):
                zattrs = dict(base[zname].attrs)
                zlabel = zattrs.get("label", b"")
                if isinstance(zlabel, bytes):
                    zlabel = zlabel.decode("ascii", errors="replace")
                if zlabel == "Zone_t":
                    out.append(f"Zone: {zname}")
                    lines = dump_full_tree(base[zname], max_depth=5)
                    out.extend(lines)
                    break
    else:
        out.append("File not found: tests/tr03_9.cgns")

    # Dump first zone of converted in full detail
    out.append("")
    out.append("=" * 80)
    out.append("CONVERTED: tests/tr03_9_converted2.cgns - first zone (full detail)")
    out.append("=" * 80)
    with h5py.File("tests/tr03_9_converted2.cgns", "r") as f:
        base = f["Base"]
        for zname in sorted(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel == "Zone_t":
                out.append(f"Zone: {zname}")
                lines = dump_full_tree(base[zname], max_depth=5)
                out.extend(lines)
                break

    # Also dump box_ansa (known working GPH conversion) for comparison
    out.append("")
    out.append("=" * 80)
    out.append("BOX_ANSA: tests/box_ansa_orig.cgns - first zone (full detail)")
    out.append("=" * 80)
    with h5py.File("tests/box_ansa_orig.cgns", "r") as f:
        base = f["Base"]
        for zname in sorted(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel == "Zone_t":
                out.append(f"Zone: {zname}")
                lines = dump_full_tree(base[zname], max_depth=5)
                out.extend(lines)
                break

    # Also dump fph2cgns output on box_ansa for comparison
    out.append("")
    out.append("=" * 80)
    out.append("BOX_ANSA converted by fph2cgns.py")
    out.append("=" * 80)
    import subprocess
    subprocess.run(["python", "fph2cgns.py", "tests/box_ansa.gph",
                    "-o", "tests/box_ansa_fph_test.cgns"],
                   capture_output=True)
    with h5py.File("tests/box_ansa_fph_test.cgns", "r") as f:
        base = f["Base"]
        for zname in sorted(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel == "Zone_t":
                out.append(f"Zone: {zname}")
                lines = dump_full_tree(base[zname], max_depth=5)
                out.extend(lines)
                break

    with open("full_structure_compare.txt", "w", encoding="utf-8") as fp:
        fp.write("\n".join(out))
    print("Wrote full_structure_compare.txt")


if __name__ == "__main__":
    main()

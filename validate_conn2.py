"""Validate connectivity and test without FlowSolution."""
import h5py
import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def validate_connectivity(path, out):
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
            out.append(f"\nZone: {zname}")

            # Find NGON_n and NFACE_n
            ngon_range = None
            nface_range = None
            ngon_conn = None
            nface_conn = None
            n_vertices = None

            zdata = zone[" data"][...]
            n_vertices = zdata[0][0]

            for cname in zone.keys():
                citem = zone[cname]
                if not isinstance(citem, h5py.Group):
                    continue
                cattrs = dict(citem.attrs)
                clabel = cattrs.get("label", b"")
                if isinstance(clabel, bytes):
                    clabel = clabel.decode("ascii", errors="replace")
                if clabel == "Elements_t":
                    data = citem[" data"][...]
                    etype = data[0]
                    erange = citem["ElementRange"][" data"][...]
                    conn = citem["ElementConnectivity"][" data"][...]
                    if etype == 22:  # NGON_n
                        ngon_range = (erange[0], erange[1])
                        ngon_conn = conn
                        out.append(f"  NGON_n: range={ngon_range}, conn shape={conn.shape}")
                        out.append(f"    conn min={conn.min()} max={conn.max()}")
                        out.append(f"    conn has_zero={np.any(conn == 0)}")
                        out.append(f"    conn has_negative={np.any(conn < 0)}")
                        # Check vertex indices are in [1, n_vertices]
                        pos_conn = conn[conn > 0]
                        if len(pos_conn) > 0:
                            out.append(f"    vertex idx range: [{pos_conn.min()}, {pos_conn.max()}]")
                            out.append(f"    n_vertices={n_vertices}")
                            if pos_conn.max() > n_vertices:
                                out.append(f"    ERROR: vertex index {pos_conn.max()} > n_vertices {n_vertices}")
                            if pos_conn.min() < 1:
                                out.append(f"    ERROR: vertex index {pos_conn.min()} < 1")
                    elif etype == 23:  # NFACE_n
                        nface_range = (erange[0], erange[1])
                        nface_conn = conn
                        out.append(f"  NFACE_n: range={nface_range}, conn shape={conn.shape}")
                        out.append(f"    conn min={conn.min()} max={conn.max()}")
                        out.append(f"    conn has_zero={np.any(conn == 0)}")
                        out.append(f"    conn has_negative={np.any(conn < 0)}")
                        # Check face indices are in NGON_n range
                        if ngon_range:
                            pos_conn = conn[conn > 0]
                            if len(pos_conn) > 0:
                                out.append(f"    face idx range: [{pos_conn.min()}, {pos_conn.max()}]")
                                if pos_conn.max() > ngon_range[1]:
                                    out.append(f"    ERROR: face index {pos_conn.max()} > ngon_max {ngon_range[1]}")
                                if pos_conn.min() < ngon_range[0]:
                                    out.append(f"    ERROR: face index {pos_conn.min()} < ngon_min {ngon_range[0]}")
            break  # only first zone


def main():
    out = []

    # Validate connectivity
    validate_connectivity("tests/tr03_9_converted2.cgns", out)
    validate_connectivity("tr03_orig.cgns", out)
    validate_connectivity("tests/box_ansa_orig.cgns", out)

    # Create a copy without FlowSolution
    out.append("\n" + "=" * 60)
    out.append("Creating copy without FlowSolution...")
    import shutil
    shutil.copy("tests/tr03_9_converted2.cgns", "tests/tr03_9_noflow.cgns")
    with h5py.File("tests/tr03_9_noflow.cgns", "r+") as f:
        base = f["Base"]
        for zname in list(base.keys()):
            zattrs = dict(base[zname].attrs)
            zlabel = zattrs.get("label", b"")
            if isinstance(zlabel, bytes):
                zlabel = zlabel.decode("ascii", errors="replace")
            if zlabel == "Zone_t":
                zone = base[zname]
                if "FlowSolution" in zone:
                    del zone["FlowSolution"]
                    out.append(f"  Removed FlowSolution from {zname}")
    out.append("Created tests/tr03_9_noflow.cgns (no FlowSolution)")

    with open("validate_conn.txt", "w", encoding="utf-8") as fp:
        fp.write("\n".join(out))
    print("Wrote validate_conn.txt")
    print("Created tests/tr03_9_noflow.cgns for testing")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Diagnose zone cell counts and surface region face ID ranges."""
import sys
from pathlib import Path

import numpy as np

from gph2cgns import _build_zone_plan, parse_gph_mesh
from gph_model import open_gph_buffer, parse_ls_cvol_ids, parse_ls_parts, parse_ls_links_summary, parse_ls_surface_regions


def diagnose(path: Path) -> None:
    print(f"=== {path.name} ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        links = parse_ls_links_summary(data)
        cvol = parse_ls_cvol_ids(data)
        parts = parse_ls_parts(data, cvol)
        regions = parse_ls_surface_regions(data)
        n_cells = links["n_cells"] if links else 0
        n_faces = links["n_faces"] if links else 0

        print(f"  n_cells={n_cells} n_faces={n_faces}", flush=True)
        if cvol is not None:
            print(f"  cvol_id: len={len(cvol)} unique={len(np.unique(cvol))} min={cvol.min()} max={cvol.max()}", flush=True)
            if len(cvol) != n_cells:
                print(f"  ** cvol_id length mismatch: {len(cvol)} != {n_cells} **", flush=True)
        else:
            print("  cvol_id: None", flush=True)

        print(f"  parts ({len(parts)}):", flush=True)
        for name, spec in parts:
            from gph_model import format_part_cvol_spec, part_cvol_cell_mask
            if cvol is not None and len(cvol) == n_cells:
                n = int(part_cvol_cell_mask(cvol, spec).sum())
            else:
                n = "?"
            print(f"    {name}: {format_part_cvol_spec(spec)} -> {n} cells", flush=True)

        if regions:
            all_ids = np.concatenate([f for _, f in regions]) if regions else np.array([])
            print(f"  surface_regions: {len(regions)} total face refs={all_ids.size}", flush=True)
            print(f"    face_id min={all_ids.min()} max={all_ids.max()}", flush=True)
            oob = (all_ids < 0) | (all_ids >= n_faces)
            print(f"    out of bounds [0,{n_faces}): {int(oob.sum())}", flush=True)
            if oob.any():
                bad = all_ids[oob]
                print(f"    bad samples: {np.unique(bad)[:20].tolist()}", flush=True)
            # 1-based hypothesis: if min==1 and max==n_faces, likely 1-based
            if all_ids.size and all_ids.min() >= 1 and all_ids.max() == n_faces:
                print("    ** likely 1-based face IDs (max==n_faces) **", flush=True)
            elif all_ids.size and all_ids.min() >= 1 and all_ids.max() == n_faces - 1:
                print("    looks 0-based", flush=True)

    mesh = parse_gph_mesh(str(path))
    plan = _build_zone_plan(mesh)
    counts = [int(m.sum()) for _, m in plan]
    print(f"  zone plan ({len(plan)} zones):", flush=True)
    for name, mask in plan:
        print(f"    {name}: {int(mask.sum())} cells", flush=True)
    if len(set(counts)) == 1 and len(plan) > 2:
        print("  ** all zones have same cell count **", flush=True)
    print(flush=True)


def main():
    root = Path(__file__).resolve().parent
    paths = [Path(p) for p in sys.argv[1:]] if len(sys.argv) > 1 else sorted(
        root.glob("*.gph"), key=lambda p: p.stat().st_size, reverse=True
    )[:8]
    for p in paths:
        if p.is_file():
            diagnose(p)


if __name__ == "__main__":
    main()

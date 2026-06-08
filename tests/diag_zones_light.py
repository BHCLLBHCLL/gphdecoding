#!/usr/bin/env python3
"""Lightweight diagnose: cvol_id, parts, surface face IDs (no full mesh parse)."""
import sys
from pathlib import Path

import numpy as np

from gph_model import (
    find_section,
    iter_data_blocks,
    open_gph_buffer,
    parse_ls_cvol_ids,
    parse_ls_links_summary,
    parse_ls_parts,
    parse_ls_string_list,
    parse_ls_surface_regions,
    section_end,
    format_part_cvol_spec,
    part_cvol_cell_mask,
)


def diagnose(path: Path) -> None:
    print(f"=== {path.name} ({path.stat().st_size // (1024*1024)} MiB) ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        links = parse_ls_links_summary(data)
        cvol = parse_ls_cvol_ids(data)
        parts = parse_ls_parts(data, cvol)
        vol_regions = parse_ls_string_list(data, "LS_VolumeRegions")
        surf = parse_ls_surface_regions(data)
        n_cells = links["n_cells"] if links else 0
        n_faces = links["n_faces"] if links else 0

        print(f"  n_cells={n_cells} n_faces={n_faces}", flush=True)
        print(f"  volume_regions ({len(vol_regions)}): {vol_regions[:6]}{'...' if len(vol_regions)>6 else ''}", flush=True)

        # LS_CvolIdOfElements blocks
        sec = find_section(data, "LS_CvolIdOfElements")
        if sec >= 0:
            sec_end = section_end(data, sec)
            cvol_blocks = [(bc, bc // 4) for _, bc in iter_data_blocks(data, sec, sec_end) if bc % 4 == 0]
            print(f"  LS_CvolIdOfElements blocks: {cvol_blocks[:8]}{'...' if len(cvol_blocks)>8 else ''}", flush=True)

        if cvol is not None:
            u = np.unique(cvol)
            print(f"  cvol_id: len={len(cvol)} match_n_cells={len(cvol)==n_cells} unique={len(u)} range=[{u.min()},{u.max()}]", flush=True)
            if len(cvol) != n_cells:
                print("  ** MISMATCH -> all zones fall back to full mesh **", flush=True)
        else:
            print("  cvol_id: None", flush=True)

        print(f"  parts ({len(parts)}):", flush=True)
        for name, spec in parts:
            n = int(part_cvol_cell_mask(cvol, spec).sum()) if cvol is not None and len(cvol) == n_cells else "?"
            print(f"    {name}: {format_part_cvol_spec(spec)} -> {n}", flush=True)

        if surf:
            all_ids = np.concatenate([f for _, f in surf])
            oob_lo = all_ids < 0
            oob_hi = all_ids >= n_faces
            oob_hi_1based = (all_ids >= 1) & (all_ids <= n_faces) & (all_ids > n_faces - 1)
            print(f"  surface_regions: {len(surf)} refs={all_ids.size} min={all_ids.min()} max={all_ids.max()}", flush=True)
            print(f"    OOB if 0-based: {int((oob_lo|oob_hi).sum())}  (max>=n_faces: {int(oob_hi.sum())})", flush=True)
            if oob_hi.any():
                print(f"    max OOB values: {all_ids[oob_hi][:5].tolist()}", flush=True)
            if all_ids.size and all_ids.min() >= 1 and int(all_ids.max()) == n_faces:
                print("    ** likely 1-based (max==n_faces) **", flush=True)
            # test 1-based fix
            if oob_hi.any() and not oob_lo.any() and all_ids.min() >= 1:
                adj = all_ids - 1
                print(f"    if subtract 1: OOB={int(((adj<0)|(adj>=n_faces)).sum())} max={adj.max()}", flush=True)
    print(flush=True)


def main():
    paths = [Path(p) for p in sys.argv[1:]]
    for p in paths:
        diagnose(p)


if __name__ == "__main__":
    main()

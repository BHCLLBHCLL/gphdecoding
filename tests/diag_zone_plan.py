#!/usr/bin/env python3
"""Simulate zone plan without full mesh parse."""
import sys
from pathlib import Path

import numpy as np

from gph2cgns import _build_zone_plan, _classify_zone_cells
from gph_model import (
    open_gph_buffer,
    parse_ls_assemblies_summary,
    parse_ls_cvol_ids,
    parse_ls_links_summary,
    parse_ls_parts,
    parse_ls_string_list,
)


def simulate(path: Path) -> None:
    print(f"=== {path.name} ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        links = parse_ls_links_summary(data)
        cvol = parse_ls_cvol_ids(data)
        parts = parse_ls_parts(data, cvol)
        regions = parse_ls_string_list(data, "LS_VolumeRegions")
        asm = parse_ls_assemblies_summary(data)
        n_cells = links["n_cells"]

        mesh = {
            "link_data": {"n_cells": n_cells},
            "cvol_id": cvol,
            "parts_with_cvol": parts,
            "volume_regions": regions,
            "assembly_info": asm,
        }
        plan = _build_zone_plan(mesh)
        counts = [int(m.sum()) for _, m in plan]
        uniq = set(counts)
        print(f"  n_cells={n_cells} cvol_ok={cvol is not None and len(cvol)==n_cells} parts={len(parts)}", flush=True)
        for name, mask in plan:
            print(f"    {name}: {int(mask.sum())}", flush=True)
        if len(uniq) == 1 and len(plan) > 2:
            print("  ** ALL SAME COUNT **", flush=True)
            if cvol is None or len(cvol) != n_cells:
                print("  cause: cvol_id missing or len mismatch", flush=True)
            elif not parts:
                print("  cause: no parts parsed", flush=True)
            else:
                # check name matching
                name_to_cvol = dict(parts)
                for zname, _ in plan:
                    if zname == "FluidRegion":
                        continue
                    m = _classify_zone_cells(zname, parts, cvol, n_cells)
                    if m.all():
                        print(f"  classify fallback (all cells): {zname}", flush=True)
    print(flush=True)


if __name__ == "__main__":
    for p in map(Path, sys.argv[1:]):
        simulate(p)

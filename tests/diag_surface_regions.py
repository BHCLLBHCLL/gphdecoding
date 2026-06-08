#!/usr/bin/env python3
"""Inspect LS_SurfaceRegions block layout."""
import sys
from pathlib import Path

import numpy as np

from gph_model import find_section, iter_data_blocks, open_gph_buffer, parse_ls_links_summary, read_i32_be, section_end


def inspect(path: Path) -> None:
    print(f"=== {path.name} ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        links = parse_ls_links_summary(data)
        n_faces = links["n_faces"] if links else 0
        sec = find_section(data, "LS_SurfaceRegions")
        sec_end = section_end(data, sec)
        blocks = list(iter_data_blocks(data, sec, sec_end))
        print(f"  n_faces={n_faces} section_len={sec_end-sec} blocks={len(blocks)}", flush=True)
        sizes = [bc for _, bc in blocks]
        from collections import Counter
        print(f"  block size histogram (top 8): {Counter(sizes).most_common(8)}", flush=True)
        # classify blocks
        i = 0
        region_idx = 0
        while i < len(blocks):
            p, bc = blocks[i]
            if bc < 256:
                raw = bytes(data[p:p+bc])
                if all(b == 0 or 32 <= b < 127 for b in raw):
                    name = raw.decode("ascii", errors="replace").strip("\x00").rstrip()
                    if i + 2 < len(blocks):
                        _, bc_i = blocks[i+1]
                        _, bc_w = blocks[i+2]
                        if bc_i == bc_w and bc_i % 4 == 0:
                            ids = np.frombuffer(data, dtype=">i4", count=bc_i//4, offset=blocks[i+1][0])
                            print(f"  region[{region_idx}] {name!r}: n={ids.size} min={ids.min()} max={ids.max()}", flush=True)
                            oob = (ids < 0) | (ids >= n_faces)
                            if oob.any():
                                print(f"    OOB={oob.sum()} samples={np.unique(ids[oob])[:5]}", flush=True)
                            region_idx += 1
                            i += 3
                            continue
            i += 1
    print(flush=True)


if __name__ == "__main__":
    for p in map(Path, sys.argv[1:]):
        inspect(p)

#!/usr/bin/env python3
"""Diagnose LS_Links conn block layout and continuation parsing."""
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from gph_model import (
    _CONN_CHUNK_BYTES,
    _conn_payload_size,
    _read_conn_continuations,
    find_section,
    iter_data_blocks,
    open_gph_buffer,
    parse_ls_links_summary,
    read_i32_be,
    section_end,
)


def analyze(path: Path) -> None:
    print(f"=== {path.name} size={path.stat().st_size} ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        summary = parse_ls_links_summary(data)
        if summary:
            for k, v in summary.items():
                print(f"  {k}: {v}", flush=True)

        sec = find_section(data, "LS_Links")
        sec_end = section_end(data, sec)
        print(f"  LS_Links section: offset={sec} end={sec_end} len={sec_end - sec}", flush=True)

        blocks = [(p, bc) for p, bc in iter_data_blocks(data, sec, sec_end) if bc > 0]
        print(f"  data blocks: {len(blocks)}", flush=True)

        block_sizes = [bc for _, bc in blocks]
        n_faces_block_size = None
        for size, count in Counter(block_sizes).most_common():
            if count >= 3 and size % 4 == 0 and size >= 4:
                n_faces_block_size = size
                break
        if n_faces_block_size is None:
            print("  ERROR: no triple block size", flush=True)
            return

        n_faces = n_faces_block_size // 4
        triples = [b for b in blocks if b[1] == n_faces_block_size][:3]
        npe_p = triples[2][0]
        npe = np.frombuffer(data, dtype=">u4", count=n_faces, offset=npe_p).astype(np.int64)
        expected = int(npe.sum())

        got, split, chunks = _conn_payload_size(data, blocks, triples, expected, sec_end)
        print(
            f"  conn expected={expected} got={got} split={split} chunks={chunks} "
            f"complete={got >= expected}",
            flush=True,
        )

        # List all non-triple I4 blocks >= 12 bytes, sorted by size
        others = [(p, bc) for p, bc in blocks if (p, bc) not in triples and bc >= 12 and bc % 4 == 0]
        others.sort(key=lambda x: -x[1])
        print("  non-triple I4 blocks (top 15):", flush=True)
        for p, bc in others[:15]:
            rel = p - sec
            entries = bc // 4
            tag = ""
            if bc == _CONN_CHUNK_BYTES:
                tag = " [1GiB]"
            elif entries == expected:
                tag = " [exact sum(npe)]"
            print(f"    bc={bc:12d} entries={entries:12d} payload_off={rel}{tag}", flush=True)

        # Walk raw bytes after primary conn block to show continuation layout
        conn_block = None
        for p, bc in blocks:
            if (p, bc) in triples:
                continue
            if bc % 4 != 0 or bc < 12:
                continue
            if bc // 4 == expected:
                conn_block = (p, bc)
                break
        if conn_block is None:
            for p, bc in blocks:
                if (p, bc) in triples:
                    continue
                if bc % 4 != 0 or bc < 12:
                    continue
                if conn_block is None or bc > conn_block[1]:
                    conn_block = (p, bc)
        if conn_block is None:
            print("  ERROR: no conn block", flush=True)
            return

        conn_p, conn_bc = conn_block
        pos = conn_p + conn_bc + 4
        print(f"  primary conn: bc={conn_bc} entries={conn_bc // 4} pos_after_trailer={pos - sec}", flush=True)

        chunk_idx = 0
        got_walk = conn_bc // 4
        while got_walk < expected and pos + 4 <= sec_end:
            bare_bc = read_i32_be(data, pos)
            need_bytes = (expected - got_walk) * 4
            std12 = read_i32_be(data, pos) == 12
            print(
                f"  cont[{chunk_idx}] at rel={pos - sec}: bare_bc={bare_bc} "
                f"need_bytes={need_bytes} std12={std12}",
                flush=True,
            )
            if bare_bc == _CONN_CHUNK_BYTES and pos + 4 + _CONN_CHUNK_BYTES <= sec_end:
                got_walk += _CONN_CHUNK_BYTES // 4
                pos += 4 + _CONN_CHUNK_BYTES
            elif bare_bc == _CONN_CHUNK_BYTES and pos + 8 <= sec_end:
                inner = read_i32_be(data, pos + 4)
                if inner == need_bytes and pos + 8 + need_bytes <= sec_end:
                    got_walk += need_bytes // 4
                    break
                elif pos + 4 + need_bytes <= sec_end:
                    got_walk += need_bytes // 4
                    break
            elif bare_bc == need_bytes and pos + 4 + need_bytes <= sec_end:
                got_walk += need_bytes // 4
                break
            elif bare_bc >= need_bytes and bare_bc % 4 == 0 and pos + 4 + bare_bc <= sec_end:
                got_walk += bare_bc // 4
                pos += 4 + bare_bc
            elif std12 and pos + 8 <= sec_end:
                bc2 = read_i32_be(data, pos + 4)
                if (
                    bc2 > 0
                    and bc2 % 4 == 0
                    and pos + 8 + bc2 + 4 <= sec_end
                    and read_i32_be(data, pos + 8 + bc2) == bc2
                ):
                    got_walk += bc2 // 4
                    pos += 8 + bc2 + 4
                else:
                    print("    -> STOP (std12 invalid)", flush=True)
                    break
            else:
                print("    -> STOP (unrecognized)", flush=True)
                break
            chunk_idx += 1
        print(f"  walk got={got_walk} expected={expected}", flush=True)

        # Sample conn values at chunk boundary if split
        if got >= expected and conn_bc < expected * 4:
            boundary = conn_bc // 4
            conn_parts = [
                np.frombuffer(data, dtype=">u4", count=conn_bc // 4, offset=conn_p).astype(np.int64).copy()
            ]
            _read_conn_continuations(data, conn_p + conn_bc + 4, sec_end, boundary, expected, conn_parts)
            conn = np.concatenate(conn_parts)[:expected]
            n_verts_sec = find_section(data, "LS_Nodes")
            # quick vertex count from summary isn't available; use max index check
            b0 = boundary - 3
            b1 = boundary + 3
            print(f"  conn sample around entry {boundary}: {conn[b0:b1].tolist()}", flush=True)
            print(f"  conn max index: {int(conn.max())}", flush=True)
    print(flush=True)


def main() -> None:
    root = Path(__file__).resolve().parent
    paths = [Path(p) for p in sys.argv[1:]] if len(sys.argv) > 1 else [
        root / "box.gph",
        root / "laptop_simplified_voxel_v9.gph",
        root / "laptop_simplified_voxel_v4.gph",
    ]
    for p in paths:
        if not p.is_file():
            print(f"MISSING: {p}", flush=True)
            continue
        analyze(p)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Validate conn vertex indices against LS_Nodes count."""
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from gph_model import (
    _read_conn_continuations,
    find_section,
    iter_data_blocks,
    open_gph_buffer,
    parse_ls_nodes_vertices,
    read_i32_be,
    section_end,
)


def load_conn(data, sec, sec_end):
    blocks = [(p, bc) for p, bc in iter_data_blocks(data, sec, sec_end) if bc > 0]
    block_sizes = [bc for _, bc in blocks]
    n_faces_block_size = None
    for size, count in Counter(block_sizes).most_common():
        if count >= 3 and size % 4 == 0 and size >= 4:
            n_faces_block_size = size
            break
    n_faces = n_faces_block_size // 4
    triples = [b for b in blocks if b[1] == n_faces_block_size][:3]
    npe_p = triples[2][0]
    npe = np.frombuffer(data, dtype=">u4", count=n_faces, offset=npe_p).astype(np.int64)
    expected = int(npe.sum())

    conn_block = None
    for p, bc in blocks:
        if (p, bc) in triples:
            continue
        if bc % 4 != 0 or bc < 12:
            continue
        if conn_block is None or bc > conn_block[1]:
            conn_block = (p, bc)
    conn_p, conn_bc = conn_block
    parts = [np.frombuffer(data, dtype=">u4", count=conn_bc // 4, offset=conn_p).astype(np.int64).copy()]
    got = len(parts[0])
    if got < expected:
        got, _, _ = _read_conn_continuations(
            data, conn_p + conn_bc + 4, sec_end, got, expected, parts,
        )
    conn = np.concatenate(parts)[:expected]
    return n_faces, npe, conn, blocks, triples


def validate(path: Path) -> None:
    print(f"=== {path.name} ===", flush=True)
    with open_gph_buffer(str(path)) as data:
        _, _, n_verts = parse_ls_nodes_vertices(data, max_preview=0)
        sec = find_section(data, "LS_Links")
        sec_end = section_end(data, sec)
        n_faces, npe, conn, blocks, triples = load_conn(data, sec, sec_end)

        print(f"  blocks ({len(blocks)}):", flush=True)
        for i, (p, bc) in enumerate(blocks):
            role = "?"
            if (p, bc) in triples:
                role = f"triple[{triples.index((p, bc))}]"
            elif bc >= 12 and bc % 4 == 0:
                role = "conn?" if bc >= 1073741824 else "other"
            print(f"    [{i}] bc={bc:12d} role={role}", flush=True)

        print(f"  n_faces={n_faces} n_vertices={n_verts} conn_len={len(conn)}", flush=True)
        print(f"  conn min={int(conn.min())} max={int(conn.max())}", flush=True)

        bad = conn >= n_verts
        n_bad = int(bad.sum())
        print(f"  indices >= n_vertices: {n_bad} ({100.0 * n_bad / len(conn):.4f}%)", flush=True)
        if n_bad:
            bad_vals = conn[bad]
            print(f"  bad value samples: {np.unique(bad_vals)[:20].tolist()}", flush=True)

        # Check chunk boundary alignment (1 GiB = 268435456 entries)
        chunk = 268435456
        if len(conn) > chunk:
            for boundary in [chunk, 2 * chunk]:
                if boundary >= len(conn):
                    continue
                b = boundary
                window = conn[b - 5 : b + 5]
                print(f"  at entry {b}: {window.tolist()}", flush=True)
                bad_at = int((conn[b : b + 1000] >= n_verts).sum())
                print(f"  bad in next 1000 after boundary: {bad_at}", flush=True)

        # face sanity: each face's nodes should be valid
        offsets = np.empty(n_faces + 1, dtype=np.int64)
        offsets[0] = 0
        np.cumsum(npe, out=offsets[1:])
        bad_faces = 0
        sample_bad = []
        for fi in range(min(n_faces, 500000)):
            sl = conn[offsets[fi] : offsets[fi + 1]]
            if (sl >= n_verts).any():
                bad_faces += 1
                if len(sample_bad) < 5:
                    sample_bad.append((fi, sl.tolist()))
        if n_faces > 500000:
            # extrapolate from sample
            rate = bad_faces / 500000
            print(f"  bad faces (first 500k): {bad_faces} (~{rate*100:.2f}%)", flush=True)
        else:
            print(f"  bad faces: {bad_faces}", flush=True)
        for fi, sl in sample_bad:
            print(f"    face {fi}: {sl}", flush=True)
    print(flush=True)


def main():
    root = Path(__file__).resolve().parent
    paths = [Path(p) for p in sys.argv[1:]] if len(sys.argv) > 1 else [
        root / "box.gph",
        root / "laptop_simplified_voxel_v9.gph",
        root / "laptop_simplified_voxel_v4.gph",
    ]
    for p in paths:
        validate(p)


if __name__ == "__main__":
    main()

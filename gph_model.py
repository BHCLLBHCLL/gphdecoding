#!/usr/bin/env python3
"""
GPH data model - parse binary to editable tree, support partial save.
"""

import struct
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np

_LARGE_GPH_BYTES = 512 * 1024 * 1024  # mmap threshold (see gph2cgns.parse_gph_mesh)
_CONN_CHUNK_BYTES = 1073741824  # 1 GiB cap per LS_Links conn payload block


def read_i32_be(data: bytes, pos: int) -> int:
    return int.from_bytes(data[pos : pos + 4], "big")


def read_f32_be(data: bytes, pos: int) -> float:
    return struct.unpack(">f", data[pos : pos + 4])[0]


def read_f64_be(data: bytes, pos: int) -> float:
    """Read a standard big-endian IEEE-754 float64."""
    return struct.unpack(">d", data[pos : pos + 8])[0]


def read_f64_wr(data: bytes, pos: int) -> float:
    """Read float64 stored as word-reversed: [lower_32bit_BE][upper_32bit_BE]."""
    lower = int.from_bytes(data[pos : pos + 4], "big")
    upper = int.from_bytes(data[pos + 4 : pos + 8], "big")
    combined = ((upper << 32) | lower).to_bytes(8, "big")
    return struct.unpack(">d", combined)[0]


def _looks_like_coords(values: list) -> bool:
    """Heuristic: do the magnitudes look like physical CFD coordinates?

    Used to auto-detect the GPH dialect (standard big-endian float64 vs the
    legacy word-reversed encoding).
    """
    if not values:
        return False
    import math
    for v in values:
        if not math.isfinite(v):
            return False
        a = abs(v)
        if a != 0.0 and (a > 1e6 or a < 1e-30):
            return False
    return True


@contextmanager
def open_gph_buffer(filepath: str):
    """Yield a bytes-like buffer for a GPH file.

    Files larger than 512 MiB are memory-mapped (as in ``gph2cgns.py``) so
    multi-gigabyte meshes such as ``laptop_simplified_voxel_v4.gph``,
    ``laptop_simplified_denser_v2_gph.gph``, or ``laptop_simplified_voxel_v6.gph``
    can be inspected without a full-RAM copy.
    """
    size = Path(filepath).stat().st_size
    if size <= _LARGE_GPH_BYTES:
        with open(filepath, "rb") as f:
            yield f.read()
        return
    import mmap
    f = open(filepath, "rb")
    try:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            yield mm
        finally:
            mm.close()
    finally:
        f.close()


def _f64_be_array(buf, offset: int, count: int) -> np.ndarray:
    return np.frombuffer(buf, dtype=">f8", count=count, offset=offset).copy()


def _f64_wr_array(buf, offset: int, count: int) -> np.ndarray:
    raw = np.frombuffer(buf, dtype=">u4", count=count * 2, offset=offset)
    lower = raw[0::2].astype(np.uint64)
    upper = raw[1::2].astype(np.uint64)
    bits = (upper << 32) | lower
    return bits.view(">f8").astype(np.float64)


def _f32_be_array(buf, offset: int, count: int) -> np.ndarray:
    return np.frombuffer(buf, dtype=">f4", count=count, offset=offset).astype(np.float64).copy()


def _read_conn_continuations(data, pos: int, sec_end: int, got: int,
                             expected: int,
                             conn_parts: Optional[list] = None) -> tuple[int, int, int]:
    """Read conn split continuations after the primary conn block.

    Very large meshes cap each conn payload at 1 GiB.  Continuation blocks use
    bare ``[I4=byte_count][payload]`` (no ``[I4=12]`` header).  Multiple full
    1 GiB chunks may appear; the final short chunk repeats the 1 GiB marker
    followed by the actual payload byte count and then the payload:
    ``[I4=1GiB][I4=need_bytes][payload]`` (e.g. ``tests/box.gph``,
    ``laptop_simplified_denser_v2_gph.gph``).

    Returns ``(new_got, final_pos, n_continuation_chunks)``.
    """
    n_continuations = 0
    while got < expected and pos + 4 <= sec_end:
        need_bytes = (expected - got) * 4
        bare_bc = read_i32_be(data, pos)

        if (bare_bc == _CONN_CHUNK_BYTES
                and pos + 4 + _CONN_CHUNK_BYTES <= sec_end):
            n = _CONN_CHUNK_BYTES // 4
            if conn_parts is not None:
                conn_parts.append(
                    np.frombuffer(data, dtype=">u4", count=n, offset=pos + 4)
                    .astype(np.int64).copy())
            got += n
            pos += 4 + _CONN_CHUNK_BYTES
            n_continuations += 1
            continue

        if (bare_bc == _CONN_CHUNK_BYTES
                and need_bytes < _CONN_CHUNK_BYTES
                and pos + 8 <= sec_end):
            inner_bc = read_i32_be(data, pos + 4)
            if (inner_bc == need_bytes
                    and pos + 8 + need_bytes <= sec_end):
                n = need_bytes // 4
                if conn_parts is not None:
                    conn_parts.append(
                        np.frombuffer(data, dtype=">u4", count=n, offset=pos + 8)
                        .astype(np.int64).copy())
                got += n
                pos += 8 + need_bytes
                n_continuations += 1
                break

        if (bare_bc == _CONN_CHUNK_BYTES
                and pos + 4 + need_bytes <= sec_end):
            n = need_bytes // 4
            if conn_parts is not None:
                conn_parts.append(
                    np.frombuffer(data, dtype=">u4", count=n, offset=pos + 4)
                    .astype(np.int64).copy())
            got += n
            pos += 4 + need_bytes
            n_continuations += 1
            break

        if (bare_bc == need_bytes
                and pos + 4 + need_bytes <= sec_end):
            n = need_bytes // 4
            if conn_parts is not None:
                conn_parts.append(
                    np.frombuffer(data, dtype=">u4", count=n, offset=pos + 4)
                    .astype(np.int64).copy())
            got += n
            n_continuations += 1
            break

        if (bare_bc >= need_bytes and bare_bc % 4 == 0
                and pos + 4 + bare_bc <= sec_end):
            n = bare_bc // 4
            if conn_parts is not None:
                conn_parts.append(
                    np.frombuffer(data, dtype=">u4", count=n, offset=pos + 4)
                    .astype(np.int64).copy())
            got += n
            pos += 4 + bare_bc
            n_continuations += 1
            continue

        if read_i32_be(data, pos) == 12 and pos + 8 <= sec_end:
            bc2 = read_i32_be(data, pos + 4)
            if (bc2 > 0 and bc2 % 4 == 0
                    and pos + 8 + bc2 + 4 <= sec_end
                    and read_i32_be(data, pos + 8 + bc2) == bc2):
                n = bc2 // 4
                if conn_parts is not None:
                    conn_parts.append(
                        np.frombuffer(data, dtype=">u4", count=n, offset=pos + 8)
                        .astype(np.int64).copy())
                got += n
                pos += 8 + bc2 + 4
                n_continuations += 1
                continue
        break
    return got, pos, n_continuations


def _conn_payload_size(data, blocks, triples, conn_total_expected: int,
                       sec_end: int) -> tuple[int, bool, int]:
    """Return ``(conn_entry_count, split, n_chunks)`` for the LS_Links conn array."""
    conn_block = None
    for p, bc in blocks:
        if (p, bc) in triples:
            continue
        if bc % 4 != 0:
            continue
        if bc // 4 == conn_total_expected:
            return conn_total_expected, False, 1
        if bc < 12:
            continue
        if conn_block is None or bc > conn_block[1]:
            conn_block = (p, bc)
    if conn_block is None:
        return 0, False, 0
    conn_p, conn_bc = conn_block
    got = conn_bc // 4
    if got >= conn_total_expected:
        return conn_total_expected, False, 1
    pos = conn_p + conn_bc + 4
    got, _, n_extra = _read_conn_continuations(
        data, pos, sec_end, got, conn_total_expected,
    )
    n_chunks = 1 + n_extra if n_extra else 1
    split = got >= conn_total_expected and got > conn_bc // 4
    return got, split, n_chunks


# ── Shared section scanners (aligned with gph2cgns.py) ───────────────────────

_SECTION_BOUNDARY_NAMES = [
    "FileRevision", "Application", "ApplicationVersion", "ReleaseDate",
    "GridType", "Dimension", "Bias", "Date", "Comments", "Cycle",
    "Unused", "Encoding", "HeaderDataEnd", "OverlapStart_0",
    "LS_CvolIdOfElements", "LS_Links", "LS_Nodes", "LS_SurfaceRegions",
    "LS_SolverUnusedRegions", "LS_VolumeRegions", "LS_Parts",
    "LS_Assemblies", "Element_InformationFlag", "OverlapEnd",
]


def find_section(data: bytes, name: str) -> int:
    """Return offset of the I4=32 marker before *name*, or -1."""
    name_padded = name.ljust(32).encode("ascii")
    idx = data.find(name_padded)
    if idx < 4:
        return -1
    if read_i32_be(data, idx - 4) == 32:
        return idx - 4
    return -1


def section_end(data: bytes, sec_start: int) -> int:
    best = len(data)
    for name in _SECTION_BOUNDARY_NAMES:
        off = find_section(data, name)
        if off > sec_start and off < best:
            best = off
    return best


def iter_data_blocks(data: bytes, sec_start: int, sec_end: int):
    """Yield ``(payload_start, byte_count)`` for each data block in a section."""
    pos = sec_start + 40
    n = len(data)
    while pos + 8 <= sec_end and pos + 8 <= n:
        if read_i32_be(data, pos) != 12:
            pos += 4
            continue
        v = read_i32_be(data, pos + 4)
        if v in (4, 8) and pos + 16 <= sec_end:
            dim0 = read_i32_be(data, pos + 8)
            dim1 = read_i32_be(data, pos + 12)
            if 0 < dim0 < 10_000_000 and 0 < dim1 < 10_000_000:
                pos += 16
                continue
        bc = v
        if bc <= 0 or pos + 8 + bc + 4 > sec_end:
            pos += 4
            continue
        payload_end = pos + 8 + bc
        if read_i32_be(data, payload_end) != bc:
            pos += 4
            continue
        yield pos + 8, bc
        pos = payload_end + 4


def parse_ls_cvol_ids(data: bytes) -> Optional["np.ndarray"]:
    """Parse LS_CvolIdOfElements -> I4[n_cells] (largest I4 block in section)."""
    sec_start = find_section(data, "LS_CvolIdOfElements")
    if sec_start < 0:
        return None
    sec_end = section_end(data, sec_start)
    best: Optional[tuple[int, int]] = None
    for p, bc in iter_data_blocks(data, sec_start, sec_end):
        if bc % 4 == 0 and bc >= 4:
            if best is None or bc > best[1]:
                best = (p, bc)
    if best is None:
        return None
    p, bc = best
    return np.frombuffer(data, dtype=">i4", count=bc // 4, offset=p).astype(np.int64).copy()


def parse_ls_nodes_vertices(
    data: bytes,
    max_preview: int = 3,
) -> tuple[Optional[list[tuple[float, float, float]]], str, int]:
    """Parse LS_Nodes -> (coord_sample, dialect_label, n_vertices).

    For large meshes only *max_preview* coordinates are materialised (for
    display); the full vertex count is always returned.  Delegates to
    :func:`parse_ls_nodes_xyz` (float32 / float64 / word-reversed).
    """
    xyz, n_vertices = parse_ls_nodes_xyz(data)
    if xyz is None or n_vertices == 0:
        return None, "", 0

    sec_start = find_section(data, "LS_Nodes")
    sec_end = section_end(data, sec_start)
    elem_hint = ls_nodes_descriptor_elem_bytes(data, sec_start, sec_end)
    if elem_hint == 4:
        dialect = "big-endian float32"
    else:
        layout = _ls_nodes_coordinate_layout(data)
        dialect = layout["dialect"] if layout else "standard BE float64"

    n_show = min(n_vertices, max_preview)
    sample = [
        (float(xyz[i, 0]), float(xyz[i, 1]), float(xyz[i, 2]))
        for i in range(n_show)
    ]
    return sample, dialect, n_vertices


def parse_ls_links_summary(data: bytes) -> Optional[dict]:
    """Return a short topology summary dict for LS_Links."""
    sec_start = find_section(data, "LS_Links")
    if sec_start < 0:
        return None
    sec_end = section_end(data, sec_start)
    blocks = [(p, bc) for p, bc in iter_data_blocks(data, sec_start, sec_end) if bc > 0]
    if not blocks:
        return None
    from collections import Counter
    block_sizes = [bc for _, bc in blocks]
    n_faces_block_size = None
    for size, count in Counter(block_sizes).most_common():
        if count >= 3 and size % 4 == 0 and size >= 4:
            n_faces_block_size = size
            break
    if n_faces_block_size is None:
        return None
    n_faces = n_faces_block_size // 4
    triples = [b for b in blocks if b[1] == n_faces_block_size][:3]
    if len(triples) < 3:
        return None
    owner_p, _ = triples[0]
    neigh_p, _ = triples[1]
    npe_p, _ = triples[2]
    npe = np.frombuffer(data, dtype=">u4", count=n_faces, offset=npe_p).astype(np.int64)
    conn_total = int(npe.sum())
    neigh_raw = np.frombuffer(data, dtype=">u4", count=n_faces, offset=neigh_p)
    boundary = int((neigh_raw == 0xFFFFFFFF).sum())
    owner = np.frombuffer(data, dtype=">u4", count=n_faces, offset=owner_p)
    n_cells = int(owner.max()) + 1
    conn_got, conn_split, conn_chunks = _conn_payload_size(
        data, blocks, triples, conn_total, sec_end,
    )
    return {
        "n_faces": n_faces,
        "n_cells": n_cells,
        "boundary_faces": boundary,
        "npe_min": int(npe.min()),
        "npe_max": int(npe.max()),
        "conn_entries": conn_total,
        "conn_got": conn_got,
        "conn_split": conn_split,
        "conn_chunks": conn_chunks,
        "conn_complete": conn_got >= conn_total,
        "polyhedral": int(npe.max()) > 3,
    }


def _ls_parts_name_blocks(
    data: bytes, sec_start: int, sec_end: int,
) -> list[tuple[str, int, int]]:
    """Return ``[(name, name_header_pos, after_trailer_pos), ...]`` in file order.

    Name blocks are ASCII (NUL/space-padded) data blocks containing at least
    one alphabetic character.  GPH files pad names to 255 bytes; FPH files
    use a smaller padding (e.g. 23 bytes).  The block size is therefore not
    hard-coded — the printable-ASCII + alpha heuristic is sufficient to
    distinguish name blocks from the surrounding cvol-descriptor metadata.
    """
    name_blocks: list[tuple[str, int, int]] = []
    for p, bc in iter_data_blocks(data, sec_start, sec_end):
        if bc <= 0 or bc > 512:
            continue
        raw = data[p : p + bc]
        if not all(b == 0 or 32 <= b < 127 for b in raw):
            continue
        name = raw.decode("ascii", errors="replace").strip("\x00").rstrip()
        if not name or not any(c.isalpha() for c in name):
            continue
        name_blocks.append((name, p - 8, p + bc + 4))
    return name_blocks


def _scan_cvol_descriptor_chain(data: bytes, start: int, end: int) -> list[int]:
    """Collect every ``[12, 4, X, 4]`` value in ``[start, end)`` in file order."""
    chain: list[int] = []
    pos = start
    while pos + 16 <= end:
        if (read_i32_be(data, pos) == 12
                and read_i32_be(data, pos + 4) == 4
                and read_i32_be(data, pos + 12) == 4):
            chain.append(read_i32_be(data, pos + 8))
        pos += 4
    return chain


# Part → either one cvol_id or a membership set (composite / background parts).
PartCvolSpec = int | frozenset[int]


def format_part_cvol_spec(spec: PartCvolSpec) -> str:
    if isinstance(spec, frozenset):
        ids = sorted(spec)
        if len(ids) <= 10:
            return "{" + ", ".join(str(i) for i in ids) + "}"
        return f"{{{ids[0]}..{ids[-1]} ... n={len(ids)}}}"
    return str(spec)


def part_cvol_cell_mask(cvol_id: "np.ndarray", spec: PartCvolSpec) -> "np.ndarray":
    """Boolean mask of cells belonging to a Part (single id or id set)."""
    if isinstance(spec, frozenset):
        if not spec:
            return np.zeros(len(cvol_id), dtype=bool)
        return np.isin(cvol_id, list(spec))
    return cvol_id == spec


def _parse_part_cvol_membership(
    data: bytes,
    start: int,
    end: int,
    actual_set: Optional[set[int]],
) -> Optional[frozenset[int]]:
    """Parse composite Part layout: ``[12,4,N,4]`` + ``I4[N]`` cvol_id list.

    Used by background parts such as ``air_domain`` in multi-region laptop
    models: ``N`` is the list length (not a cvol_id), followed by the full
    set of cvol_ids owned by that Part.
    """
    chain = _scan_cvol_descriptor_chain(data, start, end)
    if not chain:
        return None
    chain_counts = set(chain)
    for p, bc in iter_data_blocks(data, start, end):
        if bc < 8 or bc % 4 != 0:
            continue
        n = bc // 4
        if n not in chain_counts:
            continue
        vals = [int(x) for x in np.frombuffer(data, dtype=">i4", count=n, offset=p)]
        if len(vals) != n or len(set(vals)) != n:
            continue
        if actual_set is not None and not all(v in actual_set for v in vals):
            continue
        if n >= 2:
            return frozenset(vals)
    return None


def _resolve_single_part_cvol(
    chain: list[int],
    actual_set: Optional[set[int]],
) -> int:
    """Map a simple ``[1, cvol_id]`` descriptor chain to one cvol_id.

    The chain is ``[1, <cvol_id>]`` — a leading ``1`` marker followed by
    the part's cvol_id.  The last value is therefore always the cvol_id.
    It is returned even when it is absent from *actual_set* (the part
    simply has 0 cells in this mesh); falling back to an earlier chain
    value would incorrectly alias the part to another region's cells.
    """
    if not chain:
        return 1
    return int(chain[-1])


def _resolve_part_cvol_ids(
    parts: list[tuple[str, list[int]]],
    actual_set: Optional[set[int]],
) -> list[tuple[str, int]]:
    """Map each Part name to a cvol_id using descriptor chains + global validation.

    Each Part record stores a post-name descriptor chain that typically looks
    like ``[1, <cvol_id>]`` (leading ``1`` markers plus the opaque id).  The
    primary rule is therefore: pick the **last** chain value that belongs to
    the mesh's actual cvol_id set (from ``LS_CvolIdOfElements``), not merely
    the last ``[12,4,X,4]`` regardless of *X*.

    When that primary mapping is not unique / does not cover the actual set,
    fall back in order to: (a) exactly-one-candidate-per-part, (b) sequential
    ``1..N`` only when the actual set is ``{1, …, N}``, (c) raw last chain
    element, (d) best-effort primary picks.
    """
    n = len(parts)
    if n == 0:
        return []

    def _pick_last_in_set(chain: list[int]) -> Optional[int]:
        if not chain:
            return None
        if actual_set:
            for value in reversed(chain):
                iv = int(value)
                if iv in actual_set:
                    return iv
        return int(chain[-1])

    def _mapping_ok(mapping: list[tuple[str, Optional[int]]]) -> bool:
        ids = [cv for _, cv in mapping if cv is not None]
        if len(ids) != n or len(set(ids)) != n:
            return False
        if actual_set is None:
            return True
        if not set(ids) <= actual_set:
            return False
        return len(actual_set) != n or set(ids) == actual_set

    primary = [(name, _pick_last_in_set(chain)) for name, chain in parts]
    if _mapping_ok(primary):
        return [(name, int(cv)) for name, cv in primary]

    single: list[tuple[str, int]] = []
    for name, chain in parts:
        cands = [int(v) for v in chain
                 if actual_set is None or int(v) in actual_set]
        if len(cands) != 1:
            break
        single.append((name, cands[0]))
    else:
        if _mapping_ok(single):
            return single

    if actual_set is not None and actual_set == set(range(1, n + 1)):
        return [(name, idx) for idx, (name, _) in enumerate(parts, start=1)]

    raw_last = [(name, int(chain[-1])) for name, chain in parts if chain]
    if len(raw_last) == n and _mapping_ok(raw_last):
        return raw_last

    out: list[tuple[str, int]] = []
    for idx, (name, cv) in enumerate(primary, start=1):
        out.append((name, int(cv) if cv is not None else idx))
    return out


def parse_ls_parts(
    data: bytes,
    cvol_id: Optional["np.ndarray"] = None,
) -> list[tuple[str, PartCvolSpec]]:
    """Parse LS_Parts → ``[(part_name, cvol_spec), ...]`` in file order.

    ``cvol_spec`` is either a single **cvol_id** (``int``) or a ``frozenset``
    of ids for composite Parts.  Simple parts use a post-name chain
    ``[1, cvol_id]``; background parts (e.g. ``air_domain`` in multi-region
    laptop meshes) store ``[12,4,N,4]`` where *N* is the list length, then
    an ``I4[N]`` block listing every cvol_id that belongs to the Part.
    """
    sec_start = find_section(data, "LS_Parts")
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    name_blocks = _ls_parts_name_blocks(data, sec_start, sec_end)

    actual_set: Optional[set[int]] = None
    if cvol_id is not None and len(cvol_id) > 0:
        actual_set = {int(x) for x in np.unique(cvol_id)}

    out: list[tuple[str, PartCvolSpec]] = []
    for i, (name, _, after_trailer) in enumerate(name_blocks):
        scan_end = name_blocks[i + 1][1] if i + 1 < len(name_blocks) else sec_end
        membership = _parse_part_cvol_membership(
            data, after_trailer, scan_end, actual_set,
        )
        if membership is not None:
            out.append((name, membership))
            continue
        chain = _scan_cvol_descriptor_chain(data, after_trailer, scan_end)
        out.append((name, _resolve_single_part_cvol(chain, actual_set)))
    return out


def parse_ls_string_list(data: bytes, section_name: str) -> list[str]:
    sec_start = find_section(data, section_name)
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    out: list[str] = []
    for p, bc in iter_data_blocks(data, sec_start, sec_end):
        raw = data[p : p + bc]
        if all(b == 0 or 32 <= b < 127 for b in raw):
            s = raw.decode("ascii", errors="replace").strip("\x00").rstrip()
            if s:
                out.append(s)
    return out


def parse_ls_surface_regions(data: bytes) -> list[tuple[str, "np.ndarray"]]:
    """Parse LS_SurfaceRegions -> [(name, face_ids), ...]."""
    sec_start = find_section(data, "LS_SurfaceRegions")
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    blocks = list(iter_data_blocks(data, sec_start, sec_end))
    out: list[tuple[str, "np.ndarray"]] = []
    i = 0
    while i + 2 < len(blocks):
        p_n, bc_n = blocks[i]
        p_i, bc_i = blocks[i + 1]
        p_w, bc_w = blocks[i + 2]
        name_raw = data[p_n : p_n + bc_n]
        if not all(b == 0 or 32 <= b < 127 for b in name_raw):
            i += 1
            continue
        name = name_raw.decode("ascii", errors="replace").strip("\x00").rstrip()
        if not name:
            i += 1
            continue
        if bc_i > 0 and bc_i == bc_w and bc_i % 4 == 0:
            face_ids = np.frombuffer(
                data, dtype=">i4", count=bc_i // 4, offset=p_i,
            ).astype(np.int64).copy()
            out.append((name, face_ids))
            i += 3
        else:
            i += 1
    return out


def parse_ls_surface_regions_summary(data: bytes) -> list[tuple[str, int]]:
    return [(name, int(face_ids.size)) for name, face_ids in parse_ls_surface_regions(data)]


def parse_ls_assemblies_summary(data: bytes) -> dict:
    sec_start = find_section(data, "LS_Assemblies")
    empty = {"has_assemblies": False, "root_empty_prefix": None, "part_paths": {}}
    if sec_start < 0:
        return empty
    sec_end = section_end(data, sec_start)
    xml_bytes = b""
    for p, bc in iter_data_blocks(data, sec_start, sec_end):
        chunk = data[p : p + bc]
        if chunk.lstrip().startswith(b"<?xml") or b"<part" in chunk:
            xml_bytes = chunk
            break
    if not xml_bytes:
        return empty
    try:
        import xml.etree.ElementTree as ET
        root = ET.fromstring(xml_bytes.decode("utf-8", errors="replace"))
    except Exception:
        return empty
    part_paths: dict[str, Optional[str]] = {}
    has_assemblies = any(True for _ in root.iter("assembly"))

    def _walk(node, ancestors: list[str]):
        for child in node:
            if child.tag == "assembly":
                aname = child.get("name", "")
                _walk(child, ancestors + [aname] if aname else ancestors)
            elif child.tag == "part":
                pname = child.get("name", "")
                if pname:
                    part_paths[pname] = ".".join(ancestors + [pname]) if ancestors else None

    _walk(root, [])
    root_parts_count = sum(1 for part in root.findall("part") if part.get("name"))
    root_empty_prefix: Optional[str] = None
    if root_parts_count > 0:
        top_asm = next(iter(root.findall("assembly")), None)
        if top_asm is not None:
            empty_asm_names: list[str] = []
            for child in top_asm.findall("assembly"):
                if (len(child.findall("assembly")) == 0
                        and len(child.findall("part")) == 0):
                    name = child.get("name", "")
                    if name:
                        empty_asm_names.append(name)
            if len(empty_asm_names) >= root_parts_count:
                root_empty_prefix = ".".join(empty_asm_names[:root_parts_count])
    return {
        "has_assemblies": has_assemblies,
        "root_empty_prefix": root_empty_prefix,
        "part_paths": part_paths,
    }


def classify_volume_region_cells(
    zone_name: str,
    parts_with_cvol: list[tuple[str, PartCvolSpec]],
    cvol_id: Optional["np.ndarray"],
    n_cells: int,
) -> "np.ndarray":
    """Return a cell mask for a volume/part region name."""
    all_mask = np.ones(n_cells, dtype=bool)
    if zone_name == "FluidRegion":
        return all_mask
    if cvol_id is None or len(cvol_id) != n_cells or not parts_with_cvol:
        return all_mask

    name_to_cvol = {name: cv for name, cv in parts_with_cvol}
    if zone_name.startswith("@VPartRegion_"):
        rem = zone_name[len("@VPartRegion_") :].split("[", 1)[0]
        if rem in name_to_cvol:
            return part_cvol_cell_mask(cvol_id, name_to_cvol[rem])
    if zone_name.startswith("FPHPARTS."):
        candidate = zone_name[len("FPHPARTS.") :].rsplit(".", 1)[-1]
        if candidate in name_to_cvol:
            return part_cvol_cell_mask(cvol_id, name_to_cvol[candidate])
    matches = sorted(
        (p for p, _ in parts_with_cvol if p and p in zone_name),
        key=len,
        reverse=True,
    )
    if matches:
        return part_cvol_cell_mask(cvol_id, name_to_cvol[matches[0]])
    return all_mask


# Minimum plausible max |coordinate| for meter-scale CFD meshes.  float32
# payload misread as float64 produces denormal magnitudes ~1e-13 (see
# ``tests/tr03_9.fph``); real vertex coordinates are typically >= ~1e-4 m.
_COORD_MIN_ABSMAX = 1e-4
_COORD_MAX_ABSMAX = 1e6


def _score_coord_axes(axes: list["np.ndarray"]) -> float:
    """Lower score = more plausible CFD vertex coordinate axes.

    Penalises non-finite values, coordinate magnitudes outside
    ~[``_COORD_MIN_ABSMAX``, ``_COORD_MAX_ABSMAX``], a high fraction of such
    outliers (typical of wrong float32/float64 decode), and grossly mismatched
    axis scales.
    """
    score = 0.0
    axis_absmax: list[float] = []
    for ax in axes:
        arr = np.asarray(ax, dtype=np.float64)
        finite = np.isfinite(arr)
        if not finite.all():
            score += 1e30
            continue
        absv = np.abs(arr[finite])
        if absv.size == 0:
            axis_absmax.append(0.0)
            continue
        absmax = float(np.max(absv))
        axis_absmax.append(absmax)
        if absmax > _COORD_MAX_ABSMAX or (
                absmax < _COORD_MIN_ABSMAX and absmax != 0.0
        ):
            score += absmax + (1.0 / max(absmax, 1e-300))
        else:
            score += absmax
        # Per-value tiny-outlier fraction only when the axis absmax itself
        # looks misdecoded (e.g. float32 payload read as float64 → ~1e-13).
        # Do not penalise legitimate meshes whose absmax is O(1) but many
        # vertices lie near the origin (|x| < 1e-4).
        if absmax < _COORD_MIN_ABSMAX and absmax != 0.0:
            bad_frac = float(
                ((absv > _COORD_MAX_ABSMAX)
                 | ((absv < _COORD_MIN_ABSMAX) & (absv != 0.0))).mean()
            )
            if bad_frac > 0.01:
                score += 1e20 * bad_frac
        elif absmax > _COORD_MAX_ABSMAX:
            bad_frac = float((absv > _COORD_MAX_ABSMAX).mean())
            if bad_frac > 0.01:
                score += 1e20 * bad_frac
    pos = [v for v in axis_absmax if v > 0]
    if len(pos) >= 2:
        ratio = max(pos) / min(pos)
        if ratio > 1e6:
            score += ratio
    return score


def ls_nodes_descriptor_elem_bytes(
    data, sec_start: int, sec_end: int,
) -> Optional[int]:
    """Return element size (4=float32, 8=float64) from LS_Nodes type descriptors.

    Only ``[12, type, n_verts, dim1]`` with ``n_verts > 1`` are counted so
    metadata markers like ``[12, 4, 1, 1]`` do not skew the vote.
    """
    counts = {4: 0, 8: 0}
    pos = sec_start + 40
    n = len(data)
    while pos + 16 <= sec_end and pos + 16 <= n:
        if read_i32_be(data, pos) == 12:
            tc = read_i32_be(data, pos + 4)
            if tc in (4, 8):
                dim0 = read_i32_be(data, pos + 8)
                dim1 = read_i32_be(data, pos + 12)
                if dim0 > 1 and 0 < dim1 < 10_000_000:
                    counts[tc] += 1
        pos += 4
    if counts[8] > counts[4]:
        return 8
    if counts[4] > counts[8]:
        return 4
    return None


def ls_nodes_vertex_count_from_descriptors(
    data, sec_start: int, sec_end: int,
) -> Optional[int]:
    """Return vertex count from ``[12, type, n_verts, …]`` in LS_Nodes."""
    best = 0
    pos = sec_start + 40
    n = len(data)
    while pos + 16 <= sec_end and pos + 16 <= n:
        if read_i32_be(data, pos) == 12:
            tc = read_i32_be(data, pos + 4)
            if tc in (4, 8):
                dim0 = read_i32_be(data, pos + 8)
                dim1 = read_i32_be(data, pos + 12)
                if dim0 > 1 and 0 < dim1 < 10_000_000:
                    best = max(best, dim0)
        pos += 4
    return best if best > 0 else None


_COORD_SCORE_SAMPLE = 256
_ELEM_PRIOR_MISMATCH = 1e15
_F32_ON_F64_ALIGNED_PRIOR = 10.0


def parse_ls_nodes_xyz(data: bytes) -> tuple[Optional["np.ndarray"], int]:
    """Parse LS_Nodes → ``(xyz float64 N×3, n_vertices)``.

    Supports standard BE float64, word-reversed float64, and BE float32
    (FPH / ``tests/tr03_9.fph``).
    """
    sec_start = find_section(data, "LS_Nodes")
    if sec_start < 0:
        return None, 0
    sec_end = section_end(data, sec_start)

    blocks = list(iter_data_blocks(data, sec_start, sec_end))
    f_blocks = [(p, bc) for p, bc in blocks if bc >= 4 and bc % 4 == 0]
    if len(f_blocks) < 3:
        return None, 0

    sizes = [bc for _, bc in f_blocks]
    target = max(set(sizes), key=sizes.count)
    trio = [(p, bc) for p, bc in f_blocks if bc == target][:3]
    if len(trio) < 3:
        return None, 0

    bc = trio[0][1]
    elem_hint = ls_nodes_descriptor_elem_bytes(data, sec_start, sec_end)
    n_desc = ls_nodes_vertex_count_from_descriptors(data, sec_start, sec_end)

    def _ranked_score(sample_axes: list["np.ndarray"], elem_bytes: int) -> float:
        s = _score_coord_axes(sample_axes)
        if elem_hint is not None and elem_bytes != elem_hint:
            s += _ELEM_PRIOR_MISMATCH
        elif elem_hint is None and elem_bytes == 4 and bc % 8 == 0:
            s += _F32_ON_F64_ALIGNED_PRIOR
        return s

    ranked: list[tuple[float, str]] = []

    if bc % 8 == 0:
        n_f64 = bc // 8
        if n_desc is None or n_desc == n_f64:
            n_sample = min(n_f64, _COORD_SCORE_SAMPLE)
            ranked.append((
                _ranked_score([_f64_be_array(data, p, n_sample) for p, _ in trio], 8),
                "be",
            ))
            ranked.append((
                _ranked_score([_f64_wr_array(data, p, n_sample) for p, _ in trio], 8),
                "wr",
            ))

    if bc % 4 == 0 and elem_hint != 8:
        n_f32 = n_desc if n_desc is not None else bc // 4
        if n_desc is None or n_desc == bc // 4:
            n_sample = min(n_f32, _COORD_SCORE_SAMPLE)
            ranked.append((
                _ranked_score([_f32_be_array(data, p, n_sample) for p, _ in trio], 4),
                "f32",
            ))

    if not ranked:
        return None, 0

    _, kind = min(ranked, key=lambda item: item[0])
    if kind == "be":
        n_vertices = n_desc if n_desc is not None else bc // 8
        axes = [_f64_be_array(data, p, n_vertices) for p, _ in trio]
        is_wr = False
    elif kind == "wr":
        n_vertices = n_desc if n_desc is not None else bc // 8
        axes = [_f64_wr_array(data, p, n_vertices) for p, _ in trio]
        is_wr = True
    else:
        n_vertices = n_desc if n_desc is not None else bc // 4
        axes = [_f32_be_array(data, p, n_vertices) for p, _ in trio]
        is_wr = False

    xyz = np.column_stack(axes)
    if is_wr:
        xyz = xyz[:, [0, 2, 1]]
    return xyz, n_vertices


def _ls_nodes_coordinate_layout(data) -> Optional[dict]:
    sec_start = find_section(data, "LS_Nodes")
    if sec_start < 0:
        return None
    sec_end = section_end(data, sec_start)
    f64_blocks = [
        (p, bc) for p, bc in iter_data_blocks(data, sec_start, sec_end)
        if bc >= 8 and bc % 8 == 0
    ]
    if len(f64_blocks) < 3:
        return None
    sizes = [bc for _, bc in f64_blocks]
    target = max(set(sizes), key=sizes.count)
    trio = [(p, bc) for p, bc in f64_blocks if bc == target][:3]
    if len(trio) < 3:
        return None
    n_vertices = trio[0][1] // 8
    n_sample = min(n_vertices, 256)
    axes_be = [_f64_be_array(data, p, n_sample) for p, _ in trio]
    axes_wr = [_f64_wr_array(data, p, n_sample) for p, _ in trio]
    if _score_coord_axes(axes_be) <= _score_coord_axes(axes_wr):
        return {
            "blocks": trio,
            "n_vertices": n_vertices,
            "dialect": "standard BE float64",
            "word_reversed": False,
            "perm": (0, 1, 2),
        }
    return {
        "blocks": trio,
        "n_vertices": n_vertices,
        "dialect": "word-reversed float64",
        "word_reversed": True,
        "perm": (0, 2, 1),
    }


def _read_vertices_by_id(data, vertex_ids: "np.ndarray") -> tuple["np.ndarray", dict]:
    layout = _ls_nodes_coordinate_layout(data)
    if layout is None or vertex_ids.size == 0:
        return np.empty((0, 3), dtype=float), {}
    n_vertices = int(layout["n_vertices"])
    ids = np.asarray(vertex_ids, dtype=np.int64)
    ids = np.clip(ids, 0, max(n_vertices - 1, 0))
    coords_file = np.empty((ids.size, 3), dtype=float)
    for axis_idx, (p, _) in enumerate(layout["blocks"]):
        if layout["word_reversed"]:
            vals = [read_f64_wr(data, p + int(i) * 8) for i in ids]
        else:
            vals = [read_f64_be(data, p + int(i) * 8) for i in ids]
        coords_file[:, axis_idx] = vals
    coords = coords_file[:, list(layout["perm"])]
    return coords, layout


def _parse_ls_links_layout(data) -> Optional[dict]:
    sec_start = find_section(data, "LS_Links")
    if sec_start < 0:
        return None
    sec_end = section_end(data, sec_start)
    blocks = [(p, bc) for p, bc in iter_data_blocks(data, sec_start, sec_end) if bc > 0]
    if not blocks:
        return None

    from collections import Counter
    common = Counter([bc for _, bc in blocks]).most_common()
    n_faces_block_size = None
    for size, count in common:
        if count >= 3 and size % 4 == 0 and size >= 4:
            n_faces_block_size = size
            break
    if n_faces_block_size is None:
        return None
    n_faces = n_faces_block_size // 4
    triples = [b for b in blocks if b[1] == n_faces_block_size][:3]
    if len(triples) < 3:
        return None
    owner_p, _ = triples[0]
    neigh_p, _ = triples[1]
    npe_p, _ = triples[2]
    owner = np.frombuffer(data, dtype=">u4", count=n_faces, offset=owner_p).astype(np.int64)
    neigh_raw = np.frombuffer(data, dtype=">u4", count=n_faces, offset=neigh_p)
    neighbor = neigh_raw.astype(np.int64)
    neighbor[neigh_raw == 0xFFFFFFFF] = -1
    npe = np.frombuffer(data, dtype=">u4", count=n_faces, offset=npe_p).astype(np.int64)
    conn_total_expected = int(npe.sum())
    conn_block = None
    for p, bc in blocks:
        if (p, bc) in triples or bc % 4 != 0:
            continue
        if bc // 4 == conn_total_expected:
            conn_block = (p, bc)
            break
    if conn_block is None:
        for p, bc in blocks:
            if (p, bc) in triples or bc % 4 != 0 or bc < 12:
                continue
            if conn_block is None or bc > conn_block[1]:
                conn_block = (p, bc)
    if conn_block is None:
        return None
    conn_p, conn_bc = conn_block
    n_cells = int(max(
        int(owner.max()) + 1 if owner.size else 0,
        int(neighbor[neighbor >= 0].max()) + 1 if (neighbor >= 0).any() else 0,
    ))
    return {
        "n_faces": n_faces,
        "n_cells": n_cells,
        "owner": owner,
        "neighbor": neighbor,
        "npe": npe,
        "conn_p": conn_p,
        "conn_entries_expected": conn_total_expected,
        "conn_entries_available": conn_bc // 4,
    }


def build_mesh_preview(
    data,
    selected_face_ids: Optional["np.ndarray"] = None,
    selected_cell_ids: Optional["np.ndarray"] = None,
    max_faces: int = 12000,
) -> Optional[dict]:
    """Build a small polygon preview for interactive 3D display."""
    links = _parse_ls_links_layout(data)
    if not links:
        return None
    n_faces = int(links["n_faces"])
    owner = links["owner"]
    neighbor = links["neighbor"]
    npe = links["npe"]

    face_sets: list["np.ndarray"] = []
    if selected_face_ids is not None and selected_face_ids.size:
        face_sets.append(np.asarray(selected_face_ids, dtype=np.int64))
    if selected_cell_ids is not None and selected_cell_ids.size:
        cells = np.asarray(selected_cell_ids, dtype=np.int64)
        owner_in = np.isin(owner, cells)
        neigh_in = np.isin(neighbor, cells)
        face_sets.append(np.flatnonzero(owner_in | neigh_in))

    selection_active = bool(face_sets)
    if face_sets:
        face_ids = np.unique(np.concatenate(face_sets))
        face_ids = face_ids[(face_ids >= 0) & (face_ids < n_faces)]
    else:
        face_ids = np.arange(min(n_faces, max_faces), dtype=np.int64)

    if face_ids.size > max_faces:
        idx = np.linspace(0, face_ids.size - 1, max_faces, dtype=np.int64)
        face_ids = face_ids[idx]

    if face_ids.size == 0:
        return {
            "faces": [],
            "face_ids": np.empty(0, dtype=np.int64),
            "selection_active": selection_active,
            "summary": "No faces matched the selected region.",
        }

    face_offsets = np.empty(n_faces + 1, dtype=np.int64)
    face_offsets[0] = 0
    np.cumsum(npe, out=face_offsets[1:])
    conn_available = int(links["conn_entries_available"])
    valid = face_offsets[face_ids + 1] <= conn_available
    face_ids = face_ids[valid]
    if face_ids.size == 0:
        return {
            "faces": [],
            "face_ids": np.empty(0, dtype=np.int64),
            "selection_active": selection_active,
            "summary": "Selected faces are outside the available preview connectivity block.",
        }

    need_conn = int(face_offsets[face_ids[-1] + 1])
    conn = np.frombuffer(
        data, dtype=">u4", count=need_conn, offset=int(links["conn_p"]),
    ).astype(np.int64)

    face_node_ids: list["np.ndarray"] = []
    for fid in face_ids:
        lo = int(face_offsets[fid])
        hi = int(face_offsets[fid + 1])
        face_node_ids.append(conn[lo:hi])
    unique_vertices = np.unique(np.concatenate(face_node_ids)) if face_node_ids else np.empty(0, dtype=np.int64)
    coords, node_layout = _read_vertices_by_id(data, unique_vertices)
    vertex_lookup = {int(vid): coords[i] for i, vid in enumerate(unique_vertices)}
    polygons = [
        np.array([vertex_lookup[int(v)] for v in nodes if int(v) in vertex_lookup], dtype=float)
        for nodes in face_node_ids
    ]
    polygons = [poly for poly in polygons if poly.shape[0] >= 3]
    return {
        "faces": polygons,
        "face_ids": face_ids,
        "selection_active": selection_active,
        "n_faces": n_faces,
        "n_cells": int(links["n_cells"]),
        "n_vertices": int(node_layout.get("n_vertices", 0)) if node_layout else 0,
        "dialect": node_layout.get("dialect", "") if node_layout else "",
        "summary": f"{face_ids.size} faces",
    }


@dataclass
class GphNode:
    """A node in the GPH tree - can represent a section or a data item."""

    name: str
    offset: int
    size: int
    data_type: str  # "raw", "I4", "R4", "C1", "I4[]", "R4[]"
    value: Any = None  # parsed value for display/edit
    raw: bytes = b""
    children: list = field(default_factory=list)
    modified: bool = False
    parent: Optional["GphNode"] = None
    metadata: dict = field(default_factory=dict)

    def get_raw(self, doc: Optional["GphDocument"] = None) -> bytes:
        """Return raw bytes - from doc if provided (for modified data), else cached."""
        if doc is not None and hasattr(doc, "_raw_data"):
            return bytes(doc._raw_data[self.offset : self.offset + self.size])
        return self.raw

    def set_value(self, val: Any, raw: bytes) -> None:
        self.value = val
        self.raw = raw
        self.modified = True


class GphDocument:
    """In-memory GPH document with edit tracking."""

    def __init__(self):
        self.filepath: Optional[str] = None
        self._raw_data: Any = bytearray()
        self._file_handle: Any = None
        self._mmap_mode: bool = False
        self.root: Optional[GphNode] = None
        self._patches: list[tuple[int, bytes]] = []  # (offset, new_bytes)

    def close(self) -> None:
        if self._mmap_mode and self._raw_data is not None:
            self._raw_data.close()
        if self._file_handle is not None:
            self._file_handle.close()
        self._raw_data = bytearray()
        self._file_handle = None
        self._mmap_mode = False

    def load(self, filepath: str) -> bool:
        try:
            self.close()
            size = Path(filepath).stat().st_size
            if size > _LARGE_GPH_BYTES:
                import mmap
                self._file_handle = open(filepath, "rb")
                self._raw_data = mmap.mmap(
                    self._file_handle.fileno(), 0, access=mmap.ACCESS_READ,
                )
                self._mmap_mode = True
            else:
                with open(filepath, "rb") as f:
                    self._raw_data = bytearray(f.read())
            self.filepath = filepath
            self._patches = []
            self.root = self._parse()
            return True
        except Exception:
            self.close()
            return False

    def _parse(self) -> GphNode:
        data = self._raw_data
        root = GphNode("GPH File", 0, len(data), "raw", children=[])

        # ── Dynamically locate every known named section ────────────────────
        #
        # Each named section is preceded by ``[I4=32]`` and contains a
        # 32-byte ASCII (space-padded) label.  We scan for those labels to
        # build a section layout that adapts to any file size (the legacy
        # ``box.gph`` and the new ``box_ansa.gph`` differ by ~13 kB).
        candidate_names = _SECTION_BOUNDARY_NAMES
        found = []  # list of (offset, name)
        for name in candidate_names:
            padded = name.ljust(32).encode("ascii")
            idx = data.find(padded)
            if idx >= 4 and read_i32_be(data, idx - 4) == 32:
                found.append((idx - 4, name))
        found.sort(key=lambda x: x[0])

        # The fixed file-header sits before the first named section.
        first_off = found[0][0] if found else len(data)
        sections = [(0, first_off, "file_header", "Header (CRDL-FLD + dims)")]
        for i, (off, name) in enumerate(found):
            end = found[i + 1][0] if i + 1 < len(found) else len(data)
            sections.append((off, end, name, ""))

        for start, end, name, desc in sections:
            raw = data[start:end]
            node = self._create_node(name, start, raw, desc)
            self._add_binary_children(node)
            node.parent = root
            root.children.append(node)

        return root

    def _create_node(self, name: str, offset: int, raw: bytes, desc: str) -> GphNode:
        if name == "file_header" and len(raw) >= 12:
            fid = raw[4:12].decode("ascii", errors="replace")
            v1, v2, v3 = read_i32_be(raw, 12), read_i32_be(raw, 16), read_i32_be(raw, 20)
            val = f"{fid}  dims=({v1},{v2},{v3})"
            return GphNode(name, offset, len(raw), "header", value=val, raw=raw, children=[])

        if name == "FileRevision" and len(raw) >= 68:
            v = read_i32_be(raw, 64) if len(raw) > 64 else 0
            return GphNode(name, offset, len(raw), "I4", value=v, raw=raw, children=[])

        if name == "Application" and len(raw) >= 72:
            s = raw[64:72].decode("ascii", errors="replace").strip()
            return GphNode(name, offset, len(raw), "C1[8]", value=s, raw=raw, children=[])

        if name == "Dimension" and len(raw) >= 68:
            v = read_i32_be(raw, 64) if len(raw) > 64 else 0
            return GphNode(name, offset, len(raw), "I4", value=v, raw=raw, children=[])

        file_data = self._raw_data

        if name == "LS_CvolIdOfElements":
            cvol_arr = parse_ls_cvol_ids(file_data)
            if cvol_arr is None:
                return GphNode(name, offset, len(raw), "I4[]",
                               value=None, raw=raw, children=[])
            unique_cv = sorted({int(x) for x in np.unique(cvol_arr[:min(len(cvol_arr), 1_000_000)])})
            summary = (
                f"I4[{len(cvol_arr)}] cvol_ids={unique_cv[:12]}"
                f"{'...' if len(unique_cv) > 12 else ''}"
            )
            preview = cvol_arr[:1000].tolist()
            return GphNode(name, offset, len(raw), summary,
                           value=preview, raw=raw, children=[])

        if name == "LS_Nodes":
            sample, dialect, n_vertices = parse_ls_nodes_vertices(file_data)
            if n_vertices:
                elem = "R4" if "float32" in dialect else "R8"
                dtype = f"{elem}[{n_vertices},3] ({dialect})"
                return GphNode(name, offset, len(raw), dtype,
                               value=sample, raw=raw, children=[])
            return GphNode(name, offset, len(raw), "R4/R8[]", value=None, raw=raw, children=[])

        if name == "LS_Links":
            summary = parse_ls_links_summary(file_data)
            if summary:
                val = (
                    f"faces={summary['n_faces']} cells={summary['n_cells']} "
                    f"BC={summary['boundary_faces']} "
                    f"npe=[{summary['npe_min']}..{summary['npe_max']}]"
                    + (" polyhedral" if summary["polyhedral"] else "")
                    + (f" conn_split×{summary['conn_chunks']}"
                       if summary.get("conn_split") else "")
                    + (" conn_INCOMPLETE" if not summary.get("conn_complete", True) else "")
                )
                return GphNode(name, offset, len(raw), "topology", value=val,
                               raw=raw, children=[])
            arr = self._collect_data_blocks_i4(raw)
            type_str = f"I4[{len(arr)}]" if arr is not None else "I4[]"
            return GphNode(name, offset, len(raw), type_str,
                           value=arr, raw=raw, children=[])

        if name == "LS_Parts":
            cvol_arr_for_parts = parse_ls_cvol_ids(file_data)
            parts = parse_ls_parts(file_data, cvol_id=cvol_arr_for_parts)
            val = [f"{p} (cvol={format_part_cvol_spec(cv)})" for p, cv in parts]
            return GphNode(name, offset, len(raw),
                           f"parts[{len(parts)}]", value=val, raw=raw, children=[])

        if name == "LS_VolumeRegions":
            regions = parse_ls_string_list(file_data, "LS_VolumeRegions")
            return GphNode(name, offset, len(raw),
                           f"regions[{len(regions)}]", value=regions, raw=raw, children=[])

        if name == "LS_SurfaceRegions":
            regions = parse_ls_surface_regions_summary(file_data)
            val = [f"{n} ({nf} faces)" for n, nf in regions]
            return GphNode(name, offset, len(raw),
                           f"surf_regions[{len(regions)}]", value=val, raw=raw, children=[])

        if name == "LS_Assemblies":
            asm = parse_ls_assemblies_summary(file_data)
            lines = []
            if asm["root_empty_prefix"]:
                lines.append(f"root_empty_prefix={asm['root_empty_prefix']}")
            for pname, path in asm["part_paths"].items():
                lines.append(f"{pname} -> {path or '(root)'}")
            return GphNode(name, offset, len(raw), "assembly_xml",
                           value="\n".join(lines) if lines else "(no XML)",
                           raw=raw, children=[])

        return GphNode(name, offset, len(raw), "raw", value=desc, raw=raw, children=[])

    # ── Helpers ──────────────────────────────────────────────────────────────

    def _add_binary_children(self, node: GphNode) -> None:
        """Attach record-level children so the GUI can expand binary layout."""
        raw = node.raw
        if node.name == "file_header":
            self._add_file_header_children(node, raw)
            return
        if len(raw) >= 40:
            label_len = read_i32_be(raw, 0)
            self._append_child(
                node, "section_name_length", 0, 4, "I4", label_len,
            )
            label = raw[4:36].decode("ascii", errors="replace").rstrip()
            self._append_child(
                node, "section_name", 4, 32, "C1[32]", label,
            )
            self._scan_section_records(node, raw, 40)

    def _add_file_header_children(self, node: GphNode, raw: bytes) -> None:
        if len(raw) >= 4:
            n = read_i32_be(raw, 0)
            self._append_child(node, "format_id_length", 0, 4, "I4", n)
        if len(raw) >= 12:
            self._append_child(
                node, "format_id", 4, min(8, len(raw) - 4), "C1",
                raw[4:12].decode("ascii", errors="replace"),
            )
        dims = [("dim0", 12), ("dim1", 16), ("dim2", 20)]
        for name, rel in dims:
            if len(raw) >= rel + 4:
                self._append_child(node, name, rel, 4, "I4", read_i32_be(raw, rel))
        if len(raw) > 24:
            self._append_child(node, "header_payload", 24, len(raw) - 24, "raw")

    def _scan_section_records(self, node: GphNode, raw: bytes, pos: int) -> None:
        block_idx = 0
        unknown_start: Optional[int] = None

        def flush_unknown(until: int) -> None:
            nonlocal unknown_start
            if unknown_start is not None and until > unknown_start:
                self._append_child(
                    node, "unclassified", unknown_start, until - unknown_start, "raw",
                )
            unknown_start = None

        while pos + 8 <= len(raw):
            if read_i32_be(raw, pos) != 12:
                if unknown_start is None:
                    unknown_start = pos
                pos += 4
                continue
            v = read_i32_be(raw, pos + 4)
            if v in (1, 4, 8) and pos + 16 <= len(raw):
                d0 = read_i32_be(raw, pos + 8)
                d1 = read_i32_be(raw, pos + 12)
                if 0 <= d0 < 1_000_000_000 and 0 <= d1 < 1_000_000_000:
                    flush_unknown(pos)
                    val = f"type={v}, dims=({d0}, {d1})"
                    self._append_child(node, "descriptor", pos, 16, "descriptor", val)
                    pos += 16
                    continue
            bc = v
            if bc > 0 and pos + 8 + bc + 4 <= len(raw):
                payload_end = pos + 8 + bc
                if read_i32_be(raw, payload_end) == bc:
                    flush_unknown(pos)
                    block_idx += 1
                    block = self._append_child(
                        node,
                        f"data_block[{block_idx}]",
                        pos,
                        8 + bc + 4,
                        f"block[{bc} B]",
                        f"payload={bc} bytes",
                    )
                    self._append_child(block, "block_marker", 0, 4, "I4", 12)
                    self._append_child(block, "byte_count", 4, 4, "I4", bc)
                    self._append_child(block, "payload", 8, bc, "payload")
                    self._append_child(block, "byte_count_trailer", 8 + bc, 4, "I4", bc)
                    pos = payload_end + 4
                    continue
            if unknown_start is None:
                unknown_start = pos
            pos += 4
        flush_unknown(len(raw))

    @staticmethod
    def _preview_value(raw: bytes, data_type: str):
        if data_type.startswith("C1"):
            return raw.decode("ascii", errors="replace").strip("\x00").rstrip()
        if data_type == "payload":
            if raw and all(b == 0 or 32 <= b < 127 for b in raw[: min(len(raw), 128)]):
                return raw[:128].decode("ascii", errors="replace").strip("\x00").rstrip()
            return None
        return None

    def _append_child(
        self,
        parent: GphNode,
        name: str,
        rel_offset: int,
        size: int,
        data_type: str,
        value: Any = None,
    ) -> GphNode:
        raw = parent.raw[rel_offset : rel_offset + size]
        child = GphNode(
            name,
            parent.offset + rel_offset,
            size,
            data_type,
            self._preview_value(raw, data_type) if value is None else value,
            raw=raw,
            children=[],
            parent=parent,
        )
        parent.children.append(child)
        return child

    @staticmethod
    def _read_i4_array_after_label(raw: bytes) -> Optional[list[int]]:
        """Return the first I4 array stored in a section.

        Skips the 40-byte label header, any metadata descriptors of the form
        ``[12, type, dim0, dim1]`` and the data-block header (8 or 12 bytes),
        then reads ``byte_count // 4`` big-endian int32 values.
        """
        pos = 40
        n = 0
        # Find a descriptor [12, 4, n, 1] giving the array length.
        while pos + 16 <= len(raw):
            if read_i32_be(raw, pos) != 12:
                break
            tc = read_i32_be(raw, pos + 4)
            d0 = read_i32_be(raw, pos + 8)
            d1 = read_i32_be(raw, pos + 12)
            pos += 16
            if tc == 4 and d1 == 1 and d0 > 1:
                n = d0
                break
        if n == 0:
            return None
        # Locate the next data header: [12, byte_count] (8 bytes) — its
        # ``byte_count`` should equal ``n * 4``.
        expected_bc = n * 4
        # Tolerate up to ~16 bytes of additional padding/descriptor.
        for shift in (0, 4, 8, 12, 16):
            p = pos + shift
            if p + 8 + expected_bc > len(raw):
                continue
            if read_i32_be(raw, p) == 12 and read_i32_be(raw, p + 4) == expected_bc:
                p += 8
                return [read_i32_be(raw, p + i * 4) for i in range(n)]
        return None

    @staticmethod
    def _collect_data_blocks_i4(raw: bytes) -> Optional[list[int]]:
        """Walk a section, concatenating every I4 payload (preview, capped)."""
        pos = 40
        out: list[int] = []
        max_preview = 1000
        while pos + 8 <= len(raw) and len(out) < max_preview:
            if read_i32_be(raw, pos) != 12:
                pos += 4
                continue
            v = read_i32_be(raw, pos + 4)
            # Descriptor [12, type in (4,8), dim0, dim1]?
            if v in (4, 8) and pos + 16 <= len(raw):
                d0 = read_i32_be(raw, pos + 8)
                d1 = read_i32_be(raw, pos + 12)
                if 0 < d0 < 10_000_000 and 0 < d1 < 10_000_000:
                    pos += 16
                    continue
            # Data header [12, byte_count]?
            bc = v
            if bc <= 0 or pos + 8 + bc + 4 > len(raw):
                pos += 4
                continue
            payload_end = pos + 8 + bc
            if read_i32_be(raw, payload_end) != bc:
                pos += 4
                continue
            for i in range(bc // 4):
                if len(out) >= max_preview:
                    break
                out.append(read_i32_be(raw, pos + 8 + i * 4))
            pos = payload_end + 4
        return out if out else None

    def apply_patch(self, offset: int, new_bytes: bytes) -> None:
        """Record a modification for save."""
        if self._mmap_mode:
            raise RuntimeError("Cannot edit memory-mapped GPH files in place")
        self._patches.append((offset, new_bytes))
        if offset + len(new_bytes) <= len(self._raw_data):
            self._raw_data[offset : offset + len(new_bytes)] = new_bytes

    def save(self, filepath: Optional[str] = None) -> bool:
        path = filepath or self.filepath
        if not path:
            return False
        try:
            # Apply all patches (already in _raw_data if we tracked in-memory)
            with open(path, "wb") as f:
                f.write(self._raw_data)
            self.filepath = path
            self._patches = []
            self._clear_modified_flag(self.root)
            return True
        except Exception:
            return False

    def _clear_modified_flag(self, node: Optional[GphNode]) -> None:
        if node is None:
            return
        node.modified = False
        for c in node.children:
            self._clear_modified_flag(c)

    def get_data_at(self, offset: int, size: int) -> bytes:
        """Get current data (including modifications) at offset."""
        return bytes(self._raw_data[offset : offset + size])

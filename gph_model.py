#!/usr/bin/env python3
"""
GPH data model - parse binary to editable tree, support partial save.
"""

import struct
from dataclasses import dataclass, field
from typing import Any, Optional


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


def parse_ls_nodes_vertices(data: bytes) -> tuple[Optional[list[tuple[float, float, float]]], str, int]:
    """Parse LS_Nodes → (vertices, dialect_label, n_vertices)."""
    sec_start = find_section(data, "LS_Nodes")
    if sec_start < 0:
        return None, "", 0
    sec_end = section_end(data, sec_start)
    blocks = list(iter_data_blocks(data, sec_start, sec_end))
    f64_blocks = [(p, bc) for p, bc in blocks if bc >= 8 and bc % 8 == 0]
    if len(f64_blocks) < 3:
        return None, "", 0
    sizes = [bc for _, bc in f64_blocks]
    target = max(set(sizes), key=sizes.count)
    trio = [(p, bc) for p, bc in f64_blocks if bc == target][:3]
    if len(trio) < 3:
        return None, "", 0
    n_vertices = trio[0][1] // 8

    def _decode(reader):
        axes: list[list[float]] = []
        for payload_start, _ in trio:
            vals = [reader(data, payload_start + i * 8) for i in range(n_vertices)]
            axes.append(vals)
        return axes

    def _score(axes):
        score = 0.0
        for ax in axes:
            for v in ax:
                if not (v == v):  # NaN
                    score += 1e30
                    continue
                a = abs(v)
                if a != 0.0 and (a > 1e6 or a < 1e-30):
                    score += a + (1.0 / max(a, 1e-300))
                else:
                    score += a
        return score

    axes_be = _decode(read_f64_be)
    axes_wr = _decode(read_f64_wr)
    if _score(axes_be) <= _score(axes_wr):
        axes, dialect = axes_be, "standard BE float64"
        col_perm = (0, 1, 2)
    else:
        axes, dialect = axes_wr, "word-reversed float64"
        col_perm = (0, 2, 1)
    verts = [
        (axes[col_perm[0]][i], axes[col_perm[1]][i], axes[col_perm[2]][i])
        for i in range(n_vertices)
    ]
    return verts, dialect, n_vertices


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
    npe_p, _ = triples[2]
    npe = [read_i32_be(data, npe_p + i * 4) for i in range(n_faces)]
    conn_total = sum(npe)
    boundary = 0
    BNDRY = 0xFFFFFFFF
    neigh_p, _ = triples[1]
    for i in range(n_faces):
        raw = data[neigh_p + i * 4 : neigh_p + i * 4 + 4]
        if int.from_bytes(raw, "big") == BNDRY:
            boundary += 1
    owner_p, _ = triples[0]
    max_own = max(read_i32_be(data, owner_p + i * 4) for i in range(n_faces))
    n_cells = max_own + 1
    return {
        "n_faces": n_faces,
        "n_cells": n_cells,
        "boundary_faces": boundary,
        "npe_min": min(npe),
        "npe_max": max(npe),
        "conn_entries": conn_total,
        "polyhedral": max(npe) > 3 if npe else False,
    }


def parse_ls_parts(data: bytes) -> list[tuple[str, int]]:
    """Parse LS_Parts → [(part_name, cvol_id), ...] (see gph2cgns)."""
    sec_start = find_section(data, "LS_Parts")
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    name_blocks: list[tuple[str, int, int]] = []
    for p, bc in iter_data_blocks(data, sec_start, sec_end):
        if bc != 255:
            continue
        raw = data[p : p + bc]
        if not all(b == 0 or 32 <= b < 127 for b in raw):
            continue
        name = raw.decode("ascii", errors="replace").strip("\x00").rstrip()
        if not name or not any(c.isalpha() for c in name):
            continue
        name_blocks.append((name, p - 8, p + bc + 4))
    out: list[tuple[str, int]] = []
    for i, (name, _, after_trailer) in enumerate(name_blocks):
        scan_end = name_blocks[i + 1][1] if i + 1 < len(name_blocks) else sec_end
        cvol_id: Optional[int] = None
        pos = after_trailer
        while pos + 16 <= scan_end:
            if (read_i32_be(data, pos) == 12
                    and read_i32_be(data, pos + 4) == 4
                    and read_i32_be(data, pos + 12) == 4):
                cvol_id = read_i32_be(data, pos + 8)
            pos += 4
        if cvol_id is not None:
            out.append((name, int(cvol_id)))
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


def parse_ls_surface_regions_summary(data: bytes) -> list[tuple[str, int]]:
    sec_start = find_section(data, "LS_SurfaceRegions")
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    blocks = list(iter_data_blocks(data, sec_start, sec_end))
    out: list[tuple[str, int]] = []
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
            out.append((name, bc_i // 4))
            i += 3
        else:
            i += 1
    return out


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
        self._raw_data: bytearray = bytearray()
        self.root: Optional[GphNode] = None
        self._patches: list[tuple[int, bytes]] = []  # (offset, new_bytes)

    def load(self, filepath: str) -> bool:
        try:
            with open(filepath, "rb") as f:
                self._raw_data = bytearray(f.read())
            self.filepath = filepath
            self._patches = []
            self.root = self._parse()
            return True
        except Exception:
            return False

    def _parse(self) -> GphNode:
        data = bytes(self._raw_data)
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

        file_data = bytes(self._raw_data)

        if name == "LS_CvolIdOfElements":
            arr = self._read_i4_array_after_label(raw)
            unique_cv = sorted(set(arr)) if arr else []
            summary = (
                f"I4[{len(arr)}] cvol_ids={unique_cv[:12]}"
                f"{'...' if len(unique_cv) > 12 else ''}"
                if arr is not None else "I4[]"
            )
            return GphNode(name, offset, len(raw), summary,
                           value=arr, raw=raw, children=[])

        if name == "LS_Nodes":
            verts, dialect, n_vertices = parse_ls_nodes_vertices(file_data)
            if verts and n_vertices:
                dtype = f"R8[{n_vertices},3] ({dialect})"
                return GphNode(name, offset, len(raw), dtype,
                               value=verts, raw=raw, children=[])
            return GphNode(name, offset, len(raw), "R8[]", value=None, raw=raw, children=[])

        if name == "LS_Links":
            summary = parse_ls_links_summary(file_data)
            if summary:
                val = (
                    f"faces={summary['n_faces']} cells={summary['n_cells']} "
                    f"BC={summary['boundary_faces']} "
                    f"npe=[{summary['npe_min']}..{summary['npe_max']}]"
                    + (" polyhedral" if summary["polyhedral"] else "")
                )
                return GphNode(name, offset, len(raw), "topology", value=val,
                               raw=raw, children=[])
            arr = self._collect_data_blocks_i4(raw)
            type_str = f"I4[{len(arr)}]" if arr is not None else "I4[]"
            return GphNode(name, offset, len(raw), type_str,
                           value=arr, raw=raw, children=[])

        if name == "LS_Parts":
            parts = parse_ls_parts(file_data)
            val = [f"{p} (cvol_id={cv})" for p, cv in parts]
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

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
        candidate_names = [
            "FileRevision",
            "Application",
            "ApplicationVersion",
            "ReleaseDate",
            "GridType",
            "Dimension",
            "Bias",
            "Date",
            "Comments",
            "Cycle",
            "Unused",
            "Encoding",
            "HeaderDataEnd",
            "OverlapStart_0",
            "LS_CvolIdOfElements",
            "LS_Links",
            "LS_Nodes",
            "LS_SurfaceRegions",
            "Element_InformationFlag",
            "OverlapEnd",
        ]
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

        if name == "LS_CvolIdOfElements":
            # Locate the [12, type=4, dim0=n_elements, dim1=1] descriptor by
            # scanning forward from the label header; the data follows after a
            # 12-byte block header (legacy) or 8-byte block header (modern).
            arr = self._read_i4_array_after_label(raw)
            type_str = f"I4[{len(arr)}]" if arr is not None else "I4[]"
            return GphNode(name, offset, len(raw), type_str,
                           value=arr, raw=raw, children=[])

        if name == "LS_Nodes":
            # Coordinates are stored as three float64 axis blocks.  All
            # observed GPH files (legacy box.gph as well as the modern
            # ANSA / scFLOW box_ansa.gph) use the same on-disk encoding:
            #
            #   [16B descriptor: 12 / 8 / n_verts / 1]
            #   [8B  data header: 12 / byte_count]
            #   [n_verts × 8B] standard big-endian IEEE-754 float64
            #   [4B trailer: byte_count]
            #
            # axis order in the file is X, Y, Z.  A word-reversed reading with
            # a 12-byte data header is kept as defensive fallback for any
            # hypothetical GPH variant that might use the legacy mid-endian
            # encoding (none has been observed in practice — the historical
            # "word-reversed" theory was an off-by-4 alignment artifact).
            pos = 40  # skip 40B label-header
            n_vertices = 0
            while pos + 16 <= len(raw):
                if read_i32_be(raw, pos) != 12:
                    break
                tc = read_i32_be(raw, pos + 4)
                d0 = read_i32_be(raw, pos + 8)
                d1 = read_i32_be(raw, pos + 12)
                pos += 16
                if tc == 8 and d1 == 1 and d0 > 0:
                    n_vertices = d0
                    break

            def _read_axis_blocks(data_header_size: int, reader) -> list[list[float]]:
                """Return [X_vals, ?, ?] for the three coordinate blocks."""
                p = pos
                axes: list[list[float]] = []
                for _ in range(3):
                    if p + data_header_size > len(raw):
                        break
                    p += data_header_size
                    vals = []
                    for _ in range(n_vertices):
                        if p + 8 > len(raw):
                            break
                        vals.append(reader(raw, p))
                        p += 8
                    axes.append(vals)
                    # Skip the 4-byte trailing sentinel (8-byte-header dialect)
                    # or the 16-byte inter-block descriptor (12-byte-header dialect).
                    if data_header_size == 8 and p + 4 <= len(raw):
                        p += 4
                    if p + 16 <= len(raw) and read_i32_be(raw, p) == 12:
                        p += 16
                return axes

            if n_vertices > 0:
                # Try the modern 8-byte-header / standard BE float64 dialect first.
                axes_new = _read_axis_blocks(8, read_f64_be)
                axes_old = _read_axis_blocks(12, read_f64_wr)
                if (len(axes_new) == 3
                        and all(_looks_like_coords(a) for a in axes_new)):
                    # ANSA dialect: file order X, Y, Z.
                    axis_vals, swap_zy = axes_new, False
                else:
                    # Legacy SCTpre dialect: file order X, Z, Y → swap.
                    axis_vals, swap_zy = axes_old, True
                if len(axis_vals) == 3:
                    if swap_zy:
                        verts = [
                            (axis_vals[0][i], axis_vals[2][i], axis_vals[1][i])
                            for i in range(n_vertices)
                        ]
                    else:
                        verts = [
                            (axis_vals[0][i], axis_vals[1][i], axis_vals[2][i])
                            for i in range(n_vertices)
                        ]
                    return GphNode(
                        name, offset, len(raw), "R8[n,3]",
                        value=verts, raw=raw, children=[],
                    )
            return GphNode(name, offset, len(raw), "R8[]", value=None, raw=raw, children=[])

        if name == "LS_Links":
            # Concatenate every data block (owner / neighbor / npe / [face_type]
            # / conn) by scanning the section with the same data-block pattern
            # used by gph2cgns.  Returns a flat I4 array preview for the
            # viewer's Data tab.
            arr = self._collect_data_blocks_i4(raw)
            type_str = f"I4[{len(arr)}]" if arr is not None else "I4[]"
            return GphNode(name, offset, len(raw), type_str,
                           value=arr, raw=raw, children=[])

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

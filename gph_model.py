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

        # Known section layout from reverse engineering
        sections = [
            (0x0000, 0x001C, "file_header", "Header (CRDL-FLD + dims)"),
            (0x001C, 0x0078, "FileRevision", "I4"),
            (0x0078, 0x00D8, "Application", "C1[8]"),
            (0x00D8, 0x0134, "ApplicationVersion", "I4"),
            (0x0134, 0x0190, "ReleaseDate", "string"),
            (0x0190, 0x01EC, "GridType", "string"),
            (0x01EC, 0x0248, "Dimension", "I4"),
            (0x0248, 0x02A4, "Bias", "I4"),
            (0x02A4, 0x0300, "Date", "string"),
            (0x0300, 0x03A8, "Comments", "string"),
            (0x03A8, 0x04E0, "Cycle", "I4 + Unit"),
            (0x04E0, 0x0560, "Unused", "reserved"),
            (0x0560, 0x05D8, "Encoding", "string"),
            (0x05D8, 0x0600, "HeaderDataEnd", "marker"),
            (0x0600, 0x0628, "OverlapStart_0", "marker"),
            (0x0628, 0x08DC, "LS_CvolIdOfElements", "I4[135]"),
            (0x08DC, 0x26B0, "LS_Links", "I4[]"),
            (0x26B0, 0x2C60, "LS_Nodes", "R8[n,3]"),
            (0x2C60, 0x42B0, "LS_SurfaceRegions", "raw"),
            (0x42B0, 0x45A4, "Element_InformationFlag", "raw"),
            (0x45A4, len(data), "OverlapEnd", "trailer"),
        ]

        for start, end, name, desc in sections:
            raw = data[start:end]
            node = self._create_node(name, start, raw, desc)
            node.parent = root
            root.children.append(node)
            root.size += 0

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
            # Data at ~0x69C, 135 x I4
            data_start = 0x069C - offset
            if data_start >= 0 and data_start + 540 <= len(raw):
                arr = [
                    read_i32_be(raw, data_start + i * 4)
                    for i in range(min(135, (len(raw) - data_start) // 4))
                ]
                return GphNode(
                    name, offset, len(raw), "I4[135]",
                    value=arr, raw=raw, children=[],
                )
            return GphNode(name, offset, len(raw), "I4[]", value=None, raw=raw, children=[])

        if name == "LS_Nodes":
            # Coordinates are stored as three float64 axis blocks.  Two on-disk
            # dialects are supported (auto-detected at parse time):
            #
            #   Legacy SCTpre  – word-reversed float64,  12-byte data header,
            #                    axis order X, Z, Y.
            #   New ANSA/scFLOW – standard big-endian float64, 8-byte data
            #                    header, axis order X, Y, Z.
            #
            # Each block is preceded by a 16-byte descriptor [12, 8, n, 1] and
            # terminated by a 4-byte trailing length sentinel.
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
            data_start = 0x09C0 - offset
            if data_start >= 0 and data_start < len(raw):
                n = (len(raw) - data_start) // 4
                arr = [read_i32_be(raw, data_start + i * 4) for i in range(min(n, 500))]
                return GphNode(name, offset, len(raw), "I4[]", value=arr, raw=raw, children=[])
            return GphNode(name, offset, len(raw), "I4[]", value=None, raw=raw, children=[])

        return GphNode(name, offset, len(raw), "raw", value=desc, raw=raw, children=[])

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

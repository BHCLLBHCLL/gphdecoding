#!/usr/bin/env python3
"""Unit tests for LS_Links conn continuation parsing."""
import struct
import unittest

import numpy as np

from gph_model import _CONN_CHUNK_BYTES, _read_conn_continuations


def _i4(n: int) -> bytes:
    return struct.pack(">i", n)


def _u4_array(values) -> bytes:
    return b"".join(struct.pack(">I", v) for v in values)


class ConnContinuationTests(unittest.TestCase):
    def test_single_header_short_tail(self):
        """Two-chunk layout: [I4=bc][payload] (v4-style)."""
        payload = _u4_array([10, 11, 12, 13])
        need = len(payload)
        buf = bytearray(_i4(need) + payload)
        got, _, n = _read_conn_continuations(buf, 0, len(buf), 0, 4)
        self.assertEqual(got, 4)
        self.assertEqual(n, 1)

    def test_double_header_short_tail(self):
        """Three-chunk layout: final segment [1GiB marker][need_bytes][payload]."""
        tail_vals = [100, 101, 102, 103, 104]
        tail_bytes = len(tail_vals) * 4
        buf = bytearray(_i4(_CONN_CHUNK_BYTES) + _i4(tail_bytes) + _u4_array(tail_vals))
        got_start = _CONN_CHUNK_BYTES // 4 * 2
        expected = got_start + len(tail_vals)
        parts: list = []
        got, _, n = _read_conn_continuations(
            buf, 0, len(buf), got_start, expected, parts,
        )
        self.assertEqual(n, 1)
        self.assertEqual(got, expected)
        conn = np.concatenate(parts)
        self.assertEqual(conn.tolist(), tail_vals)

    def test_double_header_does_not_leak_byte_count(self):
        """Regression: byte_count must not appear as a vertex index."""
        tail_vals = [1, 2, 3]
        tail_bytes = len(tail_vals) * 4
        buf = bytearray(_i4(_CONN_CHUNK_BYTES) + _i4(tail_bytes) + _u4_array(tail_vals))
        got_start = _CONN_CHUNK_BYTES // 4 * 2
        expected = got_start + len(tail_vals)
        parts: list = []
        _read_conn_continuations(buf, 0, len(buf), got_start, expected, parts)
        conn = np.concatenate(parts)
        self.assertNotIn(tail_bytes, conn.tolist())


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""
Traverse tests/*.gph and compare volume-zone cell counts against matching
*_orig.cgns references (same stem, e.g. box_ansa.gph ↔ box_ansa_orig.cgns).

Run from repo root:
    python tests/test_volume_zone_cells.py
    python tests/test_volume_zone_cells.py -v
    python tests/test_volume_zone_cells.py -v tests/box_ansa.gph
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    import h5py
    import numpy as np
except ImportError:
    print("Error: h5py and numpy are required. Install with: pip install h5py numpy")
    sys.exit(1)

from gph2cgns import _build_zone_plan, parse_gph_mesh
from gph_model import (
    _ls_parts_name_blocks,
    _scan_cvol_descriptor_chain,
    find_section,
    open_gph_buffer,
    section_end,
)

TESTS_DIR = Path(__file__).resolve().parent


def read_cgns_zone_cells(cgns_path: Path) -> dict[str, int]:
    """Return {zone_name: n_cells} from a CGNS/HDF5 file."""
    zones: dict[str, int] = {}
    with h5py.File(cgns_path, "r") as f:
        if "Base" not in f:
            raise ValueError(f"{cgns_path}: missing 'Base' group")
        base = f["Base"]
        for name in base.keys():
            zone = base[name]
            if " data" not in zone:
                continue
            raw = zone[" data"][()]
            flat = raw.flatten()
            if flat.size < 2:
                raise ValueError(f"{cgns_path}: zone {name!r} has invalid ' data'")
            zones[name] = int(flat[1])
    return zones


def _cvol_id_histogram(cvol: np.ndarray) -> dict[int, int]:
    if cvol is None or len(cvol) == 0:
        return {}
    uniq, counts = np.unique(cvol, return_counts=True)
    return {int(u): int(c) for u, c in zip(uniq, counts)}


def _ls_parts_raw_chains(data: bytes) -> list[tuple[str, list[int]]]:
    """Return [(part_name, post-name [12,4,X,4] chain), ...] without resolving ids."""
    sec_start = find_section(data, "LS_Parts")
    if sec_start < 0:
        return []
    sec_end = section_end(data, sec_start)
    name_blocks = _ls_parts_name_blocks(data, sec_start, sec_end)
    out: list[tuple[str, list[int]]] = []
    for i, (name, _, after_trailer) in enumerate(name_blocks):
        scan_end = name_blocks[i + 1][1] if i + 1 < len(name_blocks) else sec_end
        chain = _scan_cvol_descriptor_chain(data, after_trailer, scan_end)
        out.append((name, chain))
    return out


def _explain_zone_source(
    zone_name: str,
    mesh: dict,
    n_cells: int,
    parts_with_cvol: list[tuple[str, int]],
) -> str:
    """Human-readable note on how a zone's cell mask is derived."""
    cvol_id = mesh.get("cvol_id")
    name_to_cvol = {name: cv for name, cv in parts_with_cvol}

    if zone_name == "FluidRegion":
        return "all cells (FluidRegion convention)"

    if cvol_id is None or len(cvol_id) != n_cells or not parts_with_cvol:
        return "fallback: all cells (missing/short cvol_id array or no parts)"

    if zone_name.startswith("@VPartRegion_"):
        rem = zone_name[len("@VPartRegion_") :].split("[", 1)[0]
        if rem in name_to_cvol:
            return f"cvol_id == {name_to_cvol[rem]}  (via @VPartRegion_ part {rem!r})"
        return f"@VPartRegion_ prefix; part {rem!r} not in LS_Parts"

    if zone_name.startswith("FPHPARTS."):
        candidate = zone_name[len("FPHPARTS.") :].rsplit(".", 1)[-1]
        if candidate in name_to_cvol:
            return f"cvol_id == {name_to_cvol[candidate]}  (FPHPARTS suffix part {candidate!r})"
        return f"FPHPARTS prefix; suffix part {candidate!r} not in LS_Parts"

    matches = sorted(
        (p for p, _ in parts_with_cvol if p and p in zone_name),
        key=len,
        reverse=True,
    )
    if matches:
        part = matches[0]
        return f"cvol_id == {name_to_cvol[part]}  (substring match part {part!r})"

    for part, cv in parts_with_cvol:
        if zone_name.endswith(f".{part}") or zone_name == part:
            return f"cvol_id == {cv}  (direct part zone {part!r})"

    return "fallback: all cells (no part/region name match)"


def _part_for_part_zone(zone_name: str, mesh: dict) -> tuple[str, int] | None:
    """If *zone_name* is a Part-derived zone, return (part_name, cvol_id)."""
    parts_with_cvol: list[tuple[str, int]] = mesh.get("parts_with_cvol", [])
    asm_info = mesh.get("assembly_info", {}) or {}
    part_paths: dict = asm_info.get("part_paths", {})
    root_empty_prefix = asm_info.get("root_empty_prefix")
    legacy_pa = mesh.get("part_assembly", {})

    def _zone_name_for_part(p: str) -> str:
        path = part_paths.get(p)
        if path is None:
            if root_empty_prefix:
                path = f"{root_empty_prefix}.{p}"
            else:
                a_name = legacy_pa.get(p)
                path = f"{a_name}.{p}" if a_name else p
        depth = path.count(".")
        if depth >= 2:
            return path
        return f"FPHPARTS.{path}"

    for part, cv in parts_with_cvol:
        if _zone_name_for_part(part) == zone_name:
            return part, cv
    return None


def print_partition_details(gph_path: Path, mesh: dict) -> None:
    """Print LS_CvolIdOfElements, LS_Parts chains, assemblies, and zone linkage."""
    link_data = mesh.get("link_data")
    n_cells = int(link_data["n_cells"]) if link_data else 0
    cvol = mesh.get("cvol_id")
    parts_with_cvol: list[tuple[str, int]] = mesh.get("parts_with_cvol", [])
    volume_regions = mesh.get("volume_regions", [])
    asm_info = mesh.get("assembly_info", {}) or {}

    print("  --- partition raw details ---")

    # LS_CvolIdOfElements
    print("  LS_CvolIdOfElements:")
    if cvol is None:
        print("    (section missing or unreadable)")
    else:
        print(f"    per-cell array length: {len(cvol)}")
        print(f"    mesh n_cells (LS_Links): {n_cells}")
        if len(cvol) == n_cells:
            print("    length check: OK")
        else:
            print("    length check: MISMATCH — zone masks may fall back to all cells")
        hist = _cvol_id_histogram(cvol)
        print(f"    unique cvol_ids ({len(hist)}): {sorted(hist)}")
        for cid in sorted(hist):
            print(f"      cvol_id={cid}: {hist[cid]} cells")

    # LS_VolumeRegions
    print("  LS_VolumeRegions (file order):")
    if volume_regions:
        for i, name in enumerate(volume_regions, start=1):
            print(f"    [{i}] {name!r}")
    else:
        print("    (none)")

    # LS_Parts raw + resolved
    print("  LS_Parts:")
    with open_gph_buffer(str(gph_path)) as data:
        raw_chains = _ls_parts_raw_chains(data)
    if not raw_chains and not parts_with_cvol:
        print("    (none)")
    else:
        resolved = dict(parts_with_cvol)
        for i, (name, chain) in enumerate(raw_chains, start=1):
            cv = resolved.get(name)
            chain_s = chain if chain else "(empty)"
            if cv is not None and cvol is not None and len(cvol) == n_cells:
                cells = int((cvol == cv).sum())
                in_hist = cv in _cvol_id_histogram(cvol)
                flag = "OK" if cells > 0 and in_hist else "WARN"
                print(
                    f"    [{i}] part={name!r}  chain={chain_s}  "
                    f"-> cvol_id={cv}  cells={cells}  [{flag}]"
                )
            else:
                print(
                    f"    [{i}] part={name!r}  chain={chain_s}  "
                    f"-> cvol_id={cv}  cells=(n/a)"
                )

    mapped_ids = [cv for _, cv in parts_with_cvol]
    actual_ids = sorted(_cvol_id_histogram(cvol)) if cvol is not None and len(cvol) == n_cells else []
    print("  LS_Parts <-> cvol_id mapping check:")
    if not parts_with_cvol:
        print("    no parts parsed")
    elif not actual_ids:
        print("    cannot validate (cvol array unavailable or wrong length)")
    else:
        unique_mapped = sorted(set(mapped_ids))
        print(f"    mapped cvol_ids:   {unique_mapped}")
        print(f"    actual cvol_ids:   {actual_ids}")
        if len(mapped_ids) != len(set(mapped_ids)):
            print("    WARN: duplicate cvol_id among parts")
        if set(unique_mapped) == set(actual_ids):
            print("    coverage: OK (bijection with actual set)")
        elif set(unique_mapped) <= set(actual_ids):
            print("    coverage: PARTIAL (mapped ids subset of actual, not full cover)")
        else:
            extra = sorted(set(unique_mapped) - set(actual_ids))
            print(f"    coverage: FAIL (mapped ids not in actual set: {extra})")

    # LS_Assemblies
    print("  LS_Assemblies:")
    if asm_info.get("has_assemblies"):
        prefix = asm_info.get("root_empty_prefix")
        print(f"    root_empty_prefix: {prefix!r}")
        part_paths = asm_info.get("part_paths", {})
        if part_paths:
            for part, path in part_paths.items():
                print(f"    part {part!r} -> path {path!r}")
        else:
            print("    (no part paths parsed)")
    else:
        print("    (none or not parsed)")

    # Zone plan linkage
    print("  Zone plan (cell selection logic):")
    plan = _build_zone_plan(mesh)
    part_zone_names = {
        z for z, _ in plan
        if _part_for_part_zone(z, mesh) is not None
    }
    for zone_name, mask in plan:
        count = int(mask.sum())
        if zone_name in part_zone_names:
            hit = _part_for_part_zone(zone_name, mesh)
            part, cv = hit if hit else ("?", "?")
            note = f"part {part!r}, cvol_id={cv}"
        else:
            note = _explain_zone_source(zone_name, mesh, n_cells, parts_with_cvol)
            # show which cvol_ids contribute when mask is a single-id subset
            if cvol is not None and len(cvol) == n_cells and count < n_cells:
                sub_hist = _cvol_id_histogram(cvol[mask])
                if sub_hist:
                    note += f"; cvol_ids in mask: {sub_hist}"
        print(f"    {zone_name}: {count} cells  <- {note}")

    print("  --- end partition details ---")


def gph_zone_cells(gph_path: Path) -> tuple[dict, dict]:
    """Parse GPH and return (mesh, {zone_name: n_cells})."""
    mesh = parse_gph_mesh(str(gph_path))
    link_data = mesh.get("link_data")
    if link_data is None:
        raise RuntimeError("LS_Links parse failed — cannot derive zone cell counts")
    plan = _build_zone_plan(mesh)
    return mesh, {name: int(mask.sum()) for name, mask in plan}


def compare_zone_maps(
    gph_map: dict[str, int],
    ref_map: dict[str, int],
) -> tuple[bool, list[str]]:
    """Compare two zone→cell-count maps. Return (ok, detail_lines)."""
    lines: list[str] = []
    ok = True

    gph_names = set(gph_map)
    ref_names = set(ref_map)

    only_gph = sorted(gph_names - ref_names)
    only_ref = sorted(ref_names - gph_names)
    if only_gph:
        ok = False
        lines.append(f"  zones only in GPH plan: {only_gph}")
    if only_ref:
        ok = False
        lines.append(f"  zones only in reference CGNS: {only_ref}")

    for name in sorted(gph_names & ref_names):
        gph_n = gph_map[name]
        ref_n = ref_map[name]
        if gph_n == ref_n:
            lines.append(f"  {name}: {gph_n} cells  OK")
        else:
            ok = False
            lines.append(f"  {name}: GPH={gph_n}, ref={ref_n}  MISMATCH")

    return ok, lines


def _resolve_gph_files(tests_dir: Path, patterns: list[str]) -> list[Path]:
    if patterns:
        out: list[Path] = []
        for pat in patterns:
            p = Path(pat)
            if not p.is_absolute():
                p = tests_dir / p if (tests_dir / p).exists() else ROOT / p
            if not p.is_file():
                raise FileNotFoundError(f"GPH file not found: {pat}")
            if p.suffix.lower() != ".gph":
                raise ValueError(f"Not a .gph file: {p}")
            out.append(p)
        return out
    return sorted(tests_dir.glob("*.gph"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare GPH volume-zone cell counts against *_orig.cgns references.",
    )
    parser.add_argument(
        "gph_files",
        nargs="*",
        metavar="GPH",
        help="Optional .gph file(s) to test (default: all tests/*.gph)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print LS_Parts descriptor chains, cvol_id histogram, assemblies, "
             "and per-zone cell-selection notes",
    )
    parser.add_argument(
        "--tests-dir",
        type=Path,
        default=TESTS_DIR,
        help=f"Directory containing test GPH/CGNS files (default: {TESTS_DIR})",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tests_dir = args.tests_dir.resolve()

    try:
        gph_files = _resolve_gph_files(tests_dir, args.gph_files)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1

    if not gph_files:
        print(f"No .gph files found under {tests_dir}")
        return 1

    print(f"Scanning {len(gph_files)} GPH file(s) in {tests_dir}")
    if args.verbose:
        print("Verbose partition details: ON")
    print()

    n_pass = n_fail = n_skip = n_error = 0

    for gph_path in gph_files:
        stem = gph_path.stem
        ref_path = tests_dir / f"{stem}_orig.cgns"
        if not ref_path.is_file() and gph_path.parent != tests_dir:
            ref_path = gph_path.with_name(f"{stem}_orig.cgns")

        print(f"=== {gph_path.name} ===")

        try:
            mesh, gph_map = gph_zone_cells(gph_path)
        except Exception as exc:
            print(f"  ERROR parsing GPH: {exc}")
            print("  CONCLUSION: ERROR\n")
            n_error += 1
            continue

        if args.verbose:
            print_partition_details(gph_path, mesh)
            print()

        print(f"  GPH zone plan ({len(gph_map)} zones):")
        for name, count in gph_map.items():
            print(f"    {name}: {count} cells")

        if not ref_path.is_file():
            print(f"  Reference: (no {ref_path.name})")
            print("  CONCLUSION: SKIP (no *_orig.cgns)\n")
            n_skip += 1
            continue

        print(f"  Reference: {ref_path.name}")
        try:
            ref_map = read_cgns_zone_cells(ref_path)
        except Exception as exc:
            print(f"  ERROR reading reference CGNS: {exc}")
            print("  CONCLUSION: ERROR\n")
            n_error += 1
            continue

        ok, details = compare_zone_maps(gph_map, ref_map)
        print("  Comparison:")
        for line in details:
            print(line)

        if ok:
            print("  CONCLUSION: PASS\n")
            n_pass += 1
        else:
            print("  CONCLUSION: FAIL\n")
            n_fail += 1

    print("=" * 60)
    print(f"Summary: {n_pass} passed, {n_fail} failed, {n_skip} skipped, {n_error} errors")
    if n_fail or n_error:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

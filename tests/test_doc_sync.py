"""Enforce the emission half of the sync invariant.

Every number that appears in more than one place has exactly one generator.
Tables in docs that restate artifact data live between ``sync:NAME`` markers and
are written by ``scripts/emit_doc_tables.py``; this test fails if any of them
drifts from its artifact, or if a marker pair and its builder get separated.

Hand-editing an emitted block is the failure mode this catches: the edit
survives until someone reruns the emitter, and until then two copies of the
same number disagree.
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[1]
MODULE_PATH = ROOT / "scripts" / "emit_doc_tables.py"

if not MODULE_PATH.exists():  # sdist may omit the emitter
    pytest.skip("doc table emitter not present", allow_module_level=True)

SPEC = importlib.util.spec_from_file_location("emit_doc_tables", MODULE_PATH)
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def available(path: Path) -> bool:
    return path.exists()


def test_every_emitted_block_is_current():
    """The whole point: run --check and require no drift."""
    if not (ROOT / "results" / "data").exists():
        pytest.skip("published dataset not present (sdist build)")
    drift, missing = mod.apply_blocks(write=False)
    assert not missing, "markers missing:\n" + "\n".join(missing)
    assert not drift, (
        "emitted blocks differ from their artifacts, run "
        "`uv run python scripts/emit_doc_tables.py --write`:\n" + "\n".join(drift)
    )


def test_check_mode_never_writes(tmp_path):
    """--check must be safe to run on a clean tree."""
    if not (ROOT / "results" / "data").exists():
        pytest.skip("published dataset not present (sdist build)")
    before = {
        path: path.read_bytes()
        for path, _build in mod.BLOCKS.values()
        if available(path)
    }
    mod.apply_blocks(write=False)
    for path, content in before.items():
        assert path.read_bytes() == content, f"--check modified {path}"


def test_every_block_has_both_markers_in_its_document():
    for name, (path, _build) in mod.BLOCKS.items():
        if not available(path):
            pytest.skip(f"{path} not present")
        text = path.read_text()
        assert mod.START.format(name=name) in text, f"{name}: start marker missing"
        assert mod.END.format(name=name) in text, f"{name}: end marker missing"


def test_no_document_carries_a_marker_without_a_builder():
    """An orphaned marker pair would silently stop being maintained."""
    import re

    known = set(mod.BLOCKS)
    pattern = re.compile(r"<!-- sync:([a-z0-9-]+):start -->")
    for path in list(ROOT.glob("*.md")) + list((ROOT / "docs").glob("*.md")):
        for found in pattern.findall(path.read_text()):
            assert found in known, (
                f"{path.name} has a sync block {found!r} with no builder in "
                f"scripts/emit_doc_tables.py"
            )


def test_the_invariant_covers_what_the_plan_named():
    """Initial coverage is the README table, the bracket tally, four benchmarks."""
    assert "census-coverage" in mod.BLOCKS
    assert "bracket-3-11" in mod.BLOCKS
    for name in (
        "compiled-sparse-benchmark",
        "rdm-method-benchmark",
        "exact-observable-benchmark",
        "generalized-observable-benchmark",
    ):
        assert name in mod.BLOCKS


def test_benchmark_blocks_carry_the_machine_and_control_caveats():
    """Benchmarks are artifacts, not prose, and the caveats ride with them.

    Multi-machine results are never averaged, and the equally-optimized dense
    control does not exist yet. Both statements are emitted into the footer so
    they cannot be dropped by editing the doc.
    """
    if not (ROOT / "results" / "data").exists():
        pytest.skip("published dataset not present (sdist build)")
    for name in (
        "compiled-sparse-benchmark",
        "rdm-method-benchmark",
        "exact-observable-benchmark",
        "generalized-observable-benchmark",
    ):
        _path, build = mod.BLOCKS[name]
        body = build()
        assert "machine" in body.lower(), f"{name}: no machine statement"
        assert "independently optimized dense baseline" in body, f"{name}: no control caveat"


def test_the_readme_coverage_table_reconciles_with_the_ledger():
    """The reconciliation the cleanup was asked for, asserted rather than eyeballed."""
    if not (ROOT / "results" / "data").exists():
        pytest.skip("published dataset not present (sdist build)")
    body = mod.census_coverage()
    assert "| **Total** | **799** | **643** | **156** | **799** |" in body
    assert "| (3, 8) | 38 |" in body, "the (3,8) block must appear"

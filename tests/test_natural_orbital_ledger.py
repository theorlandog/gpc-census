"""The natural-orbital ledger: every interference vertex gets an s_Q^NO bound.

results/data/states_natural_orbital.jsonl carries a rotated representative for
each of the 156 interference records, whose 1-RDM is exactly diag(lambda).
These tests re-verify the defining property from the shipped amplitudes rather
than trusting the generator, and pin the cost so it cannot be quietly forgotten.
"""
import json
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

LEDGER = DATA / "states_natural_orbital.jsonl"
SUMMARY = DATA / "natural_orbital_summary.json"


@pytest.fixture(scope="module")
def rows():
    return [json.loads(ln) for ln in LEDGER.read_text().splitlines() if ln.strip()]


@pytest.fixture(scope="module")
def summary():
    return json.loads(SUMMARY.read_text())


def _amps(r):
    return np.array(r["amplitudes_real"], dtype=complex) + 1j * np.array(
        r["amplitudes_imag"], dtype=float
    )


def test_covers_every_interference_record(rows):
    assert len(rows) == 156
    assert {r["classified"] for r in rows} == {"INTERFERENCE"}


def test_every_record_is_a_natural_orbital_representative(rows):
    """Recomputed here from the shipped amplitudes, not read off the summary."""
    from gpc_census.fiber import one_rdm

    for r in rows:
        lam = np.array([v / r["denominator"] for v in r["integer_form"]], dtype=float)
        dets = [tuple(t) for t in r["support_dets"]]
        rho = one_rdm(_amps(r), dets, len(lam))
        off = np.max(np.abs(rho - np.diag(np.diag(rho))))
        assert off < 1e-9, (r["system"], r["index"], off)
        assert np.max(np.abs(np.real(np.diag(rho)) - lam)) < 1e-9
        assert abs(np.linalg.norm(_amps(r)) - 1.0) < 1e-9


def test_every_record_passes_the_attainer_gate(rows):
    from gpc_census.validate import check_vertex_attainer

    for r in rows:
        lam = [v / r["denominator"] for v in r["integer_form"]]
        errs = check_vertex_attainer(
            _amps(r), [tuple(t) for t in r["support_dets"]], lam, d=len(lam)
        )
        assert errs == [], (r["system"], r["index"], errs)


def test_rotation_preserved_the_state(rows, summary):
    """The 2-RDM spectrum is basis-free, so it pins that nothing changed."""
    assert summary["worst_two_rdm_deviation"] < 1e-9
    for r in rows:
        assert r["two_rdm_spectrum_deviation"] < 1e-9


def test_the_sparsity_cost_is_recorded(summary):
    """The repair is not free; pinned so the tradeoff is not lost."""
    assert summary["support_grew"] == 142
    assert summary["support_shrank"] == 0
    assert summary["support_after"]["mean"] > summary["support_before"]["mean"]
    assert summary["support_after"]["max"] == 79
    assert summary["evidence"] == "NUMERICAL"


def test_the_exact_library_is_not_replaced(rows):
    """states.jsonl keeps its exact closed forms; this ledger is additive."""
    exact = [
        json.loads(ln)
        for ln in (DATA / "states.jsonl").read_text().splitlines()
        if ln.strip()
    ]
    assert len(exact) == 799
    interference = [r for r in exact if r["classified"] == "INTERFERENCE"]
    assert len(interference) == 156
    for r in interference:
        assert "closed_form" in r  # still symbolic
    for r in rows:
        assert r["evidence"].startswith("NUMERICAL")
        assert r["bounds"] == "s_Q^NO"

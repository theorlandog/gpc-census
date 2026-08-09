"""Guard test for the GPC-for-chemistry falsification study.

The study's conclusion is negative, so the numbers that carry it are exactly
the ones a future edit could quietly soften. These tests re-derive the
load-bearing ones from scratch rather than reading them back:

- the selection rule holds when the natural spectrum is nondegenerate and fails
  when it is not, which is the hypothesis closed-shell chemistry violates;
- no facet of any known system admits a one-sided omitted-weight bound;
- every reported energy is variational and respects the certified bound.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

ARTIFACT = ROOT / "results" / "data" / "chemistry_gpc_value_test.json"

cgvt = pytest.importorskip("chemistry_gpc_value_test")


@pytest.fixture(scope="module")
def payload():
    if not ARTIFACT.is_file():
        pytest.skip("chemistry_gpc_value_test.json not present")
    with ARTIFACT.open() as fh:
        return json.load(fh)


def test_artifact_is_shaped_as_documented(payload):
    for key in ("records", "experiment0_selection_rule",
                "experiment3_quasipinning", "selection_subspace_scaling"):
        assert key in payload, key
    recs = payload["records"]
    assert len(recs) > 100
    fams = {r["family"] for r in recs}
    assert fams == {"chemistry", "model"}
    # the chemistry side has to reach every active space the census covers
    systems = {(r["n"], r["d"]) for r in recs if r["family"] == "chemistry"}
    assert {(3, 6), (3, 8), (4, 8), (3, 10), (4, 10), (5, 10)} <= systems


def test_natural_orbital_basis_is_actually_natural(payload):
    """Every record must have diagonalized its own 1-RDM to machine precision."""
    for r in payload["records"]:
        assert r["diag_residual"] < 1e-9, r["name"]
        assert r["lambda_residual"] < 1e-9, r["name"]
        lam = np.array(r["lambda"])
        assert np.all(np.diff(lam) <= 1e-12), r["name"]
        assert abs(lam.sum() - r["n"]) < 1e-9, r["name"]


def test_no_known_facet_admits_a_one_sided_weight_bound(payload):
    """The exhaustive claim behind the certification NO-GO, recomputed.

    D(lambda) = sum_I |c_I|^2 D(I) is exact, so it bounds the weight outside the
    pinned subspace only if every D(I) is nonnegative. Not one row of the twelve
    known systems has that property, which is why no elementary bound exists.
    """
    fresh = cgvt.selection_subspace_census()
    stored = {row["system"]: row for row in payload["selection_subspace_scaling"]}
    assert len(fresh) == len(stored)
    # the population size is quoted in the write-up, so pin it
    assert sum(row["n_gpc_rows"] for row in fresh) == 786
    for row in fresh:
        ref = stored[row["system"]]
        assert row["rows_with_nonnegative_values"] == 0
        assert ref["rows_with_nonnegative_values"] == 0
        assert row["n_gpc_rows"] == ref["n_gpc_rows"]
        assert row["full_dim"] == ref["full_dim"]
        assert row["most_negative_row_value"] <= -1.0


def test_selection_rule_holds_exactly_when_the_spectrum_is_nondegenerate(payload):
    entries = payload["experiment0_selection_rule"]["nondegenerate"]
    assert len(entries) >= 8
    for e in entries:
        assert abs(e["slack"]) < 1e-10, e
        assert abs(e["w_out"]) < 1e-8, e
        assert e["min_gap"] > 1e-3, e
    # and re-derive one case from scratch, so the claim is not just read back
    sysg = cgvt.GpcSystem.load(3, 6)
    row = sysg.rows[sysg.gpc_rows()[0]]
    best = cgvt.minimize_slack(3, 6, row, seed=11, tries=3)
    assert abs(best["slack"]) < 1e-9
    assert abs(best["w_out"]) < 1e-7


def test_selection_rule_fails_on_the_degenerate_census(payload):
    """Every census attainer has a degenerate spectrum and most break the rule."""
    audit = payload["experiment0_selection_rule"]["census_audit"]
    assert audit["nondegenerate"] == 0
    assert audit["diagonal_gate"] == 643
    assert audit["rule_holds"] < audit["diagonal_gate"] // 2
    assert audit["rule_holds_after_tie_relabel"] < audit["diagonal_gate"]


def test_stage_chain_and_variationality(payload):
    for r in payload["records"]:
        a = r["stage_A_full"]["dim"]
        bdim = r["stage_B_sz"]["dim"]
        c = r["stage_C_pg"]["dim"]
        e = r["stage_E_gpc_eq"]["dim"]
        assert a >= bdim >= c >= e >= 1, r["name"]
        assert r["stage_B_sz"]["w_out"] <= r["stage_E_gpc_eq"]["w_out"] + 1e-9
        for st in ("stage_B_sz", "stage_C_pg", "stage_D_pauli", "stage_E_gpc_eq"):
            de = r[st]["de"]
            if np.isfinite(de):
                assert de >= -1e-9, (r["name"], st, de)
        for f in r["stage_F_gpc"]:
            assert f["dim"] <= c, r["name"]
            assert f["w_out"] >= r["stage_E_gpc_eq"]["w_out"] - 1e-9, r["name"]
            if np.isfinite(f["de"]):
                assert f["de"] >= -1e-9, (r["name"], f["de"])


def test_certified_bound_is_never_violated(payload):
    """E_S - E_0 <= w (E_max - E_0) / (1 - w), checked on every record."""
    checked = 0
    for r in payload["records"]:
        for f in r["stage_F_gpc"]:
            bound = f.get("certified_de_from_measured_w")
            if bound is None or not np.isfinite(bound) or not np.isfinite(f["de"]):
                continue
            span = r["e_max_sector"] - r["e0_sector"]
            assert f["de"] <= bound + 1e-8 * max(1.0, span), (r["name"], f["de"], bound)
            checked += 1
    assert checked > 100


def test_certified_bound_formula():
    assert cgvt.certified_energy_bound(0.0, -1.0, 3.0) == 0.0
    assert cgvt.certified_energy_bound(0.5, -1.0, 3.0) == pytest.approx(4.0)
    assert not np.isfinite(cgvt.certified_energy_bound(1.0, -1.0, 3.0))


def test_quasipinning_does_not_control_omitted_weight(payload):
    """The directed hunt must still contain a decisive counterexample.

    A usable certificate needs W_out <= C * D with a small C. The frontier
    records states whose ratio D / W_out is below 1e-3, so C would have to
    exceed 1000, and those states sit at a collapsing spectrum gap.
    """
    fr = payload["experiment3_quasipinning"]["frontier"]
    assert fr
    ratios = [f["ratio_D_over_W"] for f in fr
              if f.get("ratio_D_over_W") is not None and f["w_out"] > 1e-3]
    assert ratios
    assert min(ratios) < 1e-3
    worst = min(fr, key=lambda f: f["ratio_D_over_W"]
                if f.get("ratio_D_over_W") is not None else 1e9)
    assert worst["min_gap"] < 1e-4


def test_exact_slack_identity_on_a_record(payload):
    """D(lambda) = sum_I |c_I|^2 D(I) is what makes the negative result exact.

    Verified here on the polytope side: the stored spectrum must reproduce the
    stored nearest-facet slack from the stored row.
    """
    for r in payload["records"][:200]:
        if r["nearest_gpc_row"] is None:
            continue
        lam = np.array(r["lambda"])
        coeffs = np.array(r["nearest_gpc_row"], dtype=float)
        sysg = cgvt.GpcSystem.load(r["n"], r["d"])
        rhs = next(x["rhs"] for x in sysg.rows
                   if x["coeffs"] == r["nearest_gpc_row"])
        assert abs((rhs - coeffs @ lam) - r["nearest_gpc_slack"]) < 1e-9, r["name"]

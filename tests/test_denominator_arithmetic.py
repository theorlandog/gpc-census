"""The denominator-arithmetic artifact and its pre-registered scores.

results/data/denominator_arithmetic.json derives the state denominator in two
stages, incidence and then pinning. These tests pin the worked examples and
the scored predictions of docs/prereg_denominator_arithmetic.md, and re-derive
the v89 minting equation from the ledger with independent bookkeeping.
"""
import json
import pathlib
import sys

import pytest
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = ROOT / "tests" / "data"
if not DATA.exists():  # sdist ships results/data, not the symlink
    DATA = ROOT / "results" / "data"

ARTIFACT = DATA / "denominator_arithmetic.json"
LEDGER = DATA / "states.jsonl"


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


def test_v89_mints_five_through_a_linear_equation(artifact):
    v89 = artifact["v89"]
    assert v89["channel_class_size"] == 2
    assert v89["quadratic_terms_cancel"]
    assert v89["surviving_equation"] == "130*p9 - 1 = 0"
    assert v89["pinning_pivot"] == "130"
    assert v89["state_denominator"] == 130
    assert v89["minted_primes"] == [5]
    assert v89["ratio"] == 5
    assert v89["agrees_with_ledger"]


def test_v89_recomputed_from_the_ledger():
    """The incidence line and the silence condition, derived here from scratch."""
    rec = next(
        r
        for r in (json.loads(ln) for ln in LEDGER.read_text().splitlines() if ln.strip())
        if r["system"] == "(3,10)" and r["index"] == 89
    )
    dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
    n = len(dets)
    lam = [sp.Rational(v, rec["denominator"]) for v in rec["integer_form"]]
    p = sp.symbols(f"p0:{n}")

    rows = [[1 if m in T else 0 for T in dets] for m in range(len(lam))]
    eqs = [sp.Eq(sum(rows[m][t] * p[t] for t in range(n)), lam[m]) for m in range(len(lam))]
    eqs.append(sp.Eq(sum(p), 1))
    sol = sp.solve(eqs, p, dict=True)[0]
    free = [v for v in p if v not in sol]
    assert len(free) == 1  # a line
    t = free[0]

    # the one-hop classes, grouped here rather than imported: two determinants
    # are one hop apart when their symmetric difference is a single mode pair,
    # and that pair names the channel
    classes: dict[tuple, list] = {}
    for i, u in enumerate(dets):
        for j, v in enumerate(dets):
            if j <= i:
                continue
            diff = set(u) ^ set(v)
            if len(diff) == 2:
                classes.setdefault(tuple(sorted(diff)), []).append((i, j))
    assert list(classes) == [(3, 7)]  # a single channel at v89
    (u1, v1), (u2, v2) = classes[(3, 7)]

    from gpc_census.states import one_hop_classes  # cross-check the library

    assert one_hop_classes(dets) == classes

    def value(k):
        return sol.get(p[k], p[k])

    silence = sp.expand(value(u1) * value(v1) - value(u2) * value(v2))
    assert sp.Poly(silence, t).degree() == 1  # the quadratic terms cancelled
    root = sp.solve(sp.Eq(silence, 0), t)
    assert root == [sp.Rational(1, 130)]

    weights = [sp.Rational(sp.together(value(k)).subs(t, root[0])) for k in range(n)]
    assert sum(weights) == 1
    assert sp.ilcm(1, *[w.q for w in weights]) == 130
    ledger = [sp.Rational(w, rec["closed_form"]["den"]) for w in rec["closed_form"]["weights"]]
    assert weights == ledger


def test_v103_mints_thirteen_between_two_channels(artifact):
    v103 = artifact["v103"]
    assert v103["incidence_unimodular"]
    assert v103["incidence_solution_denominator"] == 34  # incidence mints nothing
    assert v103["incidence_kernel_dim"] == 5
    assert v103["all_channels_silent_at_certified_state"]
    assert v103["minted_primes"] == [13]
    assert v103["ratio"] == 5746

    cancellations = v103["pairwise_quadratic_cancellations"]
    assert cancellations, "expected the v89 cancellation to recur between channels"
    thirteen = [
        c
        for c in cancellations
        if c["coefficient_primes_beyond_the_spectrum_denominator"] == [13]
    ]
    assert thirteen
    assert thirteen[0]["surviving_linear_equation"] == "-13*p12 + p13 + p9 = 0"


def test_prereg_scores_are_as_recorded(artifact):
    """The pre-registered predictions, scored once and not re-editable."""
    census = artifact["census"]
    assert census["states_scored"] == 798
    assert census["excluded"] == ["(4,9) v65"]
    assert census["P_DEN_2"]["verdict"] == "FAIL"
    assert census["P_DEN_2"]["exceptions"] == 18
    assert census["P_DEN_3"]["verdict"] == "FAIL"
    assert census["P_DEN_3"]["exceptions"] == 1
    assert census["P_DEN_3"]["exception_list"] == ["(3,10) v103"]
    assert census["P_DEN_3"]["all_exceptions_are_pinned"]
    assert census["P_DEN_4"]["verdict"] == "PASS"
    assert census["P_DEN_4"]["exceptions"] == 0
    assert census["P_DEN_4"]["population"] == 736


def test_ratio_not_one_splits_by_stage(artifact):
    """The 14 states above the spectrum denominator, by which stage minted."""
    split = artifact["census"]["ratio_not_one_split"]
    assert split["status"].startswith("POST-HOC")
    assert split["count"] == 14
    assert len(split["incidence_stage"]["states"]) == 12
    assert split["incidence_stage"]["all_design_real"]
    assert split["incidence_stage"]["ratios"] == ["2"]
    assert split["pinning_stage"]["states"] == ["(3,10) v89", "(3,10) v103"]
    assert split["pinning_stage"]["all_interference"]

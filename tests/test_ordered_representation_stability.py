"""Guard the ordered representation-stability result.

Two things could rot here and they fail differently, so they are tested
differently.

The MATHEMATICS is re-derived from the published constraint tables inside this
module, without reading the artifact, so a bug in the measuring script cannot
make the tests agree with it. That covers the two theorems the write-up proves
(the ordered category is a poset, and chamber-facing rows are not permutation
stable) and the inheritance law that the constructor is invited to rely on.

The ARTIFACT is then checked for agreement with those independent derivations
and for internal consistency, and the write-up is checked to still state the
numbers it cites. A doc that drifts from its artifact is exactly the failure
mode `docs/AGENT_PLAN.md` calls the sync invariant.
"""
from __future__ import annotations

import json
import pathlib
import sys
from fractions import Fraction

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from ordered_representation_stability import (  # noqa: E402
    build_rows,
    dot,
    is_chamber_row,
    is_valid,
    load_systems,
    phi_core,
    pi_pad,
    psi_core,
    rho,
    width,
)

DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "ordered_representation_stability.json"
NOTE = ROOT / "docs" / "ordered_representation_stability.md"
PREREG = ROOT / "docs" / "prereg_ordered_representation_stability.md"


@pytest.fixture(scope="module")
def systems():
    return load_systems()


@pytest.fixture(scope="module")
def taus(systems):
    return build_rows(systems)


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


# --------------------------------------------------------------------------
# the mathematics, re-derived without the artifact
# --------------------------------------------------------------------------


def test_terminal_restriction_is_total(systems, taus):
    """Theorem: the last-mode face of the ordered slice IS the lower slice."""

    checked = 0
    for (n, d), rows in systems.items():
        if (n, d - 1) not in systems:
            continue
        lower = systems[n, d - 1]["vertices"]
        for tau in taus[n, d]:
            image = rho(tau)
            if image is None:
                continue
            checked += 1
            assert is_valid(image, lower), f"({n},{d}) row {tau} restricted invalid"
    assert checked == 724


def test_nonterminal_restriction_is_not_total(systems, taus):
    """Theorem 1's computational shadow: middle insertions do not act."""

    found = False
    for (n, d), rows in systems.items():
        if (n, d - 1) not in systems:
            continue
        lower = systems[n, d - 1]["vertices"]
        for tau in taus[n, d]:
            for j in range(d - 1):
                head = tau[:j] + tau[j + 1:]
                if len(set(head)) <= 1:
                    continue
                if not is_valid(head, lower):
                    found = True
                    break
            if found:
                break
        if found:
            break
    assert found, "some non-terminal restriction must fail, or Theorem 1 is wrong"


def test_weyl_saturation_is_not_convex(systems, taus):
    """Theorem 2: (1^(N-1), 1/2, 1/2, 0^(d-N-1)) is a non-attainable midpoint.

    Both endpoints are permutations of the Slater spectrum, so both lie in the
    saturation; the midpoint is already dominant, so membership is decided by
    the published system itself.
    """

    half = Fraction(1, 2)
    for (n, d), payload in systems.items():
        if not payload["primary"]:
            continue
        midpoint = (Fraction(1),) * (n - 1) + (half, half) + (Fraction(0),) * (d - n - 1)
        assert len(midpoint) == d
        assert sorted(midpoint, reverse=True) == list(midpoint)
        assert any(dot(tau, midpoint) > 0 for tau in taus[n, d]), (
            f"({n},{d}) midpoint unexpectedly attainable, Theorem 2 would fail")


def test_no_nonstructural_row_is_permutation_stable(systems, taus):
    """Theorem 2's consequence: S_d does not act on chamber-facing rows.

    By the rearrangement inequality the worst permutation pairs the sorted row
    with the sorted vertex, so no permutation has to be enumerated.
    """

    survivors = []
    for (n, d), payload in systems.items():
        sorted_vertices = [sorted(v, reverse=True) for v in payload["vertices"]]
        for index, tau in enumerate(taus[n, d]):
            worst = sorted(tau, reverse=True)
            if all(sum((Fraction(a) * b for a, b in zip(worst, v)), Fraction(0)) <= 0
                   for v in sorted_vertices):
                survivors.append((n, d, index, is_chamber_row(tau, n)))
    assert survivors, "expected the structural Pauli rows to survive"
    assert all(structural for *_, structural in survivors), (
        f"a nonstructural row survived permutation: {survivors}")
    assert len(survivors) == 6


def test_tight_extensions_are_sections_and_preserve_facets(systems, taus):
    """The inheritance law, re-derived: 612 of 612 with the two families named."""

    degenerate = {key for key, payload in systems.items() if payload["equalities"]}
    eligible = preserved = sections = 0
    for (n, d), rows in sorted(systems.items()):
        for target, apply_map in (
            ((n, d + 1), lambda tau, t: pi_pad(tau, systems[t]["vertices"])),
            ((n + 1, d + 1), lambda tau, t: psi_core(tau, n, systems[t]["vertices"])),
        ):
            if target not in systems:
                continue
            published = set(taus[target])
            back = rho if target == (n, d + 1) else (lambda tau: phi_core(tau, n + 1))
            for tau in taus[n, d]:
                image = apply_map(tau, target)
                assert image is not None, f"({n},{d}) tight extension undefined"
                assert is_valid(image, systems[target]["vertices"])
                sections += back(image) == tau
                if is_chamber_row(tau, n) or (n, d) in degenerate:
                    continue
                eligible += 1
                preserved += image in published
    assert sections == 620
    assert (eligible, preserved) == (612, 612)


# --------------------------------------------------------------------------
# the artifact and the write-up
# --------------------------------------------------------------------------


def test_theorem_1_hypothesis_holds(systems, taus):
    """Non-degeneracy: the slice is not trapped in {lambda_d = 0}.

    Theorem 1 needs an interior point with every coordinate positive, and the
    Weyl-averaging shortcut is unavailable because Theorem 2 rules it out, so
    the point is exhibited: the barycenter of the vertex list, which is in the
    slice by convexity.
    """

    for (n, d), payload in systems.items():
        vertices = payload["vertices"]
        centre = tuple(sum(v[i] for v in vertices) / len(vertices) for i in range(d))
        assert all(x > 0 for x in centre), f"({n},{d}) barycentre touches a coordinate wall"
        assert all(dot(tau, centre) <= 0 for tau in taus[n, d]), f"({n},{d}) barycentre outside"


def test_every_published_row_is_a_facet(systems, taus):
    """The disproof counts facets, so a redundant row would invalidate it.

    Re-derived here rather than read from the artifact: a row is a facet
    exactly when its tight vertices span one dimension below the slice.
    """

    from ordered_representation_stability import _affine_rank

    total = 0
    for (n, d), payload in systems.items():
        vertices = payload["vertices"]
        ambient = _affine_rank(vertices)
        assert ambient == (3 if (n, d) == (3, 6) else d - 1), f"({n},{d}) slice dimension"
        for tau in taus[n, d]:
            tight = [v for v in vertices if dot(tau, v) == 0]
            assert _affine_rank(tight) == ambient - 1, f"({n},{d}) row {tau} is not a facet"
            total += 1
    assert total == 907


def test_artifact_facetness_agrees(artifact):
    report = artifact["facetness"]
    assert all(entry["all_facets"] for entry in report)
    assert sum(entry["facets"] for entry in report) == 907
    assert sum(entry["redundant_or_lower_face"] for entry in report) == 0


def test_artifact_nondegeneracy_agrees(artifact):
    report = artifact["nondegeneracy"]
    assert len(report) == 15
    assert all(entry["all_coordinates_positive"] and entry["in_slice"] for entry in report)


def test_artifact_scorecard_is_pinned(artifact):
    """The scored verdicts are the record; three FAILs are findings, not bugs."""

    expected = {
        "P1_terminal_restriction_total": "CONFIRMATORY-PASS",
        "P2_nonterminal_restriction_fails": "PASS",
        "P3_tight_padding_is_a_section": "PASS",
        "P4_extension_does_not_preserve_facetness": "PASS",
        "P5_primitive_counts_do_not_vanish": "FAIL",
        "P6_width_unbounded": "PASS",
        "P7_height_unbounded": "CONFIRMATORY-FAIL",
        "P8a_holdout_3_10": "PASS",
        "P8b_holdout_4_10": "PASS",
        "P8c_holdout_5_10": "FAIL",
        "P9_composability": "PASS",
    }
    score = artifact["scorecard"]
    assert {k: score[k]["verdict"] for k in expected} == expected
    # Standing rule R2: every verdict reports scope and independent scope, and
    # independent scope is strictly smaller wherever inherited rows were scored.
    for key in expected:
        assert score[key]["independent_scope"] <= score[key]["scope"], key
        assert score[key]["independent_scope"] > 0, key


def test_prereg_scorecard_matches_the_artifact(artifact):
    """The scored table in the pre-registration is quoted, not retyped."""

    text = PREREG.read_text()
    assert "## SCORED" in text
    for key, label in (
        ("P5_primitive_counts_do_not_vanish", "P5 primitive counts do not vanish"),
        ("P8c_holdout_5_10", "P8c holdout (5,10)"),
        ("P9_composability", "P9 composability"),
    ):
        entry = artifact["scorecard"][key]
        row = next(line for line in text.splitlines() if line.startswith(f"| {label} |"))
        cells = [cell.strip() for cell in row.strip("|").split("|")]
        assert cells[1].strip("*") == entry["verdict"], label
        assert cells[2] == str(entry["scope"]), label
        assert cells[3] == str(entry["independent_scope"]), label


def test_artifact_core_replay_is_complete(artifact):
    replay = artifact["core_replay"]
    assert replay["all_reproduced"] is True
    assert replay["inherited_rows_checked"] == replay["rebuilt_exactly"] == 514
    assert replay["failures"] == []


def test_artifact_label_counts_close(artifact):
    """Every row carries exactly one label and the labels sum to the row count."""

    for system, summary in artifact["summary"].items():
        assert sum(summary["labels"].values()) == summary["rows"], system
        assert summary["independent_scope"] == summary["primitive"], system


def test_artifact_primitive_counts_do_not_stabilise(artifact):
    summary = artifact["summary"]
    assert [summary[f"(3,{d})"]["primitive"] for d in range(6, 11)] == [1, 4, 27, 21, 41]
    assert [summary[f"(3,{d})"]["rows"] for d in range(6, 11)] == [1, 4, 31, 52, 93]
    assert [summary[f"(4,{d})"]["rows"] for d in (8, 9, 10)] == [15, 60, 125]
    assert summary["(5,10)"]["rows"] == 161


def test_artifact_holdout(artifact):
    held = artifact["holdout_rank10"]
    assert (held["(3,10)"]["covered"], held["(3,10)"]["published"]) == (52, 93)
    assert (held["(4,10)"]["covered"], held["(4,10)"]["published"]) == (80, 125)
    assert (held["(5,10)"]["covered"], held["(5,10)"]["published"]) == (104, 161)
    assert all(v["coverage_percent"] < 100.0 for v in held.values())


def test_artifact_records_row_identities_not_only_counts(artifact):
    """The brief requires exact row identities and transport witnesses."""

    rows = artifact["rows"]["(3,10)"]
    assert len(rows) == 93
    assert all(len(record["tau"]) == 10 for record in rows)
    inherited = [r for r in rows if r["label"] == "PADDING"]
    assert inherited and all(r["donor"] and r["witness"] for r in inherited)
    assert all(len(entry["word"]) >= 0 for entry in artifact["cores"]["(3,10)"])


def test_write_up_quotes_the_artifact(artifact):
    """Load-bearing numbers in the note must still be the measured ones."""

    text = NOTE.read_text()
    preservation = artifact["facet_preservation"]
    replay = artifact["core_replay"]
    obstruction = artifact["symmetry_obstruction"]["global_validity"]
    nonstructural = (sum(r["rows"] for r in obstruction)
                     - sum(r["of_those_structural"] for r in obstruction))
    for needle in (
        f"{preservation['landed_on_published_facet']} of "
        f"{preservation['nonstructural_eligible']}",
        f"{replay['rebuilt_exactly']} inherited rows",
        f"0 of {nonstructural} nonstructural",
        "QUESTION ANSWER:",
        "MISSING LEMMA:",
    ):
        assert needle in text, f"write-up no longer states: {needle}"
    assert "docs/prereg_ordered_representation_stability.md" in text
    assert PREREG.exists()


def test_width_is_gauge_free_and_not_an_obstruction(artifact):
    """P6 passed, but width does not separate inherited rows from primitive ones."""

    assert width((3, 1, 1, 1)) == 1
    assert width((5, 5, 5, 5)) == 0
    # lambda_4 >= 0 reads as lambda_1 + lambda_2 + lambda_3 <= N on the trace
    # slice, so it is supported on the first three coordinates.
    assert width((0, 0, 0, -1)) == 3
    for system, summary in artifact["summary"].items():
        if summary["max_width_inherited"] is None:
            continue
        assert summary["max_width_primitive"] == summary["max_width_inherited"], system

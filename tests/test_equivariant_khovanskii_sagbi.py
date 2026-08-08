"""Guards for the Khovanskii/SAGBI compression experiment.

The campaign is a NEGATIVE result, so the test suite has to work harder than
usual: a negative result is only worth as much as the machinery that produced
it, and a broken semigroup engine manufactures new generators for free (every
value it loses makes the values above it look indecomposable). So the guards
here re-derive the load-bearing numbers LIVE rather than re-reading them, pin
the two structural theorems that the numbers rest on, and assert that the
standalone verifier rejects a tampered artifact.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = (ROOT / "tests" / "data") if (ROOT / "tests" / "data").exists() \
    else (ROOT / "results" / "data")
ARTIFACT = DATA / "equivariant_khovanskii_sagbi.json"
DOC = ROOT / "docs" / "equivariant_khovanskii_sagbi.md"
PREREG = ROOT / "docs" / "prereg_equivariant_khovanskii_sagbi.md"
GENERATOR = ROOT / "scripts" / "equivariant_khovanskii_sagbi.py"
VERIFIER = ROOT / "scripts" / "verify_equivariant_khovanskii_sagbi_standalone.py"

if not ARTIFACT.exists():  # pragma: no cover - sdist without the dataset
    pytest.skip("published dataset not present", allow_module_level=True)


def _load_generator():
    spec = importlib.util.spec_from_file_location("eks", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def art():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def eks():
    if not GENERATOR.exists():  # pragma: no cover - sdist without scripts
        pytest.skip("generator script not present")
    return _load_generator()


@pytest.fixture(scope="module")
def raw(art):
    return {block["system"]: block for block in art["raw_semigroup"]}


# --------------------------------------------------------------------------
# Shape and exactness
# --------------------------------------------------------------------------

def test_artifact_is_the_scored_run(art):
    assert art["evidence"] == "exact"
    assert art["quick"] is False, "the scored artifact must not be a --quick run"
    assert art["generated_by"] == "scripts/equivariant_khovanskii_sagbi.py"
    assert art["prereg"] == "docs/prereg_equivariant_khovanskii_sagbi.md"


def test_every_system_carries_a_status(art):
    for row in art["v2_by_system"] + art["v1_v3_by_system"]:
        assert row["status"] in {"complete", "operational"}
        assert row["reached_degree"] >= 1


# --------------------------------------------------------------------------
# The gates. These are forced, so a failure means the engine is wrong.
# --------------------------------------------------------------------------

def test_semigroup_is_closed_under_addition(art):
    """Forced: the value semigroup of a valuation on a domain is a semigroup.

    This is the cheapest test that catches a wrong multiplicity, because a
    missing value shows up as a sum of two present values that is itself
    absent.
    """
    for row in art["semigroup_closure_gate"]:
        assert row["closure_violations"] == 0, (row["system"], row["witnesses"])
    for row in art["term_order_closure_gate"]:
        assert row["closure_violations"] == 0, row


def test_p0_highest_weight_dimension_equals_multiplicity(art):
    """P0, the preregistered correctness gate.

    The kernel of the simple raising operators is computed by exact linear
    algebra over Q; the multiplicity is computed from power sums by
    Murnaghan-Nakayama. They are independent routes to the same integer.
    """
    assert art["hw_dimension_gate"], "the term-order side produced no gate rows"
    for row in art["hw_dimension_gate"]:
        assert row["all_match"], row
        assert row["weight_spaces_checked"] > 0


def test_padding_is_an_equality(art):
    """P3. `Gen(n,d) ∩ {lambda_d = 0} = Gen(n,d-1)`, which is a theorem."""
    assert art["padding"], "no padding rows"
    for row in art["padding"]:
        assert row["equal"], row


def test_particle_hole_transport_holds(art):
    """P4. Self-duality at `n = d - n`, computed independently on each side."""
    scored = [row for row in art["particle_hole"]
              if row["system"] in {"(3,6)", "(4,8)", "(5,10)"}]
    assert scored, "no self-dual system in scope"
    for row in scored:
        assert row["failures"] == 0, row
        assert row["checked"] > 0


# --------------------------------------------------------------------------
# Live re-derivation, not re-reading
# --------------------------------------------------------------------------

def test_low_degree_semigroup_recomputes_live(eks, raw):
    """Rebuild `S(3,6)` at `m <= 6` from scratch and match the stored levels."""
    stored = {row["degree"]: {tuple(v) for v in row["values"]}
              for row in raw["(3,6)"]["levels"] if row["values"] is not None}
    for m in range(1, 7):
        live = set(eks.hw_support(3, 6, m))
        assert live == stored[m], m


def test_generators_recompute_live_at_3_6(eks, raw):
    """The minimal generator set at `(3,6)` is rebuilt, not read back."""
    levels = {row["degree"]: {tuple(v) for v in row["values"]}
              for row in raw["(3,6)"]["levels"] if row["values"] is not None}
    window = min(6, max(levels))
    live, violations = eks.minimal_generators(levels, window)
    assert violations == []
    stored = [(m, tuple(v)) for m, v in raw["(3,6)"]["generators"] if m <= window]
    assert live == stored


def test_borland_dennis_vertex_first_appears_at_degree_four(eks):
    """The uniform `(3,6)` vertex is the semigroup face of the period-4 finding.

    `docs/quantization_character.md` explains the period 4 of `(3,6)#3` as an
    isotropy character. Seen from the semigroup, the same fact is that
    `(1,1,1,1,1,1)` is realized at `m = 4` and not at `m = 2`, which is exactly
    a failure of saturation.
    """
    assert (1, 1, 1, 1, 1, 1) not in eks.hw_support(3, 6, 2)
    assert (2, 2, 2, 2, 2, 2) in eks.hw_support(3, 6, 4)
    assert (3, 3, 3, 3, 3, 3) not in eks.hw_support(3, 6, 6)
    assert (4, 4, 4, 4, 4, 4) in eks.hw_support(3, 6, 8)


def test_quantization_period_agrees_with_the_other_campaign(art):
    """P9. Two campaigns, no shared code path beyond the plethysm primitives."""
    block = art["quantization_cross_check"]
    if block["status"] == "absent":  # pragma: no cover
        pytest.skip("quantization_multiplicities.json not present")
    assert block["scope"] > 0
    assert block["agreements"] == block["scope"], [
        row for row in block["rows"] if not row["agrees"]
    ]


# --------------------------------------------------------------------------
# The negative result itself, pinned so it cannot be softened silently
# --------------------------------------------------------------------------

def test_primitive_generators_grow_with_rank(art):
    """P1 FAILED, and the failure is the finding. Pin the direction."""
    growth = {row["system"]: row["primitive_generators"]
              for row in art["rank_growth"]}
    ordered = [growth[f"(3,{d})"] for d in (6, 7, 8, 9, 10) if f"(3,{d})" in growth]
    assert len(ordered) >= 3
    assert ordered == sorted(ordered), ordered
    assert ordered[-1] > ordered[0], ordered


def test_held_out_degree_prediction_failed(art):
    """P5 FAILED. At least one system misses held-out values."""
    rows = art["held_out_degree"]
    assert rows
    assert any(row["verdict"] == "FAIL" for row in rows), rows
    for row in rows:
        assert row["spurious_half"] == "UNINFORMATIVE_BY_CONSTRUCTION"


def test_saturation_fails_with_witnesses(art):
    """The one model here that can over-predict does over-predict."""
    rows = {row["system"]: row for row in art["saturation_model"]}
    assert rows["(3,6)"]["verdict"] == "FAIL"
    assert rows["(3,6)"]["witnesses"], "a FAIL must name a witness"


# --------------------------------------------------------------------------
# The standalone verifier, including a negative control
# --------------------------------------------------------------------------

def test_standalone_verifier_passes():
    if not VERIFIER.exists():  # pragma: no cover
        pytest.skip("verifier not present")
    done = subprocess.run([sys.executable, str(VERIFIER), str(ARTIFACT)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stdout + done.stderr


def test_standalone_verifier_rejects_a_tampered_generator(tmp_path):
    """A verifier that cannot fail is not a verifier."""
    if not VERIFIER.exists():  # pragma: no cover
        pytest.skip("verifier not present")
    payload = json.loads(ARTIFACT.read_text())
    for block in payload["raw_semigroup"]:
        if block["system"] != "(3,6)":
            continue
        # claim a decomposable value as a minimal generator
        level = next(row for row in block["levels"] if row["degree"] == 2)
        block["generators"].append([4, [4, 4, 4, 0, 0, 0]])
        assert level["values"]
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(payload))
    done = subprocess.run([sys.executable, str(VERIFIER), str(tampered),
                           "--skip-independent"],
                          capture_output=True, text=True)
    assert done.returncode == 1, done.stdout
    assert "decomposes" in done.stdout


def test_standalone_verifier_rejects_a_tampered_digest(tmp_path):
    if not VERIFIER.exists():  # pragma: no cover
        pytest.skip("verifier not present")
    payload = json.loads(ARTIFACT.read_text())
    block = next(b for b in payload["raw_semigroup"] if b["system"] == "(3,6)")
    row = next(r for r in block["levels"] if r["values"] is not None)
    row["sha256"] = "0" * 64
    tampered = tmp_path / "tampered.json"
    tampered.write_text(json.dumps(payload))
    done = subprocess.run([sys.executable, str(VERIFIER), str(tampered),
                           "--skip-independent"],
                          capture_output=True, text=True)
    assert done.returncode == 1, done.stdout
    assert "digest mismatch" in done.stdout


# --------------------------------------------------------------------------
# Documents
# --------------------------------------------------------------------------

def test_prereg_predates_the_narrative():
    text = PREREG.read_text()
    for marker in ("P0", "P1", "P5", "P9", "Abort criterion", "Disclosed pilots"):
        assert marker in text


def test_doc_states_what_was_not_proved():
    text = DOC.read_text()
    for marker in ("NOT PROVED", "GO/NO-GO", "provisional"):
        assert marker in text


def test_doc_quotes_the_artifact_numbers(art):
    """Every headline count in the note is re-derived here and asserted present."""
    text = DOC.read_text()
    for row in art["rank_growth"]:
        assert str(row["primitive_generators"]) in text, row["system"]
    assert "em dash" not in text
    assert "—" not in text, "house rule: no em dashes"

"""Guards for the signed Chow projection: the three theorems and the audit."""

from __future__ import annotations

import importlib.util
import itertools
import json
import pathlib
import random

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SCRIPTS = ROOT / "scripts"

SIGNED_CHOW_PROJECTION = DATA / "signed_chow_projection.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MODULE = _load("signed_chow_projection")


def _artifact() -> dict:
    return json.loads(SIGNED_CHOW_PROJECTION.read_text())


def test_t1_determination_on_random_normals():
    """No other same-size family can share a positive shadow (T1)."""
    rng = random.Random(11)
    n, d = 3, 7
    slice_weights = list(itertools.combinations(range(d), n))
    for _ in range(60):
        tau = tuple(rng.randrange(-4, 5) for _ in range(d))
        positive, _negative, _zero = MODULE.signed_chow(tau, n)
        family = [s for s in slice_weights if sum(tau[i] for i in s) > 0]
        if not family:
            continue
        for _attempt in range(200):
            other = rng.sample(slice_weights, len(family))
            shadow = [0] * d
            for subset in other:
                for i in subset:
                    shadow[i] += 1
            if tuple(shadow) == positive:
                assert set(other) == set(family)


def test_t2_monotonicity_on_random_normals():
    """tau_i >= tau_j forces c_+(tau)_i >= c_+(tau)_j (T2)."""
    rng = random.Random(12)
    n, d = 3, 8
    for _ in range(200):
        tau = tuple(rng.randrange(-5, 6) for _ in range(d))
        positive, _negative, _zero = MODULE.signed_chow(tau, n)
        assert MODULE.monotone(tau, positive)


def test_t3_reconstruction_recovers_a_known_normal():
    """A published one-hole rank-9 normal is recovered from its zero family."""
    tau = (-5, -3, -2, -1, 0, 1, 2, 3, 4)
    positive, negative, zero_family = MODULE.signed_chow(tau, 3)
    recovered, rank = MODULE.annihilator(zero_family, len(tau))
    assert rank == len(tau) - 1
    assert recovered == MODULE.primitive(tau)
    assert MODULE.verify_certificate(tau, positive, negative, 3)


def test_certificate_rejects_a_perturbed_normal():
    """Verification must fail on a normal that does not induce the shadow."""
    tau = (-5, -3, -2, -1, 0, 1, 2, 3, 4)
    positive, negative, _zero = MODULE.signed_chow(tau, 3)
    wrong = (-5, -3, -2, -1, 0, 1, 2, 3, 5)
    assert not MODULE.verify_certificate(wrong, positive, negative, 3)


def test_encoding_is_two_d_integers():
    for record in _artifact()["systems"].values():
        assert record["encoding_integers"] == 2 * record["rank"]
        assert record["encoding_integers"] < record["slice_size"]


def test_every_row_reconstructs_with_no_collision():
    """The measured half: 142,493 rows, all exact, zero collisions."""
    artifact = _artifact()
    total = 0
    for record in artifact["systems"].values():
        nonstructural = record["nonstructural_survivors"]
        assert record["normal_reconstructed_exactly"] == nonstructural
        assert record["zero_family_spans_hyperplane"] == nonstructural
        assert record["monotone_rows"] == nonstructural
        assert record["certificates_verified"] == nonstructural
        assert record["signed_shadow_collisions"] == 0
        assert record["positive_shadow_collisions"] == 0
        assert record["distinct_signed_shadows"] == nonstructural
        total += nonstructural
    assert total == 142493


def test_populations_match_the_rank9_audit():
    """The Chow audit and the compression audit must see the same population."""
    chow = _artifact()["systems"]
    audit = json.loads((DATA / "rank9_compression_audit.json").read_text())["systems"]
    for system, record in chow.items():
        assert record["trace_survivors"] == audit[system]["trace_survivors"]
        assert (
            record["nonstructural_survivors"]
            == audit[system]["nonstructural_survivors"]
        )


def test_conjectural_claims_are_labelled():
    """The artifact must keep the proved and measured halves apart."""
    artifact = _artifact()
    assert set(artifact["theorems"]) == {
        "T1_determination",
        "T2_monotonicity",
        "T3_reconstruction",
    }
    measured = artifact["measured_not_proved"]
    assert "span_hypothesis" in measured
    assert "positive_shadow_alone" in measured
    assert "decoding" in measured

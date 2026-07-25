"""Guardrail for the fixed-rho / fixed-spectrum distinction.

These are two different fibers over a vertex spectrum and conflating them has
caused repeated errors in this project. The tests below pin the one census
instance where the two verdicts DISAGREE (v103: fixed-rho rigid, fixed-spectrum
tangent 6) and one where they AGREE (v89: both rigid), and they fail if anyone
reintroduces an unqualified "fiber_dim" name anywhere in the package.
"""
import json
import pathlib
import re
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data" / "states.jsonl"


def _record(system, index):
    if not DATA.exists():
        pytest.skip("no states atlas present")
    for line in DATA.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r["system"] == system and r["index"] == index:
                return r
    raise AssertionError((system, index))


def test_v103_two_fibers_disagree():
    """v103 = (18^4,5^6)/34: fixed-rho rigid, fixed-spectrum 6-dimensional."""
    from gpc_census.fiber import (amplitudes_from_record, fixed_rho_tangent,
                                  fixed_spectrum_tangent)

    c, dets = amplitudes_from_record(_record("(3,10)", 103))
    rho_t = fixed_rho_tangent(c, dets, 10)
    spec_t = fixed_spectrum_tangent(c, dets, 10)
    assert rho_t.fiber == "fixed_rho" and spec_t.fiber == "fixed_spectrum"
    assert rho_t.tangent_dim == 0, "v103 is fixed-rho rigid (P4 scorecard sense)"
    assert spec_t.tangent_dim == 6, "v103 fixed-spectrum fiber is a 6-manifold germ"
    assert rho_t.gauge_dim == spec_t.gauge_dim == 9


def test_v89_two_fibers_agree():
    """v89 = (15^2,6^8)/26: rigid in BOTH senses (silence-rank lemma corollary)."""
    from gpc_census.fiber import (amplitudes_from_record, fixed_rho_tangent,
                                  fixed_spectrum_tangent)

    c, dets = amplitudes_from_record(_record("(3,10)", 89))
    assert fixed_rho_tangent(c, dets, 10).tangent_dim == 0
    assert fixed_spectrum_tangent(c, dets, 10).tangent_dim == 0


def test_vb_fixed_spectrum_tangent_is_a_surface():
    """v_B = (4,9) v65: the archive's reduced spectrum-fiber SURFACE, dim 2."""
    from gpc_census.fiber import amplitudes_from_record, fixed_spectrum_tangent

    c, dets = amplitudes_from_record(_record("(4,9)", 65))
    assert fixed_spectrum_tangent(c, dets, 9).tangent_dim == 2


def test_json_keys_name_their_fiber():
    from gpc_census.fiber import amplitudes_from_record, tangent_report

    c, dets = amplitudes_from_record(_record("(3,10)", 103))
    rep = tangent_report(c, dets, 10)
    assert rep["fixed_rho_tangent_dim"] == 0
    assert rep["fixed_spectrum_tangent_dim"] == 6
    for key in rep:
        if key.endswith("tangent_dim") or key.endswith("gauge_dim"):
            assert key.startswith(("fixed_rho_", "fixed_spectrum_")), key


def test_no_unqualified_fiber_dim_in_package():
    """Rule 3 of the fiber program, enforced: no generic fiber dimension name."""
    banned = re.compile(r"\bfiber_dim\b|\bfiber_kernel_dim\b|\bfiber_tangent_dim\b")
    offenders = []
    for path in sorted((ROOT / "src" / "gpc_census").rglob("*.py")):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if banned.search(line) and "banned" not in line:
                offenders.append(f"{path.relative_to(ROOT)}:{i}: {line.strip()}")
    assert not offenders, (
        "unqualified fiber dimension names found; say fixed_rho or "
        "fixed_spectrum:\n" + "\n".join(offenders))

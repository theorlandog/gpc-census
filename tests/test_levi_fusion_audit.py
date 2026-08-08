"""Guards for the exact Levi-fusion allocation arithmetic."""

from __future__ import annotations

import importlib.util
import pathlib

SCRIPT = (
    pathlib.Path(__file__).resolve().parents[1]
    / "scripts"
    / "levi_fusion_audit.py"
)
SPEC = importlib.util.spec_from_file_location("levi_fusion_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_single_sector_primitive_3_7():
    tau = (-2, 1, 1, 1, 1, -2, -2)
    result = MODULE.descriptor_from_tau(tau, 3)
    assert result["fusion_rank"] == 1
    assert result["zero_grade_dimension"] == 18
    assert result["distinct_levels"] == 2


def test_multisector_primitive_3_8():
    tau = (2, 2, -4, -1, -1, -1, -1, -1)
    result = MODULE.descriptor_from_tau(tau, 3)
    assert result["fusion_rank"] == 2
    assert result["zero_grade_dimension"] == 21
    assert sorted(
        sector["dimension"] for sector in result["active_sectors"]
    ) == [1, 20]


def test_highest_pilot_fusion_rank():
    tau = (-3, -2, -1, 0, 1, 2, -1, 1)
    result = MODULE.descriptor_from_tau(tau, 3)
    assert result["fusion_rank"] == 5
    assert result["zero_grade_dimension"] == 9
    assert result["distinct_levels"] == 6


def test_4_8_nonstructural_has_two_sectors():
    tau = (0, 0, 0, 0, 1, -1, -1, -1)
    result = MODULE.descriptor_from_tau(tau, 4)
    assert result["fusion_rank"] == 2
    assert result["zero_grade_dimension"] == 19

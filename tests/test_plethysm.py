"""Guard the plethysm inner-hull generator: it must exactly reproduce the
(3,6) vertex set at M<=6 and respect the N=2 even-pairing theorem. The full
(3,7) M<=10 gate lives in the script's own --validate (too slow for CI)."""
import pathlib
import sys
from collections import Counter
from fractions import Fraction as F

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))
import plethysm_inner_hull as ph  # noqa: E402


def test_reproduces_3_6_vertices_exactly():
    pts = ph.inner_spectra(3, 6, 6)
    dim, extreme, unverified = ph.extreme_points(pts)
    known = {
        (F(1), F(1), F(1), F(0), F(0), F(0)),
        (F(1), F(1, 2), F(1, 2), F(1, 2), F(1, 2), F(0)),
        (F(3, 4), F(3, 4), F(1, 2), F(1, 2), F(1, 4), F(1, 4)),
        (F(1, 2),) * 6,
    }
    assert dim == 3
    assert set(extreme) == known
    assert not unverified


def test_n2_components_are_evenly_paired():
    for spec in ph.inner_spectra(2, 6, 4):
        for _val, mult in Counter(x for x in spec if x != 0).items():
            assert mult % 2 == 0

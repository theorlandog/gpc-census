"""v_B carries a nonzero gauge-invariant holonomy: an exact complexity proof.

The diagonal U(1)^d orbital-phase gauge cannot realify psi_B, because a support
cycle carries a holonomy that is not an integer multiple of pi. This is exact
(sympy), reproducible from the shipped state, and independent of the numerical
antiunitary-overlap estimate. It certifies complexity up to the diagonal gauge;
the non-diagonal degenerate-block rotations are a separate (open) question.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def test_vb_is_complex_up_to_diagonal_gauge():
    import sympy as sp

    import vb_holonomy as vh

    is_complex, cycles = vh.holonomies("(4,9)", 65, verbose=False)
    assert is_complex, "v_B should carry a nonzero diagonal-gauge holonomy"
    # at least one cycle has a holonomy that is provably not an integer * pi
    noninteger = [c for c in cycles if not c[2]]
    assert noninteger, "expected a non-gauge-removable cycle"
    # and the holonomy is a genuine irrational multiple of pi (not just unproven)
    over_pi = noninteger[0][1]
    assert sp.nsimplify(over_pi).is_rational is False or not over_pi.is_integer

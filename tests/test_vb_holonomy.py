"""The shipped v_B library state is real up to diagonal gauge.

The library previously shipped a phase-carrying fiber point for
v_B = (20,14,14,14,14,4,4,4,4)/23, whose support cycle carried a holonomy
provably not an integer multiple of pi. The record now ships the real
endpoint of its wall family (t0 = 7/17 - 18*sqrt(35)/85, squared weights in
Q(sqrt(35))), so the same exact holonomy computation must find every support
cycle at an integer multiple of pi: the shipped state is real, and its one
loop is a sign, not a phase.
"""
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def test_vb_library_state_is_real_up_to_diagonal_gauge():
    import vb_holonomy as vh

    is_complex, cycles = vh.holonomies("(4,9)", 65, verbose=False)
    assert not is_complex, "shipped v_B state must carry only trivial holonomy"
    # the support still carries its one-dimensional incidence kernel
    assert len(cycles) == 1, "v_B support must carry exactly one loop"
    # and the loop holonomy is an integer multiple of pi (a sign, not a phase)
    assert all(integral for _, _, integral in cycles)

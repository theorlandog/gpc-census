"""Differential guard: the fast signed-design enumerator must agree exactly
with the reference generic enumerator.

signed_design_fast.py adds three exact prunes (lone-term, pair-magnitude,
GF(2) sign propagation) on top of signed_design_generic.py. The prunes are
sound only if they never drop a real hit and never invent one, so the guard
is equivalence on spectra that actually PRODUCE hits (design and
signed-design), not just on a zero-hit spectrum like v96. See
docs/RESEARCH.md, "v96 campaign".
"""
import importlib.util
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _load(name):
    path = ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _summary(stats):
    hits = sorted(
        (h[0], tuple(sorted((tuple(d), w) for d, w in h[1])))
        for h in stats["hits"]
    )
    return stats["skel"], hits


# each case is small enough to enumerate fully; chosen so at least the
# design/signed-design cases carry hits that must survive every prune
CASES = [
    ((1, 1, 1, 1), 2),          # disjoint supports -> DESIGN hits
    ((2, 2, 1, 1), 2),          # one signed design
    ((2, 2, 2, 2), 2),
    ((3, 3, 2, 2, 1, 1), 2),
    ((2, 2, 1, 1, 1, 1), 2),
    ((1, 1, 1, 1, 1, 1), 3),
]


def test_fast_matches_generic():
    gen = _load("signed_design_generic")
    fast = _load("signed_design_fast")
    for ints, n in CASES:
        for do_phases in (False, True):
            a = _summary(gen.run(ints, n, do_phases=do_phases))
            b = _summary(fast.run(ints, n, do_phases=do_phases))
            assert a == b, f"divergence ints={ints} N={n} phases={do_phases}"
    # at least one case must actually exercise a hit, or the guard is vacuous
    assert _summary(gen.run((1, 1, 1, 1), 2, do_phases=False))[1]

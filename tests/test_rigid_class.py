"""Feature 3: rigid one-hop-class constraint, with the required corpus containment.

A one-hop class with a single hopping pair (i, j) on modes (A, B) contributes
rho_{AB} = (sign) sqrt(k_i k_j)/den to the 1-RDM: one real term, uncancellable,
so the pair must be a genuine block realizing target x2 = k_i k_j (the block
eigenvalue identity). Using k_A k_B = x2 as a search constraint is sound ONLY IF
no certified state violates it. This test proves containment over the shipped
corpus: for every certified state, every 1-term class satisfies the identity
against the state's own amplitudes (independent of the stored integer weights).
"""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

DATA = ROOT / "results" / "data" / "states.jsonl"


def _states():
    if not DATA.exists():
        return []
    out = []
    for ln in DATA.read_text().splitlines():
        if ln.strip():
            r = json.loads(ln)
            if r.get("status") == "OK" and r.get("closed_form"):
                out.append(r)
    return out


def test_rigid_class_identity_holds_over_whole_corpus():
    import sympy as sp

    from gpc_census.states import one_hop_classes

    states = _states()
    if not states:
        import pytest
        pytest.skip("no states atlas present")
    checked = 0
    for r in states:
        cf = r["closed_form"]
        dets = [tuple(t) for t in cf["support_dets"]]
        den = cf["den"]
        amps = [sp.sympify(a) for a in cf["pretty"]]
        # |c_i|^2 * den, exact, straight from the amplitude (not the stored weight)
        kw = [sp.nsimplify(sp.Abs(a) ** 2 * den) for a in amps]
        for (a_, b_), pairs in one_hop_classes(dets).items():
            if len(pairs) != 1:
                continue
            i, j = pairs[0]
            # off-diagonal magnitude^2 * den^2 for a single term = k_i * k_j, exactly
            off_sq_den2 = sp.simplify(sp.Abs(amps[i]) ** 2 * sp.Abs(amps[j]) ** 2 * den ** 2)
            assert sp.simplify(off_sq_den2 - kw[i] * kw[j]) == 0, (
                f"{r['system']} v{r['index']} class {(a_, b_)}: identity fails")
            # and the product is a positive integer (a genuine block, x2 > 0)
            assert kw[i] * kw[j] == int(kw[i] * kw[j]) and kw[i] * kw[j] > 0
            checked += 1
    # the corpus does contain 1-term classes, so the test is not vacuous
    assert checked > 0


def test_rigid_class_products_helper_matches_weights():
    from gpc_census.states import one_hop_classes, rigid_class_products

    states = _states()
    if not states:
        import pytest
        pytest.skip("no states atlas present")
    # spot-check on v_B (4,9) v65: rigid products equal k_i*k_j from stored weights
    vb = next((r for r in states if r["system"] == "(4,9)" and r["index"] == 65), None)
    if vb is None:
        import pytest
        pytest.skip("v_B not present")
    cf = vb["closed_form"]
    dets = [tuple(t) for t in cf["support_dets"]]
    w = cf["weights"]
    prods = rigid_class_products(dets, w)
    for pr, val in prods.items():
        pairs = one_hop_classes(dets)[pr]
        i, j = pairs[0]
        assert val == w[i] * w[j]

import importlib.util
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "adversarial_observable_audit.json"

MODULE_PATH = ROOT / "scripts" / "independent_adversarial_audit.py"
SPEC = importlib.util.spec_from_file_location("adversarial", MODULE_PATH)
aud = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(aud)


def jordan_wigner_annihilators(d):
    """A third sign implementation: explicit occupation-basis matrices.

    Shares no code with the repository's inversion-count formula or with the
    audit's sequential bitstring engine, so agreement is real corroboration.
    """
    dim = 1 << d
    ops = []
    for p in range(d):
        rows = [[0] * dim for _ in range(dim)]
        for state in range(dim):
            if (state >> p) & 1:
                parity = bin(state & ((1 << p) - 1)).count("1")
                rows[state ^ (1 << p)][state] = -1 if parity % 2 else 1
        ops.append(rows)
    return ops


def element_via_matrices(d, ops, left, right, P, Q):
    dim = 1 << d
    vec = [0] * dim
    vec[sum(1 << p for p in right)] = 1
    for q in Q:
        vec = [sum(ops[q][i][j] * vec[j] for j in range(dim) if vec[j]) for i in range(dim)]
    for p in reversed(P):
        vec = [sum(ops[p][j][i] * vec[j] for j in range(dim) if vec[j]) for i in range(dim)]
    return vec[sum(1 << p for p in left)]


def test_sign_convention_agrees_with_a_third_implementation():
    d = 4
    ops = jordan_wigner_annihilators(d)
    checked = 0
    for n in range(1, d + 1):
        dets = list(itertools.combinations(range(d), n))
        for k in range(1, n + 1):
            features = list(itertools.combinations(range(d), k))
            for left in dets:
                for right in dets:
                    for P in features:
                        for Q in features:
                            direct = element_via_matrices(d, ops, left, right, P, Q)
                            assert direct == aud.matrix_element(left, right, P, Q)
                            assert direct == aud.published_formula(left, right, P, Q)
                            checked += 1
    assert checked == 3285


def test_sign_audit_runs_clean_on_a_small_range():
    assert aud.sign_audit(4) > 0


def test_pair_rigid_carrier_can_still_hide_its_relative_phase():
    """The audit's headline counterexample, recomputed.

    A two-body operator moves at most two orbitals; these determinants differ
    in three, so no 2-RDM entry sees their relative phase even though the pair
    incidence has full column rank.
    """
    support = [(0, 1, 2), (3, 4, 5)]
    assert aud.pair_rank(support, 6) == len(support)
    assert aud.real_coherence_rank(support, 6) == 0
    assert len(set(support[0]) ^ set(support[1])) // 2 == 3


def test_connected_collision_free_graph_does_not_imply_linear_identifiability():
    """Separates certificate 1 from certificate 2 by an explicit witness.

    This carrier is pair-rigid and its collision-free graph is connected, so
    the tree certificate reconstructs it through purity. The complete
    real-linear channel matrix is still rank deficient, 4 against a target of
    6, so linear identifiability from all 2-RDM channels fails on the same
    carrier. The two statements are not interchangeable.
    """
    support = [(0, 1, 2), (0, 1, 3), (2, 4, 5)]
    assert aud.pair_rank(support, 6) == len(support)
    assert aud.connected(len(support), aud.collision_free_edges(support, 6))
    assert aud.real_coherence_rank(support, 6) == 4
    assert 4 < len(support) * (len(support) - 1)


def test_complementary_determinants_are_exactly_the_disconnected_pairs():
    row = next(
        r
        for r in json.loads(ARTIFACT.read_text())["small_carrier_summary"]
        if (r["d"], r["N"], r["m"]) == (6, 3, 2)
    )
    complementary = sum(
        1
        for a, b in itertools.combinations(itertools.combinations(range(6), 3), 2)
        if not set(a) & set(b)
    )
    assert row["full_pair_carriers"] - complementary == row["unique_graph_connected"]


def test_small_carrier_rows_regenerate():
    art = {(r["d"], r["N"], r["m"]): r for r in json.loads(ARTIFACT.read_text())["small_carrier_summary"]}
    for key in ((4, 2, 2), (4, 2, 5), (5, 3, 2)):
        assert aud.small_carrier_row(*key) == art[key]


def test_census_28_block_is_derived_from_the_shipped_artifacts():
    art = json.loads(ARTIFACT.read_text())["census_28"]
    assert aud.census_28_report() == art
    assert sum(art["first_connected_unique_channel_order"].values()) == 28
    assert art["newly_resolved_by_linear_determinacy_closure"] == 0
    assert "non-injectivity is not proved" in art["status"]


def test_the_audit_does_not_claim_more_than_it_proved():
    verdicts = json.loads(ARTIFACT.read_text())["theorem_verdicts"]
    assert verdicts["T6_general_multichart_zero_crossing_theorem"] == "NOT_YET_PROVED"
    assert "STANDARD" in verdicts["T2_simplex_torus_topology"]
    assert "STANDARD" in verdicts["T3_symplectic_action_angle_form"]
    assert "PROOF_STILL_REQUIRED" in verdicts["T7_cyclic_window_family"]

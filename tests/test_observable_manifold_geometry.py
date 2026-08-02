import importlib.util
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).parents[1]
DATA = ROOT / "tests" / "data"
if not DATA.exists():
    DATA = ROOT / "results" / "data"
ARTIFACT = DATA / "observable_manifold_geometry.json"
MODULE_PATH = Path(__file__).parents[1] / "scripts" / "observable_manifold_geometry.py"
SPEC = importlib.util.spec_from_file_location("observable_manifold_geometry", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_full_support_tree_chart_is_connected():
    for m in range(2, 10):
        assert MODULE.component_count(frozenset(range(m)), MODULE.path_tree(m)) == 1
        assert MODULE.component_count(frozenset(range(m)), MODULE.star_tree(m)) == 1


def test_boundary_ambiguity_is_components_minus_one():
    edges = MODULE.path_tree(5)
    active = frozenset({0, 2, 4})
    assert MODULE.component_count(active, edges) == 3
    audit = MODULE.audit_tree(5, edges, "path")
    assert audit["worst_boundary_stratum"]["ambiguity_dimension"] == 2


def test_star_root_loss_exposes_large_chart_singularity():
    m = 7
    edges = MODULE.star_tree(m)
    active = frozenset(range(1, m))
    assert MODULE.component_count(active, edges) == m - 1


def test_shipped_artifact_matches_a_fresh_build():
    assert json.loads(ARTIFACT.read_text()) == MODULE.build(10)


def _counts(audit):
    return {
        int(size): {int(k): v for k, v in blk["ambiguity_dimension_counts"].items()}
        for size, blk in audit["boundary"].items()
    }


def test_path_and_star_strata_match_closed_form_combinatorics():
    """Independent count, derived by hand rather than by rerunning the audit.

    On a path the induced forest components are the maximal runs of
    consecutive active indices, and there are C(s-1,r-1)*C(m-s+1,r) supports
    of size s with r runs.  On a star the chart is connected exactly when the
    centre is active, so the c(A)-1 spectrum is 0 with multiplicity
    C(m-1,s-1) and s-1 with multiplicity C(m-1,s).
    """
    audits = {(a["name"], a["vertices"]): a for a in json.loads(ARTIFACT.read_text())["audits"]}
    for m in range(2, 11):
        path = {
            s: {
                r - 1: comb(s - 1, r - 1) * comb(m - s + 1, r)
                for r in range(1, s + 1)
                if comb(s - 1, r - 1) * comb(m - s + 1, r)
            }
            for s in range(1, m + 1)
        }
        assert _counts(audits[("path", m)]) == path, m

        star = {}
        for s in range(1, m + 1):
            block = {}
            if comb(m - 1, s - 1):
                block[0] = comb(m - 1, s - 1)
            if s <= m - 1 and comb(m - 1, s):
                block[s - 1] = block.get(s - 1, 0) + comb(m - 1, s)
            star[s] = block
        assert _counts(audits[("star", m)]) == star, m


def test_every_audit_covers_all_supports_at_the_projective_dimension():
    for audit in json.loads(ARTIFACT.read_text())["audits"]:
        m = audit["vertices"]
        assert audit["full_support_dimension_real"] == 2 * m - 2
        assert audit["full_support_chart_connected"] is True
        total = audit["all_nonempty_strata_ambiguity_counts"]
        assert sum(total.values()) == 2**m - 1
        for size, block in audit["boundary"].items():
            assert block["strata"] == comb(m, int(size))

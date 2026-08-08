#!/usr/bin/env python3
"""The theorem-dependency graph for rank-independent GPC compression.

Emits `results/data/infinite_wedge_tca_dependency.json`: machine-readable nodes
and arrows for the chain

    stable algebra A_N = Sym(wedge^N C^infinity)
        -> finite-d coordinate representation Sym(wedge^N C^d)
        -> highest-weight support semigroup Sigma_{N,d}
        -> rational moment cone
        -> trace slice and dominant chamber
        -> GPC polytope
        -> supporting hyperplanes / facets
        -> finite facet templates across d

Every arrow carries exactly one status from

    STANDARD THEOREM         cited, holds off the shelf
    PROVED IN THIS WORK      proved in docs/infinite_wedge_tca_bridge.md
    FINITE COMPUTATIONAL CHECK   verified on the census window only
    MISSING THEOREM          the statement the program needs and does not have
    FALSE                    refuted, with a named exact counterexample

Live verdict fields are pulled from `results/data/infinite_wedge_tca_sanity.json`
rather than typed, so the graph cannot drift from the measurement.

Run: uv run scripts/infinite_wedge_tca_dependency.py
"""
from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
SANITY = DATA / "infinite_wedge_tca_sanity.json"
OUT = DATA / "infinite_wedge_tca_dependency.json"

# Literature. `checked` records how far the citation was verified in this
# session: "abstract" means the statement was read from the publisher or search
# abstract, not from the full text, because the network policy of the drafting
# environment blocks the primary hosts. A human must confirm the precise
# statements before any of this is quoted outside the repository.
SOURCES = {
    "hadziev-grosshans": {
        "what": "C[X]^U is finitely generated for an affine G-variety X, G "
                "reductive, U a maximal unipotent subgroup",
        "where": "Hadziev 1967; Grosshans, Algebraic Homogeneous Spaces and "
                 "Invariant Theory, LNM 1673",
        "checked": "abstract",
    },
    "mumford-ness-brion": {
        "what": "the moment polytope of a projective G-variety is the closure "
                "of the normalized highest weights occurring in its homogeneous "
                "coordinate ring",
        "where": "Mumford/Ness; Brion, Sur l'image de l'application moment; "
                 "used in this project as Klyachko's practical algorithm",
        "checked": "abstract, plus the exact finite check in this campaign",
    },
    "vergne-walter": {
        "what": "the moment cone of an arbitrary finite-dimensional unitary "
                "representation of a compact connected Lie group is cut out by "
                "finitely many linear inequalities, via a variant of Ressayre's "
                "dominant pairs",
        "where": "Vergne and Walter, J. Symplectic Geom. 15 (2017); "
                 "arXiv:1410.8144",
        "checked": "abstract",
    },
    "bulois-denis-ressayre": {
        "what": "an algorithm computing the MINIMAL inequality list of the "
                "moment cone of any representation of a complex reductive "
                "group, with the fermionic cone GL_d acting on wedge^r C^d "
                "implemented explicitly",
        "where": "Bulois, Denis and Ressayre, An Algorithm to compute the "
                 "Kronecker cone and other moment cones, 2025; arXiv:2505.08812",
        "checked": "abstract",
    },
    "draisma-polynomial-functors": {
        "what": "every finite-degree polynomial functor over an infinite field "
                "is topologically noetherian, so a GL-stable closed subset is "
                "cut out by finitely many orbits of equations",
        "where": "Draisma, J. Amer. Math. Soc. 32 (2019) 691-707",
        "checked": "abstract",
    },
    "draisma-eggermont": {
        "what": "bounded Plucker varieties are defined in bounded degree; the "
                "limit in the dual infinite wedge is cut out, up to symmetry, "
                "by finitely many equations",
        "where": "Draisma and Eggermont, J. reine angew. Math. 737 (2018) "
                 "189-215",
        "checked": "abstract",
    },
    "nekrasov": {
        "what": "the dual infinite wedge is GL_infinity-equivariantly "
                "noetherian, extending the bounded Plucker case",
        "where": "Nekrasov, arXiv:2008.09531 (2020)",
        "checked": "abstract",
    },
    "nagpal-sam-snowden": {
        "what": "among unbounded twisted commutative algebras, noetherianity is "
                "known ONLY for Sym(Sym^2 C^infinity) and Sym(wedge^2 "
                "C^infinity); the degree three and higher cases, which include "
                "A_3 = Sym(wedge^3 C^infinity), are open",
        "where": "Nagpal, Sam and Snowden, Selecta Math. 22 (2016); and the "
                 "skew-commutative sequel, Selecta Math. (2019)",
        "checked": "abstract",
    },
    "fi-oi-noetherianity": {
        "what": "finitely generated FI-modules and OI-modules over noetherian "
                "rings are noetherian, and finitely generated symmetric ideals "
                "in infinite polynomial rings are finitely generated up to "
                "symmetry",
        "where": "Church-Ellenberg-Farb; Aschenbrenner-Hillar; "
                 "Hillar-Sullivant, Adv. Math. 229 (2012); Draisma, "
                 "Noetherianity up to symmetry, LNM 2108 (2014)",
        "checked": "abstract",
    },
    "burgisser-et-al": {
        "what": "membership in a moment polytope is in NP and coNP, and there "
                "is a polynomial-time weak membership oracle by scaling, for a "
                "family described as having exponentially many facets",
        "where": "Burgisser-Christandl-Mulmuley-Walter, SIAM J. Comput.; "
                 "Burgisser-Franks-Garg-Oliveira-Walter-Wigderson, FOCS 2018",
        "checked": "abstract; the exponential facet count is asserted for "
                   "moment polytopes broadly and was NOT verified for the "
                   "fermionic family specifically",
    },
    "horn-klyachko-ktw": {
        "what": "the Horn cone for sums of Hermitian matrices is cut out by the "
                "Horn inequalities indexed by triples (I,J,K) of subsets of "
                "size r for every r < n, and the essential ones are exactly "
                "those with Littlewood-Richardson coefficient 1. Their arity "
                "grows with r, hence is UNBOUNDED in n",
        "where": "Klyachko, Selecta Math. 4 (1998); Knutson-Tao, J. Amer. Math. "
                 "Soc. 12 (1999); Knutson-Tao-Woodward, J. Amer. Math. Soc. 17 "
                 "(2004)",
        "checked": "standard, from memory of the statement; the arity claim is "
                   "read off the shape of the inequalities, not off a growth "
                   "theorem",
        "why_here": "the closest moment cone whose facets are COMPLETELY known, "
                    "and it fails bounded-arity templating. Absent a reason why "
                    "the fermionic cone should be better behaved, T1 should be "
                    "expected false",
    },
    "reuvers": {
        "what": "the GPC polytope volume approaches the Pauli hypersimplex "
                "volume as the one-body dimension grows, so the generalized "
                "constraints are volumetrically negligible except in few-fermion "
                "settings",
        "where": "Reuvers, J. Math. Phys. 62, 032204 (2021); arXiv:1911.00471",
        "checked": "abstract",
    },
    "repo-census-sandwich": {
        "what": "P = O at all nine census systems, so validity on the certified "
                "vertex set is validity on the polytope",
        "where": "docs/census_inner_sandwich.md; "
                 "results/data/census_inner_sandwich.json",
        "checked": "in-repo artifact",
    },
    "repo-quantization": {
        "what": "the quantization period of a census vertex is the order of a "
                "finite isotropy character of its certified attainer; p = q at "
                "ten of fourteen measured vertices and 2q at the other four",
        "where": "docs/quantization_character.md; docs/symplectic_invariants.md; "
                 "results/data/quantization_multiplicities.json",
        "checked": "in-repo artifact",
    },
    "repo-screening-orbit": {
        "what": "screening is NOT a Weyl-orbit invariant: the dominant chamber "
                "is a fundamental domain and so is exactly what the Weyl group "
                "does not preserve, and the moment cone's Weyl invariance is a "
                "statement about the AMBIENT cone, not a chamber-facing list",
        "where": "docs/screening_orbit_structure.md",
        "checked": "in-repo artifact",
    },
}

NODES = [
    {"id": "A_N", "label": "stable algebra A_N = Sym(wedge^N C^infinity)",
     "kind": "twisted commutative algebra / polynomial functor",
     "note": "the object the compression program proposes to work with"},
    {"id": "A_N_d", "label": "finite-rank coordinate representation Sym(wedge^N C^d)",
     "kind": "graded GL_d-algebra"},
    {"id": "Sigma_stab", "label": "stable highest-weight support semigroup Sigma_N^stab",
     "kind": "sub-monoid of {partitions} x Z_{>0}",
     "definition": "{(lambda,k) : <s_lambda, h_k[e_N]> > 0}, no d anywhere"},
    {"id": "Sigma_d", "label": "finite-rank semigroup Sigma_{N,d}",
     "kind": "sub-monoid of the dominant weight lattice x Z_{>0}",
     "definition": "{(lambda,k) : ell(lambda) <= d, mult(V_lambda, "
                   "Sym^k(wedge^N C^d)) > 0}"},
    {"id": "cone_d", "label": "rational moment cone cone(Sigma_{N,d})",
     "kind": "rational polyhedral cone in R^d x R"},
    {"id": "Delta_d", "label": "trace slice and dominant chamber: the GPC polytope Delta(N,d)",
     "kind": "polytope, the k = 1 normalized slice"},
    {"id": "facets_d", "label": "supporting hyperplanes and facets of Delta(N,d)",
     "kind": "finite list of rational inequalities, the GPC rows"},
    {"id": "templates", "label": "finite facet templates across d",
     "kind": "the object rank-independent compression would need",
     "note": "a finite list of coefficient patterns plus a decision rule saying "
             "which instances are facets at rank d"},
]


def _sanity() -> dict:
    if not SANITY.exists():
        raise SystemExit(f"missing {SANITY}; run scripts/infinite_wedge_tca_sanity.py first")
    return json.loads(SANITY.read_text())


def arrows(sanity: dict) -> list[dict]:
    score = sanity["scorecard"]
    padding = sanity["padding"]["pairs"]
    insertion = sanity["insertion_post_hoc"]["pairs"]
    restriction = sanity["restriction"]["pairs"]
    templates = sanity["templates"]["per_system"]
    saturation = sanity["saturation"]

    broken_padding = [f"{p['from']} -> {p['to']}: "
                      f"{p['rows_tested'] - p['rows_still_valid']} of "
                      f"{p['rows_tested']} rows invalid"
                      for p in padding if not p["padding_preserves_validity"]]
    broken_insertion = [f"{p['from']} -> {p['to']}: "
                        f"{p['rows_tested'] - p['rows_with_some_valid_embedding']} "
                        f"of {p['rows_tested']} rows have NO valid embedding"
                        for p in insertion
                        if not p["every_row_has_some_valid_embedding"]]

    return [
        {
            "id": "A1", "from": "A_N", "to": "A_N_d",
            "statement": "evaluate the polynomial functor at C^d",
            "status": "STANDARD THEOREM",
            "sources": ["draisma-polynomial-functors"],
            "obstruction": None,
        },
        {
            "id": "A2", "from": "A_N_d", "to": "Sigma_d",
            "statement": "Sigma_{N,d} is the weight monoid of the U-invariant "
                         "subalgebra Sym(wedge^N C^d)^U, graded by degree; it is "
                         "a monoid because Sym is a domain and the product of "
                         "two highest-weight vectors is a nonzero highest-weight "
                         "vector of the sum weight, and it is FINITELY GENERATED "
                         "by Hadziev-Grosshans",
            "status": "STANDARD THEOREM",
            "sources": ["hadziev-grosshans"],
            "obstruction": None,
        },
        {
            "id": "A3", "from": "Sigma_stab", "to": "Sigma_d",
            "statement": "Sigma_{N,d} = Sigma_N^stab intersect {ell(lambda) <= d}. "
                         "The multiplicity of V_lambda in Sym^k(wedge^N C^d) is "
                         "the plethysm coefficient <s_lambda, h_k[e_N]>, which "
                         "does not depend on d, and S_lambda(C^d) = 0 exactly "
                         "when ell(lambda) > d. So the whole family is the "
                         "length truncation of ONE d-independent monoid",
            "status": "PROVED IN THIS WORK",
            "also": "FINITE COMPUTATIONAL CHECK",
            "check": {"verdict": score["P1_three_routes_agree"]["verdict"],
                      "pairs": score["P1_three_routes_agree"]["scope"],
                      "control": score["P2_convention_is_load_bearing"]["verdict"]},
            "sources": ["repo-quantization"],
            "obstruction": None,
            "note": "this is the honest content of 'one stable object'. It is "
                    "elementary, and it is why the stable object is NOT where "
                    "the difficulty lives",
        },
        {
            "id": "A4", "from": "Sigma_d", "to": "cone_d",
            "statement": "a finitely generated monoid spans a rational "
                         "polyhedral cone, and {ell <= d} is a face of the "
                         "dominant cone, so cone(Sigma_{N,d}) = cone(Sigma_N^stab) "
                         "intersect {lambda_{d+1} = 0}",
            "status": "STANDARD THEOREM",
            "sources": ["hadziev-grosshans"],
            "obstruction": None,
        },
        {
            "id": "A4rev", "from": "cone_d", "to": "Sigma_d",
            "statement": "the cone determines the semigroup, that is, "
                         "Sigma_{N,d} is saturated",
            "status": "FALSE",
            "counterexample": {
                "witnesses": [w["vertex"] for w in
                              saturation["non_saturation_witnesses"]],
                "smallest": (saturation["non_saturation_witnesses"][0]
                             if saturation["non_saturation_witnesses"] else None),
                "predicted_census_wide": saturation["predicted_census_wide_witnesses"],
            },
            "sources": ["repo-quantization"],
            "obstruction": "the quantization period exceeds the spectrum "
                           "denominator at these vertices, so (lambda_int, q) is "
                           "a dominant lattice point of the cone with "
                           "multiplicity zero",
            "consequence": "finite generation of the SEMIGROUP is strictly "
                           "stronger than finite generation of the CONE, and "
                           "only the cone is needed for facets. Any compression "
                           "theorem must say which one it is about",
        },
        {
            "id": "A5", "from": "cone_d", "to": "Delta_d",
            "statement": "the normalized k = 1 slice of cone(Sigma_{N,d}) is the "
                         "fermionic moment polytope",
            "status": "STANDARD THEOREM",
            "also": "FINITE COMPUTATIONAL CHECK",
            "check": {"systems": [
                {"system": r["system"], "max_degree": r["max_degree"],
                 "generated_points": r["generated_points"],
                 "hull_equals_polytope": r["hull_equals_polytope"]}
                for r in sanity["recovery"]["systems"]]},
            "sources": ["mumford-ness-brion", "repo-census-sandwich"],
            "obstruction": None,
        },
        {
            "id": "A6", "from": "Delta_d", "to": "facets_d",
            "statement": "the moment cone of a FIXED representation has a finite "
                         "irredundant inequality list, computable",
            "status": "STANDARD THEOREM",
            "sources": ["vergne-walter", "bulois-denis-ressayre"],
            "obstruction": None,
            "note": "this is finiteness at each fixed d and carries no "
                    "uniformity in d whatsoever",
        },
        {
            "id": "A7", "from": "facets_d", "to": "facets_d",
            "direction": "downward, rank d+1 to rank d",
            "statement": "Delta(N,d) = Delta(N,d+1) intersect {lambda_{d+1} = 0}, "
                         "because an orbital of occupation zero is annihilated by "
                         "its own annihilation operator; hence restricting the "
                         "rank d+1 system to lambda_{d+1} = 0 cuts out exactly "
                         "Delta(N,d)",
            "status": "PROVED IN THIS WORK",
            "also": "FINITE COMPUTATIONAL CHECK",
            "check": {"verdict": score["P5_restriction_is_complete"]["verdict"],
                      "pairs": score["P5_restriction_is_complete"]["pairs"],
                      "method": "exact affine Farkas over the rationals"},
            "sources": ["repo-census-sandwich"],
            "obstruction": None,
            "consequence": "compression DOWN the rank ladder is a theorem: the "
                           "published (3,10) system implies the (3,9), (3,8), "
                           "(3,7) and (3,6) systems, and (4,10) implies (4,9) "
                           "and (4,8)",
        },
        {
            "id": "A8", "from": "facets_d", "to": "facets_d",
            "direction": "upward by zero padding, rank d to rank d+1",
            "statement": "a valid row at rank d, padded with a zero coefficient, "
                         "is valid at rank d+1",
            "status": "FALSE",
            "counterexample": {"broken_pairs": broken_padding,
                               "verdict": score[
                                   "P4_zero_padding_preserves_validity_above_rank_6"
                               ]["verdict"]},
            "sources": ["repo-census-sandwich"],
            "obstruction": "a published row names orbitals by POSITION relative "
                           "to d, so padding silently re-points it at a "
                           "different orbital",
        },
        {
            "id": "A8oi", "from": "facets_d", "to": "facets_d",
            "direction": "upward by ANY order-preserving embedding",
            "statement": "for a valid row at rank d there is SOME order-"
                         "preserving injection [d] -> [d+1] carrying it to a "
                         "valid row",
            "status": "FALSE",
            "counterexample": {"broken_pairs": broken_insertion},
            "sources": ["repo-screening-orbit"],
            "obstruction": "the inequality side of the census carries no "
                           "covariant OI structure at all; it is a purely "
                           "contravariant (restriction) functor",
            "note": "POST-HOC block, added after the pre-registered zero-padding "
                    "prediction failed",
        },
        {
            "id": "A9", "from": "facets_d", "to": "templates",
            "statement": "there is a finite list of templates, and a decision "
                         "rule of bounded complexity, that together produce the "
                         "facet list of Delta(N,d) for every d",
            "status": "MISSING THEOREM",
            "sources": [],
            "obstruction": "no known theorem supplies it; see the audit arrows "
                           "below, all of which fail for the same reason",
            "measured_symptoms": {
                "rows_per_system": {r["system"]: r["published_inequalities"]
                                    for r in templates},
                "distinct_symmetric_templates": {
                    r["system"]: r["distinct_symmetric_templates"]
                    for r in templates},
                "max_symmetric_arity": {r["system"]: r["max_symmetric_arity"]
                                        for r in templates},
                "new_symmetric_templates": {
                    r["system"]: r["new_symmetric_templates_vs_previous_rank"]
                    for r in templates},
                "reading": "the symmetric template stock is far smaller than the "
                           "row count, which is encouraging, but it does not "
                           "saturate inside the census window and the templates "
                           "alone do not say WHICH instances are facets",
            },
        },
        {
            "id": "L1", "from": "templates", "to": "facets_d",
            "statement": "if every facet of the ambient Weyl-invariant moment "
                         "cone at every rank is an instance under S_d of a fixed "
                         "finite template list of maximum symmetric arity a, "
                         "then that cone has O(d^a) facets, since an orbit of an "
                         "arity-a template has at most binom(d,a) * a! members",
            "status": "PROVED IN THIS WORK",
            "sources": [],
            "obstruction": None,
            "sources_for_the_prior": ["horn-klyachko-ktw", "burgisser-et-al"],
            "consequence": "contrapositive: superpolynomial facet growth in d "
                           "REFUTES finite templating. This converts the "
                           "compression question into a facet-counting growth "
                           "question, which is the smallest decisive experiment",
            "caution": "the literature describes moment polytopes as having "
                       "exponentially many facets, but that assertion was not "
                       "verified for the fermionic family in this session",
        },
        {
            "id": "X1", "from": "draisma-polynomial-functors", "to": "A9",
            "statement": "topological noetherianity of polynomial functors "
                         "implies finite facet templates",
            "status": "FALSE",
            "kind": "scope audit",
            "sources": ["draisma-polynomial-functors"],
            "object_covered": "wedge^N C^infinity as a polynomial functor: YES",
            "finiteness_proved": "descending chains of GL-stable Zariski-closed "
                                 "subsets stabilize; equations up to orbit",
            "covers_A_N": True,
            "covers_the_weight_semigroup": False,
            "facet_orbit_finiteness_follows": False,
            "obstruction": "a facet of a moment cone is a real convex-geometric "
                           "datum, not a Zariski-closed condition, and the family "
                           "d -> Delta(N,d) is not a descending chain of closed "
                           "subsets of anything. Finitely many equations up to "
                           "symmetry says nothing about finitely many "
                           "inequalities up to symmetry",
        },
        {
            "id": "X2", "from": "draisma-eggermont", "to": "A9",
            "statement": "bounded Plucker varieties and the infinite wedge give "
                         "finite facet templates",
            "status": "FALSE",
            "kind": "scope audit",
            "sources": ["draisma-eggermont", "nekrasov"],
            "object_covered": "subvarieties of the infinite wedge closed under "
                              "the Plucker operations: YES, and this is the "
                              "closest existing result to the stable object",
            "finiteness_proved": "ideals generated in bounded degree up to "
                                 "GL_infinity",
            "covers_A_N": True,
            "covers_the_weight_semigroup": False,
            "facet_orbit_finiteness_follows": False,
            "obstruction": "same as X1. In addition the fermionic moment "
                           "polytope is the moment image of the FULL projective "
                           "space P(wedge^N C^d), not of a proper subvariety, so "
                           "there is no interesting Plucker variety in play here "
                           "at all; the equations these theorems bound are "
                           "equations of objects the GPC problem does not use",
        },
        {
            "id": "X3", "from": "fi-oi-noetherianity", "to": "A9",
            "statement": "FI or OI noetherianity gives finite facet templates",
            "status": "FALSE",
            "kind": "scope audit",
            "sources": ["fi-oi-noetherianity"],
            "object_covered": "finitely generated FI-modules and OI-modules",
            "finiteness_proved": "submodules of finitely generated modules are "
                                 "finitely generated",
            "covers_A_N": "not directly",
            "covers_the_weight_semigroup": False,
            "facet_orbit_finiteness_follows": False,
            "obstruction": "two independent failures. First, module generation "
                           "allows arbitrary R-linear combinations while a cone "
                           "needs NONNEGATIVE ones, so finite module generation "
                           "of the span of the facet normals does not bound the "
                           "extreme rays of the dual cone. Second, and fatally, "
                           "arrow A8oi shows the covariant OI action on the "
                           "inequality side does not exist, so there is no "
                           "OI-module to which the theorem could be applied",
        },
        {
            "id": "X4", "from": "nagpal-sam-snowden", "to": "A9",
            "statement": "noetherianity of the tca A_N gives finite facet "
                         "templates",
            "status": "FALSE",
            "kind": "scope audit",
            "sources": ["nagpal-sam-snowden"],
            "object_covered": "unbounded tcas; known only for Sym(Sym^2) and "
                              "Sym(wedge^2)",
            "finiteness_proved": "module- and ideal-theoretic noetherianity",
            "covers_A_N": "A_2 yes, A_N for N >= 3 OPEN",
            "covers_the_weight_semigroup": False,
            "facet_orbit_finiteness_follows": False,
            "obstruction": "modules and ideals, not cones; the same gap as X3's "
                           "first failure. Note the coincidence worth recording: "
                           "the one fermionic case where tca noetherianity IS "
                           "known, N = 2, is exactly the case where the GPC "
                           "answer is classically trivial (lambda_{2i-1} = "
                           "lambda_{2i}), and the project's main series N = 3 "
                           "sits precisely in the open range",
        },
        {
            "id": "X5", "from": "vergne-walter", "to": "A9",
            "statement": "the fixed-representation finite inequality list is "
                         "uniform in d",
            "status": "FALSE",
            "kind": "scope audit",
            "sources": ["vergne-walter", "bulois-denis-ressayre"],
            "object_covered": "the moment cone of one fixed representation",
            "finiteness_proved": "a finite irredundant inequality list at each "
                                 "fixed d",
            "covers_A_N": "one d at a time",
            "covers_the_weight_semigroup": "the cone, yes; the semigroup, no",
            "facet_orbit_finiteness_follows": False,
            "obstruction": "the theorem quantifies over one representation. "
                           "Arrows A8 and A8oi show the output is not even "
                           "monotone in d for these tables, so no uniformity can "
                           "be read off it",
        },
        {
            "id": "X6", "from": "burgisser-et-al", "to": "templates",
            "statement": "rank-20 or rank-40 usability requires the facet list",
            "status": "FALSE",
            "kind": "scope audit, and the practical alternative",
            "sources": ["burgisser-et-al", "reuvers"],
            "obstruction": "membership in a moment polytope is already in NP and "
                           "coNP, with a polynomial-time weak membership oracle "
                           "by scaling. So the physics question 'is this "
                           "spectrum admissible at rank 40' does not need a "
                           "facet list, and the compression program should say "
                           "explicitly which of the two goals it is pursuing",
        },
    ]


def main() -> int:
    sanity = _sanity()
    payload = {
        "generated_by": "scripts/infinite_wedge_tca_dependency.py",
        "note": "docs/infinite_wedge_tca_bridge.md",
        "prereg": "docs/prereg_infinite_wedge_tca_sanity.md",
        "measurement": "results/data/infinite_wedge_tca_sanity.json",
        "evidence": "mixed: cited theorems, elementary proofs, and exact finite "
                    "checks over the census window; every arrow says which",
        "status_legend": {
            "STANDARD THEOREM": "cited, holds off the shelf",
            "PROVED IN THIS WORK": "proved in docs/infinite_wedge_tca_bridge.md",
            "FINITE COMPUTATIONAL CHECK": "verified on the census window only",
            "MISSING THEOREM": "the statement the program needs and lacks",
            "FALSE": "refuted, with a named exact counterexample",
        },
        "citation_caveat": "the network policy of the drafting environment "
                           "blocked the primary hosts, so every literature entry "
                           "was verified to abstract level only. Confirm the "
                           "precise statements before quoting this outside the "
                           "repository.",
        "sources": SOURCES,
        "nodes": NODES,
        "arrows": arrows(sanity),
    }
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")
    for arrow in payload["arrows"]:
        print(f"  {arrow['id']:6s} {arrow['status']:26s} "
              f"{arrow['from']} -> {arrow['to']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Cross-realization invariant census: lattice cones and exact stabilizers.

Computes, for all 799 census vertices and their certified states:

  LATTICE       the tangent cone at each vertex in the weight lattice, with
                its primitive edge rays, Smith invariant factors, lattice
                index and unimodularity;
  SYMPLECTIC    the exact projective stabilizer `Stab_{U(d)}([psi])`, the orbit
                dimension, the commutant of the state's own 1-RDM and the
                dimension of its orbit through the state.

and scores the pre-registered predictions of
`docs/prereg_symplectic_invariants.md`.

Everything is exact. Ranks are computed over algebraic number fields, never
numerically, and lattice data is integral throughout.

Run: uv run scripts/symplectic_invariant_census.py [--resume]
Writes results/data/symplectic_invariants.json.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time
from fractions import Fraction as F


ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "results" / "data"
OUT = DATA / "symplectic_invariants.json"
CHECKPOINT = ROOT / "build" / "symplectic_invariants_progress.jsonl"

SYSTEMS = [(3, 6), (3, 7), (3, 8), (4, 8), (3, 9), (4, 9), (3, 10), (4, 10), (5, 10)]

from gpc_census.symplectic import (  # noqa: E402
    affine_dim,
    facet_supports,
    hull_lattice,
    stabilizer_invariants,
    support_invariant,
    tangent_cone_invariants,
    toral_stabilizer_dim_from_support,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resume", action="store_true", help="reuse the checkpoint file")
    args = ap.parse_args()

    states = [json.loads(line) for line in (DATA / "states.jsonl").read_text().splitlines()]
    by_key = {(r["system"], r["index"]): r for r in states}

    done: dict[str, dict] = {}
    if args.resume and CHECKPOINT.exists():
        for line in CHECKPOINT.read_text().splitlines():
            row = json.loads(line)
            done[f"{row['system']}#{row['index']}"] = row
        print(f"resuming with {len(done)} rows already computed")
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    started = time.time()
    with CHECKPOINT.open("a") as ck:
        for n, d in SYSTEMS:
            system = f"({n},{d})"
            verts = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())
            points = [[F(x) for x in v["spectrum"]] for v in verts]
            dim = affine_dim(points)
            facets = facet_supports(points, n, d, dim)
            lattice = hull_lattice(points, d)
            assert len(lattice) == dim, (system, len(lattice), dim)
            print(f"{system}: dim {dim}, {len(points)} vertices, {len(facets)} facets",
                  flush=True)

            for i, v in enumerate(verts):
                key = f"{system}#{i}"
                record = by_key[(system, i)]
                if key in done:
                    # Only the exact stabilizer block is expensive enough to cache.
                    # The purely combinatorial fields are recomputed every run so a
                    # change to their definition cannot be masked by a stale
                    # checkpoint.
                    cached = done[key]
                    cdets = record["closed_form"]["support_dets"]
                    cached["support_invariant"] = support_invariant(cdets, d)
                    cached["toral_stab_dim_predicted"] = (
                        toral_stabilizer_dim_from_support(cdets, d))
                    rows.append(cached)
                    continue
                assert list(record["integer_form"]) == list(v["integer_form"]), key
                t0 = time.time()
                cone = tangent_cone_invariants(points, i, facets, lattice, dim)
                stab = stabilizer_invariants(record)
                dets = record["closed_form"]["support_dets"]
                row = {
                    "system": system,
                    "index": i,
                    "n": n,
                    "d": d,
                    "dim": dim,
                    "classified": record["classified"],
                    "integer_form": list(record["integer_form"]),
                    "denominator": record["denominator"],
                    "support_size": len(dets),
                    "support_invariant": support_invariant(dets, d),
                    "evidence": "exact",
                    "cone": cone,
                    "stabilizer": stab,
                    "toral_stab_dim_predicted": toral_stabilizer_dim_from_support(dets, d),
                    "secs": round(time.time() - t0, 2),
                }
                rows.append(row)
                ck.write(json.dumps(row) + "\n")
                ck.flush()
                print(f"  {key} stab {stab['stab_dim']:>3} orbit {stab['orbit_dim']:>3} "
                      f"cone {'simplicial' if cone['simplicial_tangent_cone'] else 'non-simp'}"
                      f" idx {cone['lattice_index']} ({row['secs']}s)", flush=True)

    payload = {
        "generated_by": "scripts/symplectic_invariant_census.py",
        "prereg": "docs/prereg_symplectic_invariants.md",
        "evidence": "exact",
        "note": (
            "Lattice-cone and exact-stabilizer invariants of every census vertex and "
            "its certified state. Stabilizer numbers are UNREDUCED and FIXED-RHO: the "
            "commutant is that of the state's own 1-RDM. No first-order fiber tangent "
            "dimension is reported; see the prereg for why such a number is vacuous "
            "at these ranks."
        ),
        "rows": rows,
        "scorecard": scorecard(rows),
        "wall_secs": round(time.time() - started, 1),
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n")
    print(f"\nwrote {OUT.relative_to(ROOT)} ({len(rows)} rows)")
    for name, verdict in payload["scorecard"].items():
        if isinstance(verdict, dict) and "verdict" in verdict:
            print(f"  {name}: {verdict['verdict']}")
    return 0


def scorecard(rows: list[dict]) -> dict:
    """Score the pre-registered predictions from the rows alone."""
    classes = {r["support_invariant"] for r in rows}
    scope = {"scope": len(rows), "independent_scope": len(classes),
             "independent_distinct_states": len(classes)}

    diagonal = sum(1 for r in rows if r["stabilizer"]["rho_is_diagonal"])
    p1 = {"verdict": "PASS" if diagonal == 645 else "FAIL",
          "predicted_diagonal": 645, "measured_diagonal": diagonal,
          "measured_non_diagonal": len(rows) - diagonal, **scope}

    slater, p2_bad = [], []
    for r in rows:
        den = r["denominator"] or 1
        if set(r["integer_form"]) <= {0, den}:
            slater.append(r)
            want = r["n"] ** 2 + (r["d"] - r["n"]) ** 2
            if r["stabilizer"]["stab_dim"] != want or \
                    r["stabilizer"]["orbit_dim"] != 2 * r["n"] * (r["d"] - r["n"]):
                p2_bad.append(f"{r['system']}#{r['index']}")
    p2 = {"verdict": "PASS" if slater and not p2_bad else "FAIL",
          "slater_rows": len(slater), "violations": p2_bad,
          "scope": len(slater), "independent_scope": len(slater),
          "independent_distinct_states": len(slater)}

    stabs = [r["stabilizer"]["stab_dim"] for r in rows]
    p3 = {"verdict": "PASS" if min(stabs) > 1 else "FAIL",
          "min_stab_dim": min(stabs), "max_stab_dim": max(stabs),
          "rows_at_minimum": sum(1 for s in stabs if s == min(stabs)), **scope}

    p4_bad = [f"{r['system']}#{r['index']}" for r in rows
              if r["stabilizer"]["toral_stab_dim"] != r["toral_stab_dim_predicted"]]
    p4 = {"verdict": "PASS" if not p4_bad else "FAIL",
          "violations": p4_bad, "violation_count": len(p4_bad), **scope}

    simplicial = sum(1 for r in rows if r["cone"]["simplicial_tangent_cone"])
    simple = sum(1 for r in rows if r["cone"]["simple"])
    unimodular = sum(1 for r in rows if r["cone"]["unimodular"])
    mismatch = [f"{r['system']}#{r['index']}" for r in rows
                if r["cone"]["simple"] != r["cone"]["simplicial_tangent_cone"]]
    p5 = {"verdict": "PASS" if (simplicial == 183 and unimodular < 92) else "FAIL",
          "predicted_simplicial": 183, "measured_simplicial": simplicial,
          "measured_simple": simple, "simple_vs_simplicial_mismatch": mismatch,
          "measured_unimodular": unimodular, "unimodular_bound": 92, **scope}

    groups: dict[tuple[str, str], set[int]] = {}
    for r in rows:
        groups.setdefault((r["system"], r["support_invariant"]), set()).add(
            r["stabilizer"]["stab_dim"])
    split = [f"{k[0]}:{sorted(v)}" for k, v in groups.items() if len(v) > 1]
    p7 = {"verdict": "PASS" if not split else "FAIL",
          "classes": len(groups), "classes_with_split_stabilizer": len(split),
          "examples": split[:10], **scope}

    contained = [f"{r['system']}#{r['index']}" for r in rows
                 if r["stabilizer"]["stab_dim"] != r["stabilizer"]["commutant_stab_dim"]]
    matrix = determination_matrix(rows)
    mult_bad = [f"{r['system']}#{r['index']}" for r in rows
                if not r["stabilizer"]["commutant_dim_matches_multiplicities"]]
    return {
        "P1_rho_diagonal_count": p1,
        "P2_slater_vertices_are_grassmannian": p2,
        "P3_stabilizer_never_trivial": p3,
        "P4_toral_stabilizer_is_combinatorial": p4,
        "P5_cone_types": p5,
        "P7_stabilizer_constant_on_support_classes": p7,
        "internal_checks": {
            "stab_contained_in_commutant_violations": contained,
            "commutant_dim_vs_multiplicity_violations": mult_bad,
        },
        "determination_matrix": matrix,
    }


# Which presentation each field belongs to, so the matrix reads as a statement
# about presentations and not just about columns.
SOURCES = {
    "support_invariant": ("combinatorial", lambda r: r["support_invariant"]),
    "support_size": ("combinatorial", lambda r: r["support_size"]),
    "verdict_class": ("combinatorial", lambda r: r["classified"]),
    "spectrum": ("spectral", lambda r: tuple(r["integer_form"])),
    "spectrum_multiplicities": ("spectral",
                                lambda r: tuple(r["stabilizer"]["spectrum_multiplicities"])),
    "denominator": ("spectral", lambda r: r["denominator"]),
    "cone_simplicial": ("lattice", lambda r: r["cone"]["simplicial_tangent_cone"]),
    "cone_lattice_index": ("lattice", lambda r: r["cone"]["lattice_index"]),
    "cone_invariant_factors": ("lattice", lambda r: tuple(r["cone"]["invariant_factors"])),
    "active_facets": ("lattice", lambda r: r["cone"]["active_facets"]),
}

TARGETS = {
    "stab_dim": ("symplectic", lambda r: r["stabilizer"]["stab_dim"]),
    "toral_stab_dim": ("symplectic", lambda r: r["stabilizer"]["toral_stab_dim"]),
    "commutant_dim": ("symplectic", lambda r: r["stabilizer"]["commutant_dim"]),
    "commutant_orbit_dim": ("symplectic", lambda r: r["stabilizer"]["commutant_orbit_dim"]),
    "cone_lattice_index": ("lattice", lambda r: r["cone"]["lattice_index"]),
    "cone_simplicial": ("lattice", lambda r: r["cone"]["simplicial_tangent_cone"]),
}


def determination_matrix(rows: list[dict]) -> dict:
    """For each (A, B): does the value of A pin the value of B?

    Reported WITHIN each system, because `stab_dim` lives in `u(d)` and cannot be
    compared across ranks, and then aggregated. `determined` counts the classes
    of A on which B is constant; `rows_determined` weights those classes by size,
    which is the honest figure when one huge class carries the count.
    """
    out: dict[str, dict] = {}
    systems = sorted({r["system"] for r in rows})
    for a_name, (a_kind, a_get) in SOURCES.items():
        for b_name, (b_kind, b_get) in TARGETS.items():
            if a_name == b_name:
                continue
            classes = split = rows_ok = 0
            examples: list[str] = []
            for system in systems:
                buckets: dict[object, set] = {}
                sizes: dict[object, int] = {}
                for r in rows:
                    if r["system"] != system:
                        continue
                    key = a_get(r)
                    buckets.setdefault(key, set()).add(b_get(r))
                    sizes[key] = sizes.get(key, 0) + 1
                for key, values in buckets.items():
                    classes += 1
                    if len(values) == 1:
                        rows_ok += sizes[key]
                    else:
                        split += 1
                        if len(examples) < 5:
                            examples.append(
                                f"{system} {a_name}={key!r} -> "
                                f"{sorted(values, key=repr)}")
            out[f"{a_name}->{b_name}"] = {
                "from": a_kind,
                "to": b_kind,
                # A source whose classes are all singletons determines everything
                # for free and carries no information. `spectrum` is exactly that
                # here: all 799 vertices have distinct spectra.
                "vacuous": classes == len(rows),
                "classes": classes,
                "classes_determined": classes - split,
                "rows_determined": rows_ok,
                "rows": len(rows),
                "determines": split == 0,
                "examples_of_failure": examples,
            }
    return out


if __name__ == "__main__":
    raise SystemExit(main())

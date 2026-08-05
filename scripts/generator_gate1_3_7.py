#!/usr/bin/env python3
"""Generator ladder gate 1: exhaustive (3,7) by the inner/outer closure.

Scores `docs/prereg_generator_gate1_3_7.md`. The route, from
`docs/fixed_n3_generator_methodology.md`: if every H-row is proved valid, the
H-polyhedron is vertex-enumerated without omission, and every enumerated
H-vertex has an exact attainable state, then P subset H by validity and
H subset P by convexity, so P = H.

Nothing here is read off a manifest and believed. The four Ressayre
certificates are replayed from their recorded (h, z) and evaluation point,
with the weights, level sets, roots and tangent matrix rebuilt from scratch
and the determinant recomputed by exact Fraction elimination. Each attaining
state's 1-RDM spectrum is recomputed from its support and amplitudes rather
than read from the record's own integer form.

Predictions G6 and G7 are the controls: they check that the four rows are
doing the cutting, individually and together, so a PASS cannot come from the
ordered slice alone.
"""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

from gpc_census import exact_polytope as XP

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
MANIFEST = DATA / "stage1_ressayre_3_7.json"
STATES = DATA / "states.jsonl"
VERTICES = DATA / "vertices" / "vertices_3_7.json"
DEFAULT_OUT = DATA / "generator_gate1_3_7.json"

N, D = 3, 7


# ----------------------------------------------------------- exact replay


def weights(n: int, d: int):
    return [
        tuple(1 if i in occ else 0 for i in range(d))
        for occ in itertools.combinations(range(d), n)
    ]


def root_action(det, p, q):
    """a_p^dag a_q on an ordered occupation tuple, sign by inversion count."""
    if q not in det or (p in det and p != q):
        return None
    lst = [p if x == q else x for x in det]
    inv = sum(
        1
        for i in range(len(lst))
        for j in range(i + 1, len(lst))
        if lst[i] > lst[j]
    )
    return tuple(sorted(lst)), (-1) ** inv


def exact_det(matrix) -> int:
    """Gaussian elimination over Fraction; a different algorithm from Bareiss."""
    rows = [[F(v) for v in row] for row in matrix]
    size = len(rows)
    result = F(1)
    for col in range(size):
        piv = next((r for r in range(col, size) if rows[r][col]), None)
        if piv is None:
            return 0
        if piv != col:
            rows[col], rows[piv] = rows[piv], rows[col]
            result = -result
        result *= rows[col][col]
        inv = 1 / rows[col][col]
        rows[col] = [v * inv for v in rows[col]]
        for r in range(col + 1, size):
            f = rows[r][col]
            if f:
                rows[r] = [a - f * b for a, b in zip(rows[r], rows[col])]
    assert result.denominator == 1
    return int(result)


def replay_certificate(cert) -> dict:
    """Rebuild every object from (h, z) and recompute the determinant."""
    h, z = list(cert["h"]), cert["z"]
    ws = weights(N, D)
    on, below = [], []
    for i, w in enumerate(ws):
        value = sum(a * b for a, b in zip(w, h))
        if value == z:
            on.append(i)
        elif value < z:
            below.append(i)
    roots = [
        (i, j) for i in range(D) for j in range(i + 1, D) if h[j] - h[i] < 0
    ]
    trace_ok = len(below) == len(roots)
    result = {
        "h": h,
        "z": z,
        "on_matches_record": on == list(cert["on_indices"]),
        "below_matches_record": below == list(cert["below_indices"]),
        "roots_match_record": [list(r) for r in roots]
        == [list(r) for r in cert["negative_roots"]],
        "trace_condition": trace_ok,
        "below": len(below),
        "roots": len(roots),
    }
    if not trace_ok:
        result["determinant"] = None
        result["certified"] = False
        return result
    ev = list(cert["determinant_evaluation"])
    occ = [tuple(i for i, x in enumerate(ws[k]) if x) for k in range(len(ws))]
    below_row = {below[r]: r for r in range(len(below))}
    matrix = [[0] * len(roots) for _ in below]
    for col, (i, j) in enumerate(roots):
        for var, wi in enumerate(on):
            act = root_action(occ[wi], i, j)
            if act is None:
                continue
            target = tuple(1 if t in act[0] else 0 for t in range(D))
            ti = ws.index(target)
            if ti in below_row:
                matrix[below_row[ti]][col] += act[1] * ev[var]
    value = exact_det(matrix)
    result["determinant"] = value
    result["determinant_matches_record"] = value == cert["determinant_value"]
    result["certified"] = value != 0
    return result


# ------------------------------------------------------- polytope assembly


def ordered_slice_rows():
    """Structural rows of the ordered trace-N slice, as (coeffs, rhs) <= form."""
    rows = []
    for i in range(D):
        row = [F(0)] * D
        row[i] = F(1)
        rows.append((list(row), F(1)))          # l_i <= 1
        rows.append(([-v for v in row], F(0)))  # -l_i <= 0
    for i in range(D - 1):
        row = [F(0)] * D
        row[i] = F(-1)
        row[i + 1] = F(1)
        rows.append((row, F(0)))                # l_{i+1} - l_i <= 0
    return rows


def ressayre_rows(manifest) -> list[tuple[list[F], F]]:
    """The retained rows as `coeffs . lambda <= rhs` on the ordered slice."""
    rows = []
    for row in manifest["reduction"]["retained"]:
        constraint = row["constraint"]
        rows.append(([F(c) for c in constraint["coeffs"]], F(constraint["rhs"])))
    return rows


def enumerate_vertices(extra_rows) -> list[tuple[F, ...]]:
    rows = [(list(c), r) for c, r in ordered_slice_rows() + list(extra_rows)]
    eqs = [([F(1)] * D, F(N))]
    verts = XP.brute_force_vertices(rows, eqs, D)
    return sorted(tuple(v) for v in verts)


# ---------------------------------------------------------- attainability


def spectrum_of_state(record) -> tuple[F, ...] | None:
    """Recompute the 1-RDM spectrum from the record's own support.

    Deliberately does not read `integer_form`: that field is what we are
    checking, so using it would make G3 circular.
    """
    support = record.get("support")
    if not support:
        return None
    dets, amps = [], []
    for entry in support:
        det, re_part, im_part = entry
        dets.append(tuple(sorted(det)))
        amps.append(complex(re_part, im_part))
    import numpy as np

    gamma = np.zeros((D, D), dtype=complex)
    index = {d: i for i, d in enumerate(dets)}
    for a, det in enumerate(dets):
        for p in range(D):
            for q in range(D):
                act = root_action(det, p, q)
                if act is None:
                    continue
                target, sgn = act
                if target in index:
                    gamma[p, q] += np.conj(amps[index[target]]) * sgn * amps[a]
    ev = sorted(np.linalg.eigvalsh(gamma), reverse=True)
    return tuple(ev)


def load_states():
    out = {}
    for line in STATES.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec["system"] == f"({N},{D})":
            out[rec["index"]] = rec
    return out


# ------------------------------------------------------------------ build


def build() -> dict:
    manifest = json.loads(MANIFEST.read_text())
    reference = json.loads(VERTICES.read_text())
    ref_set = sorted(
        tuple(F(x) for x in v["spectrum"]) for v in reference
    )
    rows = ressayre_rows(manifest)

    # G4: replay every certificate
    replays = []
    for row, retained in zip(rows, manifest["reduction"]["retained"]):
        cert = retained.get("ressayre_certificate") or retained.get("certificate")
        replays.append(replay_certificate(cert))
    g4 = all(r["certified"] and r["trace_condition"] for r in replays)

    # G1, G2
    verts = enumerate_vertices(rows)
    g1 = len(verts) == 10
    g2 = verts == ref_set

    # G6, G7: the controls
    slice_only = enumerate_vertices([])
    g6 = len(slice_only) > len(verts)
    drops = []
    for k in range(len(rows)):
        subset = rows[:k] + rows[k + 1 :]
        count = len(enumerate_vertices(subset))
        drops.append({"dropped_row": k, "vertices": count, "load_bearing": count > len(verts)})
    g7 = all(d["load_bearing"] for d in drops)

    # G3
    states = load_states()
    attain = []
    for i, v in enumerate(ref_set):
        rec = next(
            (r for r in states.values()
             if tuple(F(x, r["denominator"]) for x in r["integer_form"]) == v),
            None,
        )
        entry = {"vertex": [str(x) for x in v]}
        if rec is None:
            entry.update({"state": None, "ok": False})
        else:
            spec = spectrum_of_state(rec)
            match = spec is not None and all(
                abs(float(a) - b) < 1e-9 for a, b in zip(v, spec)
            )
            entry.update({
                "index": rec["index"],
                "tierB": rec.get("tierB"),
                "recomputed_spectrum_matches": match,
                "derived": rec.get("classified") is not None,
                "ok": rec.get("tierB") in {"EXACT", "EXACT-CONSTR"} and match,
            })
        attain.append(entry)
    g3 = len(attain) == len(ref_set) and all(a["ok"] for a in attain)

    g5 = all([g1, g2, g3, g4])
    return {
        "schema_version": 1,
        "system": f"({N},{D})",
        "prereg": "docs/prereg_generator_gate1_3_7.md",
        "route": "independent exact inner/outer closure",
        "scorecard": {
            "G1_ten_vertices": "PASS" if g1 else "FAIL",
            "G2_matches_reference": "PASS" if g2 else "FAIL",
            "G3_every_vertex_attained": "PASS" if g3 else "FAIL",
            "G4_every_row_certified_on_replay": "PASS" if g4 else "FAIL",
            "G5_closure_exact": "PASS" if g5 else "FAIL",
            "G6_rows_load_bearing_together": "PASS" if g6 else "FAIL",
            "G7_each_row_load_bearing": "PASS" if g7 else "FAIL",
        },
        "vertex_count": len(verts),
        "reference_vertex_count": len(ref_set),
        "ordered_slice_only_vertex_count": len(slice_only),
        "single_row_drops": drops,
        "certificate_replays": replays,
        "attainability": attain,
        "conclusion": (
            "P = H for (3,7): every H-row proved valid on replay, the "
            "H-polyhedron vertex-enumerated exactly, and every vertex carries "
            "an exact attaining state whose spectrum was recomputed here."
            if g5 else
            "closure NOT established; see the scorecard"
        ),
        "nonclaims": [
            "This is the closure route, not the enumeration route. Candidate "
            "exhaustiveness at (3,7) is bypassed, not proved, and the "
            "manifest's system_completeness field is unchanged.",
            "Nothing here bears on (3,8) or any higher rank.",
            "The ten (3,7) states are not ten independent observations; "
            "several are transported or padded from (3,6).",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    payload = build()
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for name, verdict in payload["scorecard"].items():
        print(f"  {verdict:4s} {name}")
    print(f"vertices {payload['vertex_count']} "
          f"(reference {payload['reference_vertex_count']}, "
          f"ordered slice alone {payload['ordered_slice_only_vertex_count']})")
    print(f"wrote {args.out}")
    return 0 if payload["scorecard"]["G5_closure_exact"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

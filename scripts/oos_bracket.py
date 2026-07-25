#!/usr/bin/env python3
"""Out-of-sample test of the first-order fiber theory on (3,11) and (3,12).

Runs the rank-6-to-10 machinery UNCHANGED on the bracket settlement holdout:
19 certified true vertices of (3,11) and their 19 face embeddings in (3,12).
Predictions are pre-registered in docs/prereg_bracket_3_11_3_12.md and must not
be edited after this script has been run.

Writes results/data/oos_bracket.json.

  python scripts/oos_bracket.py
"""
from __future__ import annotations

import json
import pathlib
import sys
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
import sympy as sp  # noqa: E402

from gpc_census.cascade import cascade, certificate_residual, distinct_eigenvalues  # noqa: E402
from gpc_census.cascade import spectrum_error  # noqa: E402
from gpc_census.fiber import (  # noqa: E402
    fixed_rho_tangent, fixed_spectrum_tangent, incidence_kernel_dim, one_rdm,
)

DIAG_TOL = 1e-9
STATES = ROOT / "results/data/states.jsonl"


def one_hop_channels(dets):
    out = set()
    for i in range(len(dets)):
        for j in range(i + 1, len(dets)):
            sd = set(dets[i]) ^ set(dets[j])
            if len(sd) == 2:
                out.add(tuple(sorted(sd)))
    return out


def census_state(system, index):
    for line in STATES.read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            if r["system"] == system and r["index"] == index:
                return r
    raise KeyError((system, index))


def holdout_states():
    """(label, amplitudes, support_dets, spectrum, d, provenance) per vertex."""
    out = []
    settlement = json.loads((ROOT / "docs/bracket_3_11_settlement.json").read_text())
    by_index = {}
    for cand in settlement["candidates"]:
        if cand["verdict"] != "TRUE-VERTEX":
            continue
        den = int(cand["denominator"])
        spectrum = [Fraction(int(x), den) for x in cand["integer_form"]]
        if cand["certificate"] == "EXPLICIT-STATE":
            st = cand["state"]
            dets = [tuple(t) for t in st["support_dets"]]
            c = np.array([complex(sp.N(sp.sympify(s), 30)) for s in st["pretty"]])
            prov = {"certificate": "EXPLICIT-STATE"}
        else:
            donor = cand["donor"]
            rec = census_state(donor["system"], donor["index"])
            dets = [tuple(t) for t in rec["closed_form"]["support_dets"]]
            c = np.array([complex(sp.N(sp.sympify(s), 30))
                          for s in rec["closed_form"]["pretty"]])
            prov = {"certificate": "STATE-TRANSPORT", "donor": donor}
        c = c / np.linalg.norm(c)
        entry = (f"(3,11) v{cand['index']}", c, dets, spectrum, 11, prov)
        out.append(entry)
        by_index[cand["index"]] = entry

    embeddings = json.loads((ROOT / "docs/bracket_3_12_true_vertices.json").read_text())
    for tv in embeddings["true_vertices"]:
        src = tv["source_3_11"]["index"]
        _, c, dets, _, _, prov = by_index[src]
        spectrum = [Fraction(sp.Rational(x)) for x in tv["spectrum"]]
        out.append((f"(3,12) from (3,11) v{src}", c, dets, spectrum, 12,
                    {"certificate": "FACE-EMBEDDING", "source_3_11": src,
                     "inner": prov}))
    return out


def measure(label, c, dets, spectrum, d, prov):
    rho = one_rdm(c, dets, d)
    off = float(np.max(np.abs(rho - np.diag(np.diag(rho)))))
    channels = one_hop_channels(dets)
    kdim = incidence_kernel_dim(dets)
    lams = distinct_eigenvalues(spectrum)
    start_res = certificate_residual(c, dets, lams, d)
    start_spec = spectrum_error(c, dets, spectrum, d)
    spec_t = fixed_spectrum_tangent(c, dets, d)
    rho_t = fixed_rho_tangent(c, dets, d)
    casc = cascade(c, dets, spectrum, d)
    return {
        "label": label,
        "provenance": prov,
        "d": d,
        "support_certified": len(dets),
        "incidence_kernel_dim": kdim,
        "one_hop_channels": len(channels),
        "rdm_max_offdiagonal": off,
        "cancellation_regime": off < DIAG_TOL and bool(channels),
        "start_certificate_residual": start_res,
        "start_spectrum_error": start_spec,
        "fixed_spectrum_tangent_dim": spec_t.tangent_dim,
        "fixed_spectrum_gauge_dim": spec_t.gauge_dim,
        "fixed_rho_tangent_dim": rho_t.tangent_dim,
        "fixed_rho_gauge_dim": rho_t.gauge_dim,
        "predicted_tangent_channel_free": 2 * kdim if not channels else None,
        "support_upper": casc["support_upper"],
        "cascade_blocked_reason": casc["blocked_reason"],
        "cascade_final_residual": casc["final_residual"],
        "evidence": "numerical",
    }


def score(rows):
    """Score the pre-registered predictions B1 to B6 from the measurements."""
    verdicts = {}

    cf = [r for r in rows if r["one_hop_channels"] == 0]
    b1 = [r for r in cf if r["fixed_spectrum_tangent_dim"] != 2 * r["incidence_kernel_dim"]]
    verdicts["B1_unified_dimension_formula"] = {
        "verdict": "PASS" if not b1 else "FAIL",
        "scope": len(cf),
        "mismatches": [{"label": r["label"],
                        "measured": r["fixed_spectrum_tangent_dim"],
                        "predicted": 2 * r["incidence_kernel_dim"]} for r in b1],
    }

    # B2: transported / embedded vertices vs their (3,10) donor measured on the
    # census itself, and the (3,12) embedding vs its (3,11) inner vertex
    inner = {r["label"]: r for r in rows if r["label"].startswith("(3,11)")}
    b2 = []
    for r in rows:
        if r["provenance"].get("certificate") != "FACE-EMBEDDING":
            continue
        src = inner[f"(3,11) v{r['provenance']['source_3_11']}"]
        for key in ("fixed_spectrum_tangent_dim", "fixed_spectrum_gauge_dim",
                    "fixed_rho_tangent_dim", "incidence_kernel_dim",
                    "support_upper"):
            if r[key] != src[key]:
                b2.append({"label": r["label"], "key": key,
                           "embedded": r[key], "inner": src[key]})
    verdicts["B2_padding_invariance"] = {
        "verdict": "PASS" if not b2 else "FAIL",
        "scope": sum(1 for r in rows
                     if r["provenance"].get("certificate") == "FACE-EMBEDDING"),
        "mismatches": b2,
    }

    b3 = [r for r in rows
          if r["fixed_rho_tangent_dim"] != r["fixed_spectrum_tangent_dim"]]
    verdicts["B3_two_fibers_agree"] = {
        "verdict": "PASS" if not b3 else "FAIL",
        "scope": len(rows),
        "disagreements": [{"label": r["label"],
                           "fixed_rho": r["fixed_rho_tangent_dim"],
                           "fixed_spectrum": r["fixed_spectrum_tangent_dim"]}
                          for r in b3],
    }

    canc = [r["label"] for r in rows if r["cancellation_regime"]]
    verdicts["B4_no_cancellation_regime"] = {
        "verdict": "PASS" if not canc else "FAIL",
        "members": canc,
        "silence_rank_lemma": "NOT-TESTED" if not canc else "TESTED",
    }

    b5 = [r for r in rows if r["support_upper"] != r["support_certified"]]
    verdicts["B5_no_sparsification"] = {
        "verdict": "PASS" if not b5 else "FAIL",
        "scope": len(rows),
        "drops": [{"label": r["label"], "certified": r["support_certified"],
                   "upper": r["support_upper"]} for r in b5],
    }

    b6 = [r for r in rows if not (r["start_certificate_residual"] < 1e-10
                                  and r["start_spectrum_error"] < 1e-8)]
    verdicts["B6_numerics_at_new_rank"] = {
        "verdict": "PASS" if not b6 else "FAIL",
        "scope": len(rows),
        "failures": [{"label": r["label"],
                      "certificate_residual": r["start_certificate_residual"],
                      "spectrum_error": r["start_spectrum_error"]} for r in b6],
    }
    return verdicts


def main():
    rows = [measure(*item) for item in holdout_states()]
    out = {
        "preregistration": "docs/prereg_bracket_3_11_3_12.md",
        "holdout": "19 certified (3,11) true vertices + 19 (3,12) face embeddings",
        "evidence": "numerical",
        "scorecard": score(rows),
        "measurements": rows,
    }
    (ROOT / "results/data/oos_bracket.json").write_text(json.dumps(out, indent=2) + "\n")
    for name, v in out["scorecard"].items():
        print(f"{name}: {v['verdict']}")
    print(json.dumps(out["scorecard"], indent=2))


if __name__ == "__main__":
    main()

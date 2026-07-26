#!/usr/bin/env python3
"""Natural-orbital representatives for every interference record.

WHY. 154 of the 799 shipped records fail the diagonality gate: they attain the
spectrum in a non-natural-orbital basis, so their support bounds s_Q^free and
not s_Q^NO (results/data/attainer_gate_audit.json). This produces, for every
interference record, a representative whose 1-RDM IS exactly diag(lambda), so
that a natural-orbital support bound exists for every vertex.

HOW. Diagonalize the record's 1-RDM as rho = U diag(w) U^H with w sorted
descending to match lambda, then push the state through the N-th compound of
U^T, which is the orbital change of basis on determinant amplitudes. The
resulting state has 1-RDM exactly diag(lambda) by construction. Note the
convention: the correct matrix is the TRANSPOSE U^T and not the conjugate
transpose. The two coincide for real rotations, which is why an incorrect
convention still appears to work on the real records and fails only on the
complex ones.

THE COST, stated plainly. Rotating into natural orbitals destroys sparsity.
Support grows on 141 of the 156 records and is unchanged on 15, never
shrinking; the mean goes from 6.8 to 11.7 and the maximum from 14 to 57. And
the rotated amplitudes are NUMERICAL: the shipped library's exact symbolic
closed forms do not survive the rotation, and re-recognizing them is a
research task, not a mechanical one.

NON-UNIQUENESS, which bounds what these numbers mean. When lambda is
degenerate the natural-orbital basis is not unique: any unitary acting within
an eigenspace of rho preserves diagonality. The support recorded here comes
from numpy eigh's arbitrary choice of eigenbasis, so it is an UPPER bound on
s_Q^NO for that choice and is NOT minimized over the eigenspace freedom. The
two records that already satisfy the gate (v89 and v103) are passed through
unrotated for exactly this reason: rotating them by an arbitrary basis of
their own degenerate eigenspaces scrambles them for nothing, which is how an
earlier version of this ledger reported v103 at support 79 instead of 14.

So this ledger does NOT replace states.jsonl. It is an additional artifact,
and the two answer different questions:

    states.jsonl                  EXACT closed forms, sparse supports,
                                  spec(rho) = lambda; for 154 interference
                                  records rho is not diagonal, so the support
                                  bounds s_Q^free.
    states_natural_orbital.jsonl  NUMERICAL amplitudes, larger supports,
                                  rho = diag(lambda) exactly; the support
                                  bounds s_Q^NO.

Keeping both is also what the originating handoff asked for: "Keep the current
record as a documented spectrum-attainer artifact, do not silently overwrite
it."

Usage:
  python scripts/natural_orbital_ledger.py           # write the ledger
  python scripts/natural_orbital_ledger.py --check
  python scripts/natural_orbital_ledger.py --summary
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gpc_census.fiber import amplitudes_from_record, one_rdm  # noqa: E402
from gpc_census.gauge import compound_matrix, two_rdm_spectrum  # noqa: E402

DATA = ROOT / "results" / "data"
OUT = DATA / "states_natural_orbital.jsonl"
SUMMARY = DATA / "natural_orbital_summary.json"
AMP_TOL = 1e-11


def to_natural_orbitals(rec):
    """Rotate a record into its natural-orbital basis. Returns a new record."""
    amps, dets = amplitudes_from_record(rec)
    dets = [tuple(int(o) for o in t) for t in dets]
    lam = np.array([v / rec["denominator"] for v in rec["integer_form"]], dtype=float)
    d = len(lam)
    n = len(dets[0])
    c = np.asarray(amps, dtype=complex)

    rho = one_rdm(c, dets, d)
    before_off = float(np.max(np.abs(rho - np.diag(np.diag(rho)))))
    before_diag_err = float(np.max(np.abs(np.real(np.diag(rho)) - lam)))

    # A record that already IS a natural-orbital representative must be left
    # alone. Rotating it is not a no-op: when lambda is degenerate, eigh
    # returns an arbitrary orthonormal basis of each eigenspace rather than
    # the identity, so "rotating" an already-diagonal rho scrambles the state
    # inside its degenerate blocks and destroys sparsity for nothing. This is
    # how v103 was blown from support 14 to 79 in the first version of this
    # ledger, despite already passing the gate.
    already_natural = before_off < 1e-9 and before_diag_err < 1e-9
    full = list(combinations(range(d), n))
    if already_natural:
        index = {t: i for i, t in enumerate(full)}
        rotated = np.zeros(len(full), dtype=complex)
        for amp, t in zip(c, dets):
            rotated[index[t]] = amp
    else:
        w, u = np.linalg.eigh(rho)
        u = u[:, np.argsort(w)[::-1]]  # eigenvectors ordered to match lambda
        rotated = compound_matrix(u.T, full, dets) @ c

    keep = [i for i, a in enumerate(rotated) if abs(a) > AMP_TOL]
    support = [list(full[i]) for i in keep]
    kept = rotated[keep]

    rho_new = one_rdm(rotated, full, d)
    off = float(np.max(np.abs(rho_new - np.diag(np.diag(rho_new)))))
    diag_err = float(np.max(np.abs(np.real(np.diag(rho_new)) - lam)))
    truncation = float(np.linalg.norm(rotated) - np.linalg.norm(kept))

    # basis-free fingerprint: the rotation must not have changed the state
    same_state = float(
        np.max(np.abs(two_rdm_spectrum(c, dets, d) - two_rdm_spectrum(kept, support, d)))
    )

    return {
        "system": rec["system"],
        "index": rec["index"],
        "classified": rec["classified"],
        "integer_form": rec["integer_form"],
        "denominator": rec["denominator"],
        "support_dets": support,
        "amplitudes_real": [float(np.real(a)) for a in kept],
        "amplitudes_imag": [float(np.imag(a)) for a in kept],
        "support": len(support),
        "support_in_states_jsonl": len(dets),
        "max_offdiagonal_before": before_off,
        "already_natural_orbital": bool(already_natural),
        "max_offdiagonal_after": off,
        "diagonal_error_after": diag_err,
        "truncation_norm_loss": truncation,
        "two_rdm_spectrum_deviation": same_state,
        "bounds": "s_Q^NO",
        "evidence": "NUMERICAL: rotation is floating point, no exact closed form",
        "provenance": "natural-orbital rotation of the states.jsonl record",
    }


def build():
    records = [
        json.loads(ln)
        for ln in (DATA / "states.jsonl").read_text().splitlines()
        if ln.strip()
    ]
    rows = [to_natural_orbitals(r) for r in records if r["classified"] == "INTERFERENCE"]
    grew = sum(1 for r in rows if r["support"] > r["support_in_states_jsonl"])
    same = sum(1 for r in rows if r["support"] == r["support_in_states_jsonl"])
    summary = {
        "generated_by": "scripts/natural_orbital_ledger.py",
        "ledger": "results/data/states_natural_orbital.jsonl",
        "records": len(rows),
        "all_diagonal_after": bool(
            all(r["max_offdiagonal_after"] < 1e-9 for r in rows)
        ),
        "all_diagonal_equals_lambda": bool(
            all(r["diagonal_error_after"] < 1e-9 for r in rows)
        ),
        "worst_offdiagonal_after": max(r["max_offdiagonal_after"] for r in rows),
        "worst_diagonal_error_after": max(r["diagonal_error_after"] for r in rows),
        "worst_truncation_norm_loss": max(abs(r["truncation_norm_loss"]) for r in rows),
        "worst_two_rdm_deviation": max(r["two_rdm_spectrum_deviation"] for r in rows),
        "already_natural_orbital": sum(1 for r in rows if r["already_natural_orbital"]),
        "support_grew": grew,
        "support_unchanged": same,
        "support_shrank": len(rows) - grew - same,
        "support_before": {
            "min": min(r["support_in_states_jsonl"] for r in rows),
            "max": max(r["support_in_states_jsonl"] for r in rows),
            "mean": sum(r["support_in_states_jsonl"] for r in rows) / len(rows),
        },
        "support_after": {
            "min": min(r["support"] for r in rows),
            "max": max(r["support"] for r in rows),
            "mean": sum(r["support"] for r in rows) / len(rows),
        },
        "evidence": "NUMERICAL",
        "reading": (
            "Every interference record now has a natural-orbital "
            "representative whose 1-RDM is exactly diag(lambda), so a "
            "s_Q^NO bound exists for every vertex. The cost is sparsity and "
            "exactness: support grows on 142 of 156 records and never shrinks, "
            "and the rotated amplitudes are numerical because the exact closed "
            "forms do not survive the rotation. states.jsonl is therefore NOT "
            "replaced; it remains the exact, sparse, spectrum-attaining "
            "library, and the two ledgers bound different quantities. "
            "CAVEAT on s_Q^NO: when lambda is degenerate the natural-orbital "
            "basis is NOT unique, since any unitary within an eigenspace of "
            "rho preserves diagonality. The support here is the one produced "
            "by numpy eigh's arbitrary choice of eigenbasis, so it is an "
            "UPPER bound on s_Q^NO that depends on that choice and is not "
            "minimized over the eigenspace freedom. Records that already "
            "satisfied the gate are passed through unrotated rather than "
            "scrambled by an arbitrary basis of their own eigenspaces."
        ),
    }
    return rows, summary


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--summary", action="store_true")
    a = ap.parse_args()

    rows, summary = build()
    if a.summary:
        print(json.dumps(summary, indent=2))
        return 0

    ledger = "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows)
    text = json.dumps(summary, indent=1) + "\n"
    if a.check:
        stale = (
            not OUT.exists()
            or OUT.read_text() != ledger
            or not SUMMARY.exists()
            or SUMMARY.read_text() != text
        )
        if stale:
            print("STALE: natural-orbital ledger", file=sys.stderr)
            return 1
        print("natural-orbital ledger is current")
        return 0

    OUT.write_text(ledger)
    SUMMARY.write_text(text)
    print(f"wrote {OUT} ({summary['records']} records) and {SUMMARY}")
    print(f"  all diagonal after rotation: {summary['all_diagonal_after']} "
          f"(worst {summary['worst_offdiagonal_after']:.2e})")
    print(f"  diagonal equals lambda:      {summary['all_diagonal_equals_lambda']} "
          f"(worst {summary['worst_diagonal_error_after']:.2e})")
    print(f"  same state (2-RDM spectrum): worst deviation "
          f"{summary['worst_two_rdm_deviation']:.2e}")
    print(f"  support grew {summary['support_grew']}, unchanged "
          f"{summary['support_unchanged']}, shrank {summary['support_shrank']}")
    print(f"  support mean {summary['support_before']['mean']:.1f} -> "
          f"{summary['support_after']['mean']:.1f}, "
          f"max {summary['support_before']['max']} -> {summary['support_after']['max']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

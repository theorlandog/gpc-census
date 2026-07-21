"""Certify uncertified vertices by transporting a certified sibling's state
along the census's structural maps.

Trailing-zero padding and frozen-core lift (lambda -> (1, lambda)) transport
STATES, not just verdicts: an eigenvalue-0 orbital is empty (in no determinant),
an eigenvalue-1 orbital is frozen (in every determinant), so stripping them from
a certified state yields the state of the reduced spectrum, and re-adding them
yields the state of any spectrum with the same core (the eigenvalues strictly
between 0 and 1). Two vertices with the same core are transport-equivalent.

For every uncertified vertex that shares a core with a certified donor, this
strips the donor to its core state and re-adds the target's frozen and empty
orbitals, then accepts only if verify_exact certifies the transported state
against the target spectrum (exact characteristic-polynomial identity). Merges
the accepted records into --out (tierB EXACT, provenance in a `transport` field).

    uv run scripts/transport_states.py            # dry run, lists what transports
    uv run scripts/transport_states.py --apply    # merge into results/data/states.jsonl
"""
from __future__ import annotations

import argparse
import json
import pathlib
from fractions import Fraction

DATA = pathlib.Path(__file__).resolve().parents[1] / "results" / "data"
CERT = {"EXACT", "EXACT-CONSTR"}


def parse_sys(s):
    a, b = s.strip("()").split(",")
    return int(a), int(b)


def core_of(spectrum):
    return tuple(sorted(str(f) for f in (Fraction(s) for s in spectrum) if 0 < f < 1))


def core_state(cf):
    """Donor closed form -> (core determinants, amplitudes), frozen orbitals removed."""
    import sympy as sp

    dets = [list(t) for t in cf["support_dets"]]
    amps = [sp.sympify(p) for p in cf["pretty"]]
    used = sorted({o for det in dets for o in det})
    setd = [set(det) for det in dets]
    ones = [o for o in used if all(o in dd for dd in setd)]   # frozen (eigenvalue 1)
    keep = [o for o in used if o not in ones]                 # core orbitals
    relabel = {o: i for i, o in enumerate(keep)}
    cdets = [[relabel[o] for o in det if o in relabel] for det in dets]
    return cdets, amps


def transport(cf, n_t, d_t, spectrum_t):
    """Return (support_dets, amps) for the target if verify_exact certifies, else None."""
    from gpc_census.exactify import verify_exact

    cdets, amps = core_state(cf)
    ones = sum(1 for s in spectrum_t if Fraction(s) == 1)
    tdets = [list(range(ones)) + [o + ones for o in det] for det in cdets]
    spec = [Fraction(s) for s in spectrum_t]
    if verify_exact(n_t, d_t, spec, tdets, amps):
        return tdets, amps
    return None


def ph_transport(cf, n_t, d_t, spectrum_t):
    """Particle-hole transport: the complement bijection T -> T^c carries a state
    for spectrum mu to one for 1 - mu (the 1-RDM maps rho -> I - rho). Only useful
    when the dual system d - N equals N (self-dual, d = 2 N), so the dual vertex is
    in the same census system. Accepted only if verify_exact certifies it."""
    import sympy as sp

    from gpc_census.exactify import verify_exact

    dets = [list(t) for t in cf["support_dets"]]
    amps = [sp.sympify(p) for p in cf["pretty"]]
    full = list(range(d_t))
    cdets = [[o for o in full if o not in det] for det in dets]
    spec = [Fraction(s) for s in spectrum_t]
    if verify_exact(n_t, d_t, spec, cdets, amps):   # complement bijection, moduli carry
        return cdets, amps
    return None


def main() -> int:
    import sympy as sp

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DATA / "states.jsonl"))
    ap.add_argument("--apply", action="store_true", help="merge transports into --out")
    args = ap.parse_args()
    out = pathlib.Path(args.out)

    recs = [json.loads(ln) for ln in out.read_text().splitlines() if ln.strip()]
    byidx = {(r["system"], r["index"]): r for r in recs}

    vspec, vcore = {}, {}
    for r in recs:
        n, d = parse_sys(r["system"])
        v = json.loads((DATA / "vertices" / f"vertices_{n}_{d}.json").read_text())[r["index"]]
        vspec[(r["system"], r["index"])] = (v["spectrum"], v["denominator"], v["integer_form"])
        vcore[(r["system"], r["index"])] = core_of(v["spectrum"])

    def speckey(spectrum):
        return tuple(sorted(str(Fraction(s)) for s in spectrum))

    donors: dict = {}          # core -> certified donor keys (padding / lift)
    spec_donor: dict = {}      # spectrum-key -> certified donor key (for PH duals)
    for key, r in byidx.items():
        if r.get("tierB") in CERT and isinstance(r.get("closed_form"), dict):
            donors.setdefault(vcore[key], []).append(key)
            spec_donor[speckey(vspec[key][0])] = key

    cleared = 0
    for key, r in list(byidx.items()):
        if r.get("tierB") in CERT:
            continue
        n_t, d_t = parse_sys(key[0])
        spectrum_t, den_t, intform = vspec[key]
        # candidate donors: same-core (padding/lift), then PH dual (self-dual system)
        cands = [(dk, "transport", transport) for dk in donors.get(vcore[key], [])]
        if d_t == 2 * n_t:  # self-dual system: the PH dual lives in the same system
            dual = spec_donor.get(speckey([str(1 - Fraction(s)) for s in spectrum_t]))
            if dual is not None:
                cands.append((dual, "ph-dual", ph_transport))
        for dk, via, fn in cands:
            res = fn(byidx[dk]["closed_form"], n_t, d_t, spectrum_t)
            if res is None:
                continue
            tdets, amps = res
            cf = {"den": den_t, "weights": byidx[dk]["closed_form"]["weights"],
                  "pretty": [str(a) for a in amps], "support_dets": tdets}
            support = []
            for det, a in zip(tdets, amps):
                c = complex(sp.N(a, 30))
                support.append([det, abs(c), float(sp.arg(a).evalf())])
            rec = {"system": key[0], "index": key[1], "classified": r["classified"],
                   "integer_form": intform, "denominator": den_t, "status": "OK",
                   "tierB": "EXACT", "support": support, "closed_form": cf,
                   "transport": {"from_system": dk[0], "from_index": dk[1], "via": via}}
            byidx[key] = rec
            cleared += 1
            print(f"  {key[0]} idx {key[1]} <- {dk[0]} idx {dk[1]} ({via})", flush=True)
            break

    print(f"transport certifies {cleared} vertices")
    if args.apply and cleared:
        merged = sorted(byidx.values(), key=lambda r: (r["system"], r["index"]))
        out.write_text("\n".join(json.dumps(r) for r in merged) + "\n")
        n_cert = sum(1 for r in merged if r.get("tierB") in CERT)
        print(f"merged into {out}: {n_cert}/{len(merged)} certified")
        print(f"validate: uv run scripts/validate_states.py {out}")
    elif cleared:
        print("dry run; pass --apply to merge")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

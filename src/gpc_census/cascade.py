"""Sparsification cascade on the FIXED-SPECTRUM fiber.

The certified support of a census state is an UPPER BOUND on the minimal
support s_Q of its vertex, never a minimality claim. This module tries to
lower that bound constructively: it walks the fixed-spectrum fiber towards its
boundary strata by driving one amplitude at a time to zero, keeping the
spectrum exact to numerical precision the whole way.

Method (the one that worked at v103, docs/unified_fiber_dimension.md section 3):

1. Parameterize by t = |c_j|^2 for the currently smallest amplitude j.
2. Pin the orbital-phase gauge by projecting out the gauge directions, so the
   solver cannot drift along the U(1)^d orbit instead of deforming.
3. Solve the SMOOTH certificate prod_i (rho - lambda_i I) = 0, one factor per
   DISTINCT eigenvalue (the minimal polynomial). The eigenvalue residual is
   the wrong objective here: it stalls at degenerate spectra, which is most of
   the census.
4. Continue t down to 0. If the certificate holds at t = 0, drop j from the
   support and repeat on the next smallest amplitude.

Everything this module produces is NUMERICAL. A round that reaches residual
1e-15 is evidence that a sparser exact-spectrum state exists, not a proof;
Smale alpha-theory (Task 3 of the fiber program) is what upgrades it. A round
that floors at 1e-6 or above means blocked ON THIS PATH: the recorded outcome
is the path taken, never a minimality claim about s_Q.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .fiber import gauge_generators

# a round counts as landed only well below this; floors at or above
# FLOOR_TOL are reported as blocked paths
SOLVED_TOL = 1e-10
FLOOR_TOL = 1e-6
# a landed round must also reproduce the eigenvalue MULTIPLICITIES, not just
# the eigenvalue set that the minimal polynomial pins
SPECTRUM_TOL = 1e-8


def rdm_tensor(dets, d):
    """A[p, q] with rho[p, q] = conj(c) @ A[p, q] @ c, canonical det ordering."""
    idx = {t: i for i, t in enumerate(dets)}
    a = np.zeros((d, d, len(dets), len(dets)), dtype=float)
    for ti, t in enumerate(dets):
        for qi, q in enumerate(t):
            rest = [u for u in t if u != q]
            for p in range(d):
                if p in rest:
                    continue
                s = tuple(sorted(rest + [p]))
                si = idx.get(s)
                if si is None:
                    continue
                a[p, q, si, ti] += (-1) ** qi * (-1) ** s.index(p)
    return a


def rdm(c, a):
    return np.einsum("i,pqij,j->pq", np.conj(c), a, c)


def distinct_eigenvalues(spectrum, tol=1e-9):
    """Distinct target eigenvalues, descending: the minimal polynomial roots."""
    out = []
    for x in sorted((float(v) for v in spectrum), reverse=True):
        if not out or abs(x - out[-1]) > tol:
            out.append(x)
    return out


def minpoly(rho, lams):
    """prod_i (rho - lambda_i I) over DISTINCT eigenvalues: 0 iff spec = lambda.

    Using one factor per distinct eigenvalue (not per eigenvalue with
    multiplicity) keeps the certificate smooth at degenerate spectra, where the
    per-eigenvalue residual has a non-smooth minimum and stalls.
    """
    d = rho.shape[0]
    out = np.eye(d, dtype=complex)
    for lam in lams:
        out = out @ (rho - lam * np.eye(d))
    return out


def _certificate(c, a, lams):
    q = minpoly(rdm(c, a), lams)
    iu = np.triu_indices(q.shape[0])
    return np.concatenate([q[iu].real, q[iu].imag, [float(np.vdot(c, c).real) - 1.0]])


def certificate_residual(c, dets, lams, d):
    """max |prod_i (rho - lambda_i I)| together with the norm error."""
    return float(np.max(np.abs(_certificate(np.asarray(c, dtype=complex),
                                            rdm_tensor(dets, d), lams))))


def spectrum_error(c, dets, spectrum, d):
    """max |sorted eig(rho) - sorted lambda|: the MULTIPLICITY check.

    The minimal-polynomial certificate constrains the eigenvalue SET only. On
    most census spectra the trace and rank pin the multiplicities anyway, but
    not on all of them, so every accepted cascade round is gated on this too:
    without it a round could silently land on a different vertex that shares
    the distinct values.
    """
    a = rdm_tensor([tuple(t) for t in dets], d)
    w = np.linalg.eigvalsh(rdm(np.asarray(c, dtype=complex), a))
    tgt = np.sort(np.array([float(v) for v in spectrum]))[::-1]
    return float(np.max(np.abs(np.sort(w)[::-1] - tgt)))


def _solve_at(v0, a, lams, gauge, vref, j, t, gauge_weight=1.0):
    m = a.shape[2]

    def resid(v):
        c = v[:m] + 1j * v[m:]
        return np.concatenate([
            _certificate(c, a, lams),
            gauge_weight * (gauge @ (v - vref)),
            [abs(c[j]) ** 2 - t],
        ])

    sol = least_squares(resid, v0, method="lm", xtol=3e-16, ftol=3e-16, gtol=3e-16,
                        max_nfev=20000)
    c = sol.x[:m] + 1j * sol.x[m:]
    cert = float(np.max(np.abs(_certificate(c, a, lams))))
    return sol.x, cert


def _schedule(t0, steps=8):
    """Geometric descent to zero; the last step is t = 0 exactly."""
    if t0 <= 0:
        return [0.0]
    out = [t0 * (0.35 ** k) for k in range(1, steps + 1)]
    return [x for x in out if x > 1e-14] + [0.0]


def drive_to_zero(c0, dets, lams, d, j, steps=8):
    """Continue the fixed-spectrum fiber driving |c_j|^2 to zero.

    Returns (amplitudes, certificate residual) at the endpoint. A residual
    below SOLVED_TOL means the sparser state was reached (numerically).
    """
    a = rdm_tensor(dets, d)
    gauge = gauge_generators(c0, dets, d)
    m = len(dets)
    v = np.concatenate([np.real(c0), np.imag(c0)])
    vref = v.copy()
    best = (v, float("inf"))
    for t in _schedule(float(abs(c0[j]) ** 2), steps=steps):
        v, cert = _solve_at(v, a, lams, gauge, vref, j, t)
        if t == 0.0:
            best = (v, cert)
    v, cert = best
    return v[:m] + 1j * v[m:], cert


def cascade(amplitudes, support_dets, spectrum, d=None, max_rounds=None,
            rng=None, retries=3):
    """Greedy smallest-first sparsification along the fixed-spectrum fiber.

    Repeatedly drives the smallest surviving amplitude to zero. Where the
    greedy branch stalls, retries from random points of the fiber reached by
    first deforming (a random gauge-orthogonal least-squares re-solve) and then
    cascading, before recording a floor.

    Returns a dict with the starting support size, the best support size
    reached (support_upper), the residual at each landed round, and the floor
    that stopped it. All numbers are numerical, never certified.
    """
    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    if d is None:
        d = max(o for t in dets for o in t) + 1
    lams = distinct_eigenvalues(spectrum)
    rng = np.random.default_rng(0 if rng is None else rng)
    start = len(dets)
    rounds = []
    dropped = []
    floor = None
    drift = False
    limit = start - 1 if max_rounds is None else max_rounds

    for _ in range(limit):
        if len(dets) <= 1:
            break
        order = list(np.argsort(np.abs(c)))
        landed = False
        # greedy smallest-first, then the next few smallest as fallback targets
        for j in order[: min(len(order), 3)]:
            c_out, new_dets, info = _try_drop(c, dets, lams, spectrum, d, int(j))
            if info["ok"]:
                dropped.append(dets[int(j)])
                dets, c = new_dets, c_out
                rounds.append({"dropped_det": list(dropped[-1]),
                               "support_after": len(dets),
                               "certificate_residual": info["certificate_residual"],
                               "spectrum_error": info["spectrum_error"]})
                landed = True
                break
            if info["reason"] == "multiplicity_drift":
                drift = True
            else:
                cert = info["certificate_residual"]
                floor = cert if floor is None else min(floor, cert)
        if landed:
            continue
        # blocked on the greedy branch: deform to other fiber points and retry
        landed, c, dets, extra, floor = _retry_from_fiber_points(
            c, dets, lams, spectrum, d, rng, retries, floor)
        rounds.extend(extra)
        if not landed:
            break

    return {
        "support_start": start,
        "support_upper": len(dets),
        "dropped_dets": [list(t) for t in dropped],
        "rounds": rounds,
        "final_residual": certificate_residual(c, dets, lams, d),
        "final_spectrum_error": spectrum_error(c, dets, spectrum, d),
        "final_support_dets": [list(t) for t in dets],
        "final_amplitudes_real": [float(x) for x in np.real(c)],
        "final_amplitudes_imag": [float(x) for x in np.imag(c)],
        "floor_residual": floor,
        "blocked_on_path": _blocked_reason(len(dets), floor, drift) in (
            "certificate_floor", "multiplicity_drift"),
        "blocked_reason": _blocked_reason(len(dets), floor, drift),
        "multiplicity_drift_seen": drift,
        "method": "greedy-smallest-first minimal-polynomial continuation, "
                  "gauge-pinned, scipy least_squares lm",
        "evidence": "numerical",
        "fiber": "fixed_spectrum",
    }


def _blocked_reason(size, floor, drift):
    """Why the cascade stopped. Never a minimality claim about s_Q.

    certificate_floor: the continuation could not zero the amplitude on this
        path (the floor value says how far it got).
    multiplicity_drift: an amplitude COULD be zeroed, but the resulting 1-RDM
        has the right eigenvalue set with the wrong multiplicities, i.e. the
        path leaves this vertex's fiber for another stratum.
    """
    if size <= 1:
        return "support_exhausted"
    if floor is not None and floor >= FLOOR_TOL:
        return "certificate_floor"
    if drift:
        return "multiplicity_drift"
    return "certificate_floor" if floor is not None else "no_attempt"


def polish(c, dets, lams, d):
    """Re-solve the certificate on a support after a drop.

    Driving |c_j|^2 to zero leaves c_j at roughly the solver tolerance, so
    deleting it and renormalizing perturbs the 1-RDM at first order in that
    leftover. Polishing on the REDUCED support turns the endpoint into an
    honest solution of the smaller system instead of a nearby approximation.
    """
    a = rdm_tensor(dets, d)
    gauge = gauge_generators(c, dets, d)
    m = len(dets)
    v0 = np.concatenate([np.real(c), np.imag(c)])

    def resid(v):
        cc = v[:m] + 1j * v[m:]
        return np.concatenate([_certificate(cc, a, lams), gauge @ (v - v0)])

    sol = least_squares(resid, v0, method="lm", xtol=3e-16, ftol=3e-16,
                        gtol=3e-16, max_nfev=20000)
    cc = sol.x[:m] + 1j * sol.x[m:]
    return cc, float(np.max(np.abs(_certificate(cc, a, lams))))


def _try_drop(c, dets, lams, spectrum, d, j):
    """One drop attempt. Returns (amplitudes, dets, info) with info['ok']."""
    c_new, cert = drive_to_zero(c, dets, lams, d, j)
    if cert >= SOLVED_TOL:
        return None, None, {"ok": False, "reason": "certificate_floor",
                            "certificate_residual": cert}
    keep = [i for i in range(len(dets)) if i != j]
    new_dets = [dets[i] for i in keep]
    c_out = c_new[keep] / np.linalg.norm(c_new[keep])
    c_out, cert = polish(c_out, new_dets, lams, d)
    if cert >= SOLVED_TOL:
        return None, None, {"ok": False, "reason": "certificate_floor",
                            "certificate_residual": cert}
    serr = spectrum_error(c_out, new_dets, spectrum, d)
    if serr >= SPECTRUM_TOL:
        return None, None, {"ok": False, "reason": "multiplicity_drift",
                            "certificate_residual": cert, "spectrum_error": serr}
    return c_out, new_dets, {"ok": True, "certificate_residual": cert,
                             "spectrum_error": serr}


def _retry_from_fiber_points(c, dets, lams, spectrum, d, rng, retries, floor):
    """Deform to a random nearby fiber point, then retry the smallest drop."""
    a = rdm_tensor(dets, d)
    gauge = gauge_generators(c, dets, d)
    m = len(dets)
    v0 = np.concatenate([np.real(c), np.imag(c)])
    for _ in range(retries):
        step = rng.normal(scale=0.05, size=2 * m)

        def resid(v, step=step):
            cc = v[:m] + 1j * v[m:]
            return np.concatenate([
                _certificate(cc, a, lams),
                gauge @ (v - v0),
                0.3 * (v - (v0 + step)),
            ])

        sol = least_squares(resid, v0 + step, method="lm", xtol=1e-14,
                            ftol=1e-14, gtol=1e-14, max_nfev=8000)
        c1 = sol.x[:m] + 1j * sol.x[m:]
        if float(np.max(np.abs(_certificate(c1, a, lams)))) > SOLVED_TOL:
            continue
        for j in list(np.argsort(np.abs(c1)))[:2]:
            c_out, new_dets, info = _try_drop(c1, dets, lams, spectrum, d, int(j))
            if info["ok"]:
                return True, c_out, new_dets, [{
                    "dropped_det": list(dets[int(j)]),
                    "support_after": len(new_dets),
                    "certificate_residual": info["certificate_residual"],
                    "spectrum_error": info["spectrum_error"],
                    "branch": "random-fiber-restart",
                }], floor
            if info["reason"] == "certificate_floor":
                cert = info["certificate_residual"]
                floor = cert if floor is None else min(floor, cert)
    return False, c, dets, [], floor

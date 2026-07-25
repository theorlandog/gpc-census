"""Fixed-1-RDM fiber sampler.

Purpose: decide whether the census-wide `projective_stabilizer_dim >= 2`
observation is driven by EXTREMALITY or merely by SPARSITY.

Control construction.  For a vertex lambda with degeneracy pattern (m_k), set

    mu(t) = (1-t) * lambda  +  t * (N/d) * (1,...,1)

Adding a multiple of the all-ones vector preserves every equality among
coordinates, so mu(t) has EXACTLY the same degeneracy pattern as lambda.  Both
endpoints lie in the moment polytope, so mu(t) does too, and for t>0 it is a
proper convex combination of two distinct points -- hence not a vertex.

We then sample the fixed-rho fiber over mu(t) and compare stabilizer dimensions
at matched degeneracy pattern.  If the vertex keeps more symmetry than interior
points with the same pattern, extremality is the driver.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from gpc_census.levi import determinant_basis, annihilation_matrices, fiber_and_orbit


def _pack(psi):
    return np.concatenate((psi.real, psi.imag))


def _unpack(x, H):
    return x[:H] + 1j * x[H:]


def sample_fiber(target, n, d, n_samples=8, seed=0, tol=1e-12, maxiter=4000):
    """Return states psi (unit norm) with one-RDM ~= diag(target).

    Minimises  f(psi) = || rho(psi) - diag(target) ||_F^2  on the unit sphere
    with an analytic gradient; different random starts land on different
    points of the fiber.
    """
    basis = determinant_basis(n, d)
    H = len(basis)
    ann = annihilation_matrices(n, d, basis)
    T = np.diag(np.asarray(target, dtype=float)).astype(complex)
    rng = np.random.default_rng(seed)

    def fg(x):
        psi = _unpack(x, H)
        nrm = np.linalg.norm(psi)
        if nrm == 0:
            return 1e9, np.zeros_like(x)
        psi = psi / nrm
        B = [a @ psi for a in ann]
        rho = np.array([[np.vdot(B[p], B[q]) for q in range(d)] for p in range(d)])
        G = rho - T
        f = float(np.sum(np.abs(G) ** 2))
        # grad wrt psi*  =  2 sum_p A_p^dag ( sum_q conj(G[p,q]) B_q )
        grad = np.zeros(H, dtype=complex)
        for p in range(d):
            acc = np.zeros(B[0].shape, dtype=complex)
            for q in range(d):
                acc += np.conj(G[p, q]) * B[q]
            grad += ann[p].conj().T @ acc
        grad *= 2.0
        # project off the radial direction (norm constraint) and rescale
        grad = (grad - psi * np.vdot(psi, grad)) / nrm
        return f, _pack(grad)

    out = []
    for s in range(n_samples):
        x0 = _pack((rng.normal(size=H) + 1j * rng.normal(size=H)))
        r = minimize(fg, x0, jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "ftol": tol, "gtol": tol})
        psi = _unpack(r.x, H)
        psi /= np.linalg.norm(psi)
        out.append((float(r.fun), psi))
    return basis, out


def fiber_stabilizers(target, n, d, n_samples=8, seed=0, resid_tol=1e-8):
    """Sample the fiber and report stabilizer dims for the converged samples."""
    basis, samples = sample_fiber(target, n, d, n_samples=n_samples, seed=seed)
    rows = []
    for resid, psi in samples:
        if resid > resid_tol:
            continue
        info = fiber_and_orbit(psi, basis, n, d)
        info["residual"] = resid
        info["support"] = int(np.sum(np.abs(psi) > 1e-8))
        rows.append(info)
    return rows


def interior_point(lam_int, den, n, d, t):
    """mu(t) = (1-t)*lambda + t*(N/d)*ones ; preserves the degeneracy pattern."""
    lam = np.asarray(lam_int, dtype=float) / float(den)
    return (1.0 - t) * lam + t * (float(n) / float(d)) * np.ones(d)

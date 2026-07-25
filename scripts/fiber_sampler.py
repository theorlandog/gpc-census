"""Fixed-SPECTRUM fiber sampler.  (v2 -- bug fix)

BUG IN v1.  The objective was || rho(psi) - diag(target) ||^2, which pins the
NATURAL-ORBITAL BASIS as well as the spectrum.  The census certifies
spec(rho) = lambda, not rho = diag(lambda): 154 of 156 interference states have
non-diagonal 1-RDMs in the computational basis (off-diagonal norm up to 0.314),
exactly the stage-3b off-diagonal-target states.  The two diagonal exceptions
are (3,10) v89 and v103.  Pointing v1 at an interference vertex therefore
sampled the wrong fiber.

FIX.  Objective  f(psi) = sum_i ( eig_i(rho) - lambda_i )^2  with both sorted
descending.  Basis-invariant.  Gradient via the standard eigenvalue perturbation
d eig_i = <v_i| d rho |v_i>, giving

    grad_{psi*} f = 2 sum_p A_p^dag M B_p-ish;  concretely with
    W = V diag(2(w - lambda)) V^dag  (V eigenvectors of rho):
    grad = 2 * sum_{p,q} W[p,q]^* A_p^dag (A_q psi)   ... = dGamma(W^dag) psi.

At eigenvalue crossings (degenerate lambda) the eigenvalues are not smooth;
sorting both sides descending gives a subgradient that works in practice, and
converged samples are verified a posteriori by |spec(rho) - lambda| directly.

The falsification result of v1 (projective_stabilizer_dim is representative-
dependent) STANDS: it used (3,6) DESIGN vertices, whose certified 1-RDMs are
diagonal, so v1's objective was correct there.
"""
from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from gpc_census.levi import determinant_basis, annihilation_matrices, fiber_and_orbit


def _pack(psi):
    return np.concatenate((psi.real, psi.imag))


def _unpack(x, H):
    return x[:H] + 1j * x[H:]


def _rho_and_lowered(psi, ann, d):
    B = [a @ psi for a in ann]
    rho = np.array([[np.vdot(B[p], B[q]) for q in range(d)] for p in range(d)])
    return rho, B


def spectrum_residual(psi, ann, d, target_sorted):
    rho, _ = _rho_and_lowered(psi, ann, d)
    w = np.sort(np.linalg.eigvalsh(rho))[::-1]
    return float(np.max(np.abs(w - target_sorted)))


def sample_spectrum_fiber(target, n, d, n_samples=8, seed=0, tol=1e-14,
                          maxiter=4000):
    """States psi with sorted spec(rho(psi)) ~= sorted target.  Basis-free."""
    basis = determinant_basis(n, d)
    H = len(basis)
    ann = annihilation_matrices(n, d, basis)
    lam = np.sort(np.asarray(target, dtype=float))[::-1]
    rng = np.random.default_rng(seed)

    def fg(x):
        psi = _unpack(x, H)
        nrm = np.linalg.norm(psi)
        if nrm == 0:
            return 1e9, np.zeros_like(x)
        psi = psi / nrm
        rho, B = _rho_and_lowered(psi, ann, d)
        w, V = np.linalg.eigh(rho)
        o = np.argsort(w)[::-1]
        w, V = w[o], V[:, o]
        diff = w - lam
        f = float(np.sum(diff ** 2))
        # W = V diag(2*diff) V^dag ; grad_{psi*} = sum_pq conj(W[p,q]) A_p^dag B_q
        W = (V * (2.0 * diff)) @ V.conj().T
        grad = np.zeros(H, dtype=complex)
        for p in range(d):
            acc = np.zeros(B[0].shape, dtype=complex)
            for q in range(d):
                acc += np.conj(W[p, q]) * B[q]
            grad += ann[p].conj().T @ acc
        grad = (grad - psi * np.vdot(psi, grad)) / nrm
        return f, _pack(grad)

    out = []
    for s in range(n_samples):
        x0 = _pack(rng.normal(size=H) + 1j * rng.normal(size=H))
        r = minimize(fg, x0, jac=True, method="L-BFGS-B",
                     options={"maxiter": maxiter, "ftol": tol, "gtol": tol})
        psi = _unpack(r.x, H)
        psi /= np.linalg.norm(psi)
        out.append((spectrum_residual(psi, ann, d, lam), psi))
    return basis, out


def spectrum_fiber_stabilizers(target, n, d, n_samples=8, seed=0,
                               resid_tol=1e-8):
    basis, samples = sample_spectrum_fiber(target, n, d,
                                           n_samples=n_samples, seed=seed)
    rows = []
    for resid, psi in samples:
        if resid > resid_tol:
            continue
        info = fiber_and_orbit(psi, basis, n, d)
        info["spectrum_residual"] = resid
        info["support"] = int(np.sum(np.abs(psi) > 1e-8))
        rows.append(info)
    return rows


def interior_point(lam_int, den, n, d, t):
    lam = np.asarray(lam_int, dtype=float) / float(den)
    return (1.0 - t) * lam + t * (float(n) / float(d)) * np.ones(d)

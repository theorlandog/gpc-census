"""Structural validation invariants for vertex sets. Zero dependencies.

These checks caught two real errors during the original census; every
pipeline output should pass through them.
"""
from __future__ import annotations

from fractions import Fraction

Vertex = tuple[Fraction, ...]


def check_physical(verts: list[Vertex], n: int) -> list[str]:
    """Trace, ordering, and 0 <= lambda <= 1 for every vertex."""
    errs = []
    for v in verts:
        if sum(v) != n:
            errs.append(f"trace != {n}: {v}")
        if any(v[i] < v[i + 1] for i in range(len(v) - 1)):
            errs.append(f"not descending: {v}")
        if v and (v[0] > 1 or v[-1] < 0):
            errs.append(f"occupation outside [0,1]: {v}")
    return errs


def check_embedding(lower: list[Vertex], higher: list[Vertex]) -> list[str]:
    """Every higher-rank vertex with a trailing zero must be a lower-rank vertex."""
    lowset = set(lower)
    errs = []
    for v in higher:
        if v[-1] == 0 and v[:-1] not in lowset:
            errs.append(f"trailing-zero vertex not in lower rank: {v}")
    return errs


def check_vertex_attainer(amplitudes, support_dets, spectrum, d=None,
                          library=None, claimed_new=False, tol=1e-9):
    """The acceptance gate for a state offered as attaining a vertex.

    Four conditions, the first three unconditional and the fourth required only
    of a state claimed to be a NEW attainer:

    1. the 1-RDM is diagonal,
    2. its diagonal is lambda,
    3. the state has unit norm,
    4. if claimed new, its 2-RDM spectrum DIFFERS from the library
       representative's.

    Conditions 1 and 2 are what make the state a natural-orbital representative,
    so that its support bounds s_Q^NO rather than the free-basis minimum: a
    state passing only spec(rho) = lambda attains the vertex but in some other
    basis, and its support is not comparable with s_I. Condition 4 is the cheap
    basis-free test that separates a new state from the library state written in
    rotated orbitals; the 2-RDM spectrum is invariant under every one-body
    unitary, so equality means no new state has been found. Pass `library` as
    (amplitudes, support_dets) to enable it.

    Returns a list of failure strings, empty when the state is accepted.
    """
    import numpy as np

    from .fiber import one_rdm
    from .gauge import two_rdm_spectrum

    dets = [tuple(int(o) for o in t) for t in support_dets]
    c = np.asarray(amplitudes, dtype=complex)
    lam = [float(x) for x in spectrum]
    if d is None:
        d = len(lam)
    errs = []
    rho = one_rdm(c, dets, d)
    offdiag = float(np.max(np.abs(rho - np.diag(np.diag(rho)))))
    if offdiag > tol:
        errs.append(f"1-RDM not diagonal: max off-diagonal {offdiag:.3e}; this "
                    "state attains the spectrum in another basis, so its "
                    "support bounds s_Q^free and not s_Q^NO")
    diag_err = float(np.max(np.abs(np.real(np.diag(rho)) - np.array(lam))))
    if diag_err > tol:
        errs.append(f"diagonal is not lambda: max deviation {diag_err:.3e}")
    norm_err = abs(float(np.linalg.norm(c)) - 1.0)
    if norm_err > tol:
        errs.append(f"norm != 1: deviation {norm_err:.3e}")
    if claimed_new:
        if library is None:
            errs.append("claimed new but no library representative given to "
                        "compare 2-RDM spectra against")
        else:
            lib_c, lib_dets = library
            dev = float(np.max(np.abs(
                two_rdm_spectrum(c, dets, d)
                - two_rdm_spectrum(np.asarray(lib_c, dtype=complex),
                                   [tuple(int(o) for o in t) for t in lib_dets], d))))
            if dev <= tol:
                errs.append(f"2-RDM spectrum matches the library state to "
                            f"{dev:.3e}: this is the SAME state in different "
                            "orbitals, a gauge statement, not a new attainer")
    return errs


def check_selfdual(verts: list[Vertex], n: int, d: int) -> list[str]:
    """For half filling (2n == d) the vertex set must be complement-closed."""
    if 2 * n != d:
        return []
    vset = set(verts)
    errs = []
    for v in verts:
        dual = tuple(sorted((1 - x for x in v), reverse=True))
        if dual not in vset:
            errs.append(f"complement not a vertex: {v}")
    return errs

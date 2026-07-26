# Natural orbitals are not the sparsest basis: 48 certified counterexamples

Short companion note, NOT part of the polytope paper. It is written for the
quantum-chemistry reading of the same data: the folklore that natural orbitals
give the most compact CI expansion is known not to be a theorem, and the census
data now supply certified counterexamples with definite gaps.

Artifacts: `results/data/orbit_check.json`, `results/data/gauge_min.jsonl`,
`results/data/cascade_exact.json`. Nothing here needed a new computation; it is
a reading of results produced for `docs/gauge_minimization.md` and
`docs/cascade_census.md`.

## The statement

For 48 vertices of the fermionic occupation polytopes of ranks 8 to 10 we have,
side by side:

- a certified extremal state psi with 1-RDM exactly diag(lambda), of Slater
  support s, where s is EXACTLY the minimum over psi's entire natural-orbital
  gauge family U(m_1) x ... x U(m_k) (an attained gauge-invariant rank bound,
  integer arithmetic, no search involved); and
- an exactly certified state phi with spec(rho_phi) = lambda and Slater support
  s' < s, whose 1-RDM is not diagonal, with rational squared weights and an
  exact characteristic-polynomial identity.

Gaps: s - s' = 1 for 38 of them and 2 for 10 of them. By system: (3,8) 4,
(3,9) 5, (3,10) 12, (4,9) 7, (4,10) 7, (5,10) 13.

Numerically phi is psi in rotated orbitals: their 2-RDM spectra, which are
invariant under every one-body unitary, agree to ~1e-16 in all 48 cases, and
phi is reached from psi by a continuous path. Under that identification the
counterexample is the strong one: the SAME state has a strictly sparser
determinant expansion in a non-natural orbital basis than in any of its natural
orbital bases. Without it, the pair still shows that the vertex is attainable
with fewer determinants than any natural-orbital representative of psi achieves.

## Why the "any of its natural orbital bases" is the load-bearing part

A natural-orbital basis is not unique when occupations are degenerate, and every
one of these spectra is degenerate. The comparison would be worthless if s were
merely the support of some arbitrarily chosen natural-orbital representative,
since a different choice inside the degenerate blocks might do better. It does
not, and that is certified rather than searched: the compound of a
block-diagonal unitary preserves each determinant's block profile, hence the
flattening rank of each profile class, so

    support of any gauge representative >= sum of the class flattening ranks,

and for these 48 the library support attains it. So the gap is not an artifact
of a bad gauge choice; it is a gap against the whole gauge family.

## Precision about what is exact and what is not

- EXACT: the value s and its minimality over the gauge family (integer ranks of
  exactly built matrices); the state phi, its support s', its rational squared
  weights, and spec(rho_phi) = lambda (rational and radical arithmetic).
- NUMERICAL: the identification of phi with psi up to a one-body rotation
  (2-RDM spectra to ~1e-16; at the denominator-34 vertex the connecting unitary
  is solved for explicitly, residual 7.5e-14).
- NOT CLAIMED: that s is the minimum over ALL natural-orbital states attaining
  lambda (that is s_Q^NO, a vertex-level quantity bounded below only by s_I),
  nor that s' is the free-basis minimum. Both are upper-bound statements about
  the states in hand.

## The one case that is not in the 48

The denominator-34 vertex of (3,10), which has the largest gap (support 14 in
the natural-orbital basis against 10 in a rotated one), is deliberately excluded
from the count. Its natural-orbital support of 14 is NOT certified minimal over
its gauge family: the lower bound there is 9 and the search closed no part of
the interval. Its gap of 4 is real but weaker in kind, so it is quoted
separately and never inside the 48.

## The orbitals ship with the counts

A referee in this field will ask for the orbitals, not just the numbers. Every
one of the 63 same-state endpoints, the 48 counterexamples among them, carries
its explicit one-body rotation in `results/data/orbit_check.json`: the matrix U
with Gamma(U) psi_library = psi_sparse, as real and imaginary parts, with its
residual and unitarity error. All 63 solves verified, worst residual 1.7e-15,
worst unitarity error ~4e-15. So each counterexample is reproducible end to end:
take the library state, apply the shipped U, read off the shorter expansion.

## What a write-up would still need

The identification of the two states remains numerical (residuals ~1e-15,
above). Making it exact would mean recognizing U itself in closed form, which
has not been attempted; the natural route is the same integer-relation
machinery already used on the weights. Nothing else is missing.

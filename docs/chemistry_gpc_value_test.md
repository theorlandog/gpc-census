# Is GPC geometry useful for quantum chemistry? A falsification test

Artifact: `results/data/chemistry_gpc_value_test.json`.
Generator: `scripts/chemistry_gpc_value_test.py`.
Guard test: `tests/test_chemistry_gpc_value_test.py`.
Plots: `plots/`.

This document answers one question and does not hedge it. Two products were on
trial:

1. **Certified active-space sufficiency.** Does quasipinning bound the omitted
   wavefunction norm, or the energy error, in a way a chemist could rely on?
2. **GPC-guided CI/CSF reduction.** Do generalized Pauli constraints remove
   enough determinants to matter, after particle number, Sz, S^2, point group
   and occupation-based freezing have already been applied?

Both were tested against baselines that were chosen to be hard to beat, and
both fail. The verdict is at the end, with the reasons.

## What was built, and how it was checked

`scripts/chemistry_gpc_value_test.py` is self-contained: a spin-orbital full-CI
engine, natural-orbital diagonalization, exterior-power transport of the CI
vector into the natural-orbital basis, evaluation of every published GPC row,
construction of the pinned selection-rule subspace, and the baselines.

Checks that had to pass before any number below was believed:

- the CI energy reproduces PySCF's FCI to machine precision (1e-15 Eh on H4
  in STO-3G);
- the rotated CI vector has an exactly diagonal 1-RDM whose diagonal is the
  sorted occupation vector (`diag_residual`, `lambda_residual` in every
  record, all below 1e-9);
- the rotated integrals reproduce the same energy for the rotated state, so
  state and Hamiltonian are never in different bases;
- the analytic gradient of the facet slack agrees with finite differences to
  7e-10.

Active spaces were chosen so that `(N, d)` has a trusted system in this
repository, with `d = 2 * n_active_orbitals`. Odd `d` is unreachable with
collinear spin, so it is covered by spinless model Hamiltonians, which are
labelled `model` and never quoted as chemistry evidence.

## Experiment 0: the selection rule needs a hypothesis chemistry violates

Everything downstream rests on one claim: if the occupation vector saturates a
valid row `coeffs . lambda <= rhs`, then in the natural-orbital basis every
determinant carrying amplitude satisfies
`D(I) = rhs - sum_{i in I} coeffs_i = 0`.

The claim is **true when the natural spectrum is nondegenerate**, and the proof
is the variational argument: a pinned state minimizes `D(lambda)` over the unit
sphere, first-order perturbation theory gives `d lambda_i = d gamma_ii` only if
`lambda_i` is simple, and stationarity then forces `D(n) psi = mu psi` with
`mu = D(lambda) = 0`. Measured: pushing a random state onto a facet by
minimizing the slack lands at slack below 1e-15 with omitted weight below 1e-14,
at **17 of 17** facets tested across `(3,6)`, `(3,7)`, `(3,8)` and `(4,8)`, with
minimum spectrum gaps from 8e-3 to 8e-2.

The hypothesis is not decorative. Replayed against the census of certified
extremal states, the rule holds at only **49 of 643** attainers that pass the
diagonality gate, and still fails at **421 of 643** after searching every
relabeling inside degenerate blocks. Every one of the 643 has a degenerate
spectrum, so the census cannot supply a single nondegenerate test case.

This matters for chemistry for an exact reason: **a closed-shell state has an
exactly twofold degenerate spin-orbital spectrum**, `lambda = (n_1/2, n_1/2,
n_2/2, n_2/2, ...)`. The regime where the selection rule is a theorem is
precisely the regime a singlet is never in. Throughout this study the omitted
weight is therefore reported as the **minimum over all relabelings inside
degenerate blocks**, which is generous to the method: a user without the exact
state cannot know which tie-break to take.

## Experiment 4 (stated early, because it is the sharpest result): no rigorous bound exists

The identity that any certification would have to use is exact:

    D(lambda) = sum_I |c_I|^2 D(I),

because `lambda_i = sum_{I contains i} |c_I|^2` in the natural-orbital basis. It
bounds the weight outside the pinned subspace **only if every `D(I)` is
nonnegative**, in which case `W_out <= D / min_{I not in S} D(I)` is a proof.

Checked exhaustively over the twelve systems the repository knows, `(3,6)`
through `(7,10)`, **686 rows in total: not one has all `D(I)` nonnegative.**
Every facet admits determinants with `D(I) <= -1`, so weight can sit outside the
pinned subspace while pushing the slack *down*. There is no elementary bound to
be had at any facet of any known system. (`selection_subspace_scaling` in the
artifact; recomputed from scratch by the guard test.)

That the loophole is real, and not merely unclosed, is shown in Experiment 3.

## Experiment 5, part one: the flagship (3,6) pinning is a theorem about Sz

The Borland-Dennis system is where quasipinning is usually demonstrated. In this
study every three-electron doublet ground state in a three-orbital active space
came out **exactly** pinned, at slack 1e-16, at every geometry of H3 and for 25
of 25 random spin-free Hamiltonians. That is not a GPC effect. It is three lines
of linear algebra:

With `n_alpha = 2` and `n_beta = 1` in three spatial orbitals, the state is a
pure state of one **alpha hole** and one **beta particle**, that is a `3 x 3`
bipartite state. Schmidt's theorem then forces

    lambda^alpha = 1 - lambda^beta  (elementwise),

verified to 8e-16 on 300 random states of the sector, and the Borland-Dennis
slack collapses to

    D = max(1 - 2 p_1, 0),

`p_1` the largest Schmidt coefficient, verified to 8e-16. Any state whose
leading configuration carries at least half the weight is exactly pinned; 99.3
per cent of random states of the sector already are. The famous
three-determinant carrier is the Schmidt decomposition written in fermionic
notation.

So at `(3,6)` the reduction is real (9 determinants in the Sz sector down to 3)
but it is **entirely attributable to fixed particle number and fixed Sz**. It is
also confined to a three-orbital active space, where a full CI is free.

## Experiment 5, part two: exact pinning that removes nothing

At larger systems the free saturation does not disappear, it becomes empty.
Sampling random states of each `(N, d, Ms)` sector and asking which GPC rows they
saturate exactly:

- `(3,8)` and `(5,8)` at `Ms = 1/2`: a fifth to a third of random states are
  exactly pinned to some facet, and in **100 per cent** of those cases the
  selection rule keeps the entire symmetry sector. The pinning is real and the
  information content is zero. H3 in a CAS(3,4) is exactly this case: slack
  3e-15, and the GPC subspace is all 24 determinants of the sector.
- `(4,8)`, `(4,10)`, `(6,10)` at `Ms = 0`, the closed-shell chemistry cases: no
  row is saturated by any random state, and the median distance to the nearest
  facet is 0.39 to 0.68. Singlets are far from the boundary as a matter of
  kinematics; physical states get closer only by being nearly single-reference.

The mechanism behind the vacuous cases is visible in the rows: once the alpha
natural orbitals land at sorted positions `P`, the row `sum_{i in P} lambda_i <=
n_alpha` is a genuine published facet that says nothing but "the alpha
occupations sum to `n_alpha`". Pinning to it is the Sz constraint wearing a
polytope's clothes.

## Experiment 6: would a spin-adapted (GPC-spin) polytope change the answer?

This was tested rather than argued. For four electrons in four spatial orbitals
with `S = 0`, the exact body of achievable spatial occupation vectors was mapped
by maximizing `c . n` over the 20-dimensional CSF space along sampled
directions, and compared against the body the existing census already implies
through the doubling map `lambda = (n_i/2, n_i/2)` (an exact LP over the
published rows).

The fixed-`S` body **is** strictly smaller, but barely: the support gaps have
median about 1e-8 and a maximum of about 2e-2 across sampled directions, against
physical quasipinning distances of 1e-2 to 9e-1 in the same system. The measured
gaps are upper bounds, since the fixed-`S` support is computed by local
maximization, so the true body is if anything closer to the prediction.

A spin polytope would fix one real defect: spatial occupations are generically
nondegenerate, so the selection rule would regain its hypothesis. It does not
touch either binding failure. It does not create a one-sided bound (Experiment
4 is a statement about determinant values, not about which polytope), and it
cannot move a boundary sitting 1e-2 away by an amount of order 1e-8 to 1e-2.
And the correlated regime it would be sold for is exactly where spatial
occupations cluster near 1 and the degeneracy returns.

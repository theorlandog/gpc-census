# Pre-registration: Theorem D and Theorem C' (Track B1, continued)

Standing rules in `docs/prereg_template.md`. Committed together with
`scripts/shared_loop_rationality.py` BEFORE the script is run (rule R1).

Target: remove the checked hypothesis from Theorem C
(`docs/polygon_diagonals.md`) by PROVING it, and lift its conclusion off class
size 3. Theorem C reduced (R1) at large classes to a rational diagonal and
verified the diagonal exists on 17 of 17 census classes, which was POST-HOC.
This work order turns that observation into a theorem plus a purely
combinatorial residue.

## Scope, re-derived

`results/data/polygon_diagonals.json` at `88018a8`: 17 classes of size >= 3,
16 of size 3 and 1 of size 4; 16 distinct states; 11 transport classes (R2).
Of the 54 within-class pairs, 26 have rational `t_ab` and 28 do not. A triage
pass before this pre-registration found that 19 pairs have one-hop-linked
A-side determinants and that the loop-vector identity below holds on all 19;
that triage is what motivated the theorems and is NOT scored here.

## The two theorems, stated before measuring

**Theorem D (shared loop).** Let `a`, `b` be pairs of the one-hop class
`(A, B)` with `u_b = u_a - C + E`. Since both `v` are the same determinants
with `A` swapped for `B`, also `v_b = v_a - C + E`, so `(u_a, u_b)` and
`(v_a, v_b)` are both pairs of the class `(C, E)`. Their Lemma 1 loop vectors,

    e_{v_a} - e_{u_a} - e_{v_b} + e_{u_b}   and   e_{u_b} - e_{u_a} - e_{v_b} + e_{v_a},

are identical, so the two classes SHARE the loop and its holonomy may be read
in either. If `(C, E)` has size 2, Theorem A applies there and yields

    t_ab = s_a s_b (|rho_CE|^2 - p_{u_a} p_{u_b} - p_{v_a} p_{v_b}) / (2 sigma_1 sigma_2),

which is rational under (H1) and (H2) because the monomial
`M = p_{u_a} p_{u_b} p_{v_a} p_{v_b}` equals `(r_a r_b)^2` and its square root
cancels identically.

**Theorem C' (rank-2 Gram, any size).** The sides `z_a` of a class lie in the
real plane, so the Gram matrix `G = (t_ab)`, `t_aa = r_a^2`, has rank at most
2. A class of size `m` therefore has `m - 2` shape parameters once the closure
condition `sum_a z_a = rho_AB` is imposed, and `m - 2` rational diagonals
determine it. Concretely, if the pairs admit a RATIONAL CHAIN, an ordering
whose partial sums `S_k` all have `|S_k|^2` rational, then each consecutive
triangle `(S_{k-1}, z_k, S_k)` has rational squared sides, so both the cosine
and the sine of its angle lie in `Q(sqrt(N_k), sqrt(N_k (N_k - g_k^2)))` with
`N_k = |S_{k-1}|^2 r_k^2` and `g_k = Re(conj(S_{k-1}) z_k)` rational. Composing
angles keeps cosines and sines inside the compositum, so every `cos(Phi_ab)`
lies in a MULTIQUADRATIC field and its Galois group is elementary abelian of
exponent 2: C4 and D4 are excluded at every class size.

Neither theorem is a prediction. What is pre-registered is whether the census
satisfies their hypotheses and whether the machine check of Theorem D's
explicit formula reproduces the independently measured `t_ab`.

## Predictions, committed before the run

| id | prediction | expectation |
|---|---|---|
| P-D-1 | the shared-loop identity holds on every one-hop-linked pair | PASS |
| P-D-2 | for every pair whose linking class has size 2, Theorem D's formula reproduces the measured `t_ab` exactly | PASS |
| P-D-3 | every pair with a size-2 linking class has rational `t_ab` | PASS |
| P-D-4 | every class of size >= 3 contains at least one pair with a size-2 linking class | **FAIL expected**, one class short |
| P-D-5 | every class of size >= 3 admits a rational chain | PASS |

Reasoning, so a hit is not read as a lucky guess:

- **P-D-1** is combinatorial and was already checked on 19 pairs in triage. It
  is re-registered because the script recomputes it from the determinants
  rather than reusing the triage, and a failure would mean the orientation
  convention is wrong.
- **P-D-2** is the substance of Theorem D. The formula is built from the
  LINKING class alone and compared against a `t_ab` computed in the ORIGINAL
  class, so agreement is a genuine cross-check and not an identity by
  construction. A mismatch would most likely be a fermionic sign error in
  `sigma_1 sigma_2` versus `s_a s_b`, which is exactly the kind of convention
  bug this repository has been bitten by before.
- **P-D-3** follows from P-D-2 and is separated so that a sign error shows up
  as one failure rather than two.
- **P-D-4 is pre-registered as a FAILURE**, which is the honest position:
  triage already showed (3,10) v103 channel (1,6) has no size-2 linking class,
  its only one-hop-linked pair sitting in the size-4 class. Registering the
  expected exception in advance is the difference between a scoped theorem and
  a claim quietly trimmed to fit. Per rule R4 this is ALREADY FALSIFIED by data
  in hand and is recorded as scoped, not discovered.
- **P-D-5** is the hypothesis of Theorem C'. It is weaker than P-D-4 because a
  chain may use diagonals that are rational for other reasons, and (3,10) v103
  is fully rational, so it should have a chain even without a size-2 link.

## What would REFUTE the programme

A class of size >= 3 with NO rational chain. Theorem C' would then say nothing
about it, Theorem C would not apply either, and the multiquadratic conclusion
at that class would have no proof, only the measured degree.

## Nonclaims

- Theorem D proves the diagonal is rational; it does not prove that a size-2
  linking class EXISTS. That is the combinatorial residue and it is open.
- Theorem C' concludes multiquadratic, which is what 2-elementarity needs. It
  does not by itself give the literal (R1) bound `deg_Q(cos^2) <= 2` at large
  class size, since a long chain can enlarge the field; the census degrees are
  1 and 2 throughout, which is measured, not proved.
- All arithmetic is on the SHIPPED representative from `states.jsonl`, never
  the numerical natural-orbital ledger.

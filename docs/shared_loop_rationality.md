# Theorem D and Theorem C': the induction on class size

Status: **PROVED**, two theorems, with the census consequence upgraded from
checked to proved on **16 of 17 large classes**. Theorem C
(`docs/polygon_diagonals.md`) needed a rational diagonal and verified one
existed; Theorem D proves the diagonal is rational, and Theorem C' lifts the
conclusion off class size 3. The remaining class is named, and the naive
combinatorial conjecture that would have closed it is **REFUTED** by that same
class.

Pre-registration: `docs/prereg_shared_loop_rationality.md`, committed with the
script at `7e112cc` before the run. Certifying artifacts:
`scripts/shared_loop_rationality.py`,
`results/data/shared_loop_rationality.json`,
`tests/test_shared_loop_rationality.py`.

## Theorem D (shared loop), PROVED

> **Theorem D.** Assume (H1) and (H2). Let `a`, `b` be pairs of the one-hop
> class `(A, B)` whose A-side determinants are one hop apart,
> `u_b = u_a - C + E` with `C, E` outside `{A, B}`. Then:
>
> 1. `v_b = v_a - C + E`, so `(u_a, u_b)` and `(v_a, v_b)` are both pairs of
>    the one-hop class `(C, E)`;
> 2. the Lemma 1 loop of `{a, b}` in `(A, B)` and the Lemma 1 loop of that pair
>    in `(C, E)` are the SAME integer vector;
> 3. if `(C, E)` has size 2 then
>
>        t_ab = s_a s_b (|rho_CE|^2 - p_{u_a} p_{u_b} - p_{v_a} p_{v_b}) / (2 sigma_1 sigma_2)
>
>    which is RATIONAL.

*Proof.* (1) Both `v` are their `u` with `A` replaced by `B`, and that
substitution commutes with replacing `C` by `E`, so `v_b = v_a - C + E`. Since
`C in u_a, v_a` and `E in u_b, v_b`, both determinant pairs belong to `(C, E)`.

(2) Lemma 1 gives the `(A, B)` loop as `e_{v_a} - e_{u_a} - e_{v_b} + e_{u_b}`
and the `(C, E)` loop as `e_{u_b} - e_{u_a} - e_{v_b} + e_{v_a}`. These are the
same formal sum, so the two classes SHARE the loop and its holonomy may be read
in either.

(3) At size 2, Theorem A applied to `(C, E)` gives

    cos(Phi) = (|rho_CE|^2 - p_{u_a} p_{u_b} - p_{v_a} p_{v_b}) / (2 sigma_1 sigma_2 sqrt(M)),
    M = p_{u_a} p_{u_b} p_{v_a} p_{v_b}.

By (2) this same `Phi` is the `(A, B)` holonomy `Phi_ab`. And
`r_a r_b = sqrt(p_{u_a} p_{v_a} p_{u_b} p_{v_b}) = sqrt(M)`, the identical
monomial. So in `t_ab = s_a s_b r_a r_b cos(Phi_ab)` the two square roots
cancel identically, leaving a quotient of rationals. ∎

**Why this matters.** Theorem C's hypothesis stops being an observation. The
rational diagonal is Theorem A, read through a shared loop from a smaller
class, and the cancellation is exact rather than numerical.

## Theorem C' (rank-2 Gram at any size), PROVED

> **Theorem C'.** The sides `z_a` of a one-hop class lie in the real plane, so
> the Gram matrix `G = (t_ab)` with `t_aa = r_a^2` has RANK AT MOST 2. If the
> pairs admit a RATIONAL CHAIN, an ordering whose partial sums `S_k` all have
> `|S_k|^2` rational, then every `cos(Phi_ab)` lies in a multiquadratic field
> and its Galois group is elementary abelian of exponent 2. C4 and D4 are
> excluded at EVERY class size.

*Proof.* Along a chain, `|S_k|^2 = |S_{k-1}|^2 + r_k^2 + 2 g_k` with
`g_k = Re(conj(S_{k-1}) z_k)`, so each `g_k` is rational. The angle between
`S_{k-1}` and `z_k` therefore has

    cos = g_k / sqrt(N_k),   sin = sqrt(N_k - g_k^2) / sqrt(N_k),   N_k = |S_{k-1}|^2 r_k^2,

with `N_k` rational, so both lie in `Q(sqrt(N_k), sqrt(N_k(N_k - g_k^2)))`,
whose radicands are rational. Every `z_a` differs from `z_1` by a sum of such
angles, and cosines and sines of sums stay in the compositum of these fields.
Finally `t_ab = r_a r_b cos(Phi_ab)` with `(r_a r_b)^2` rational adjoins at
most one more square root of a rational. Every radicand is rational, so the
compositum is multiquadratic. ∎

**What replaces the quadrilateral picture.** Theorem C was proved from one
vanishing 3x3 Gram determinant. The general statement is that the whole Gram
matrix has rank at most 2 for the same reason, and `m - 2` is simply the shape
count after the closure condition. Proposition B tried to reach large classes
by iterating a reduction and only ever produced a 2-tower, which may close
dihedrally; C' produces a multiquadratic field directly because every radicand
it adjoins is rational.

## The induction, and how far it reaches

    Theorem A  (size <= 2)                     base case
        |
        | shared loop, Theorem D
        v
    rational diagonals at a larger class
        |
        | rank-2 Gram, Theorem C'
        v
    multiquadratic holonomy field, any size

Theorem C' needs `m - 2` rational diagonals. Theorem D supplies one per pair
carrying a size-2 link. Where D supplies enough, the class is proved end to end
from Theorem A with **no checked arithmetic hypothesis at all**.

| | classes | proved end to end | exactly at the minimum `m-2` |
|---|---|---|---|
| size 3 | 16 | 15 | 15 |
| size 4 | 1 | 1 | 1 |
| **total** | **17** | **16** | **16** |

Every one of the 16 sits EXACTLY at `m - 2`: never a spare diagonal. The size-4
class needs two and has exactly two.

## Scorecard

Scope 17 classes, 16 distinct states, **11 transport classes** (rule R2).

| id | prediction | expectation | verdict |
|---|---|---|---|
| P-D-1 | the shared-loop identity holds on every one-hop-linked pair | PASS | **PASS**, 19 of 19 |
| P-D-2 | Theorem D's formula reproduces the measured `t_ab` exactly | PASS | **PASS**, 17 of 17 |
| P-D-3 | every size-2-linked pair has rational `t_ab` | PASS | **PASS**, 17 of 17 |
| P-D-4 | every class of size >= 3 has a size-2 link | **FAIL expected** | **FAIL**, 16 of 17 |
| P-D-5 | every class of size >= 3 admits a rational chain | PASS | **PASS**, 17 of 17 |

P-D-2 is the one that carries the theorem, and it is a genuine cross-check: the
predicted value is built from the LINKING class `(C, E)` alone, while `t_ab` is
computed inside `(A, B)`. A fermionic sign error in `sigma_1 sigma_2` against
`s_a s_b` would show up here, and does not.

**A measurement bug caught before scoring.** The first run reported (3,10) v98
as having no rational chain. That was wrong: the chain search skipped orderings
with `order[0] > order[-1]` to "dedupe reverses", but reversing an ordering
changes which partial sums appear, so the reverse of a chain is generally not a
chain. With every ordering searched, P-D-5 is 17 of 17. The bug is recorded
because it would have manufactured a refutation of the programme's own
hypothesis.

## The refuted conjecture, which is the useful negative

P-D-4 was pre-registered as an expected failure and failed. The class

    (3,10) v103 channel (1,6)

has NO pair with a size-2 link: its only one-hop-linked pair sits in the size-4
class (3,9) of the same state. So the clean combinatorial conjecture that would
have closed the induction,

> every one-hop class of size >= 3 contains a pair linked by a class of size
> at most 2,

is **FALSE**, and the census already contains its counterexample. Recording
this matters because that conjecture was the obvious next thing to try to
prove, and it is dead.

What saves the class is that all three of its diagonals are rational anyway, so
Theorem C applies with a checked hypothesis. And the reason is visible: its
linking class (3,9) has every pairwise `t` rational, which suggests the honest
generalisation of Theorem D replaces "size 2" with "all pairwise `t` rational",
turning the induction into a mutual recursion between classes whose
well-foundedness is the new open question.

## What this does and does not license

**Does.** Multiquadratic holonomy fields, hence 2-elementary Galois groups with
C4 and D4 excluded, at 16 of the 17 large census classes, PROVED from Theorem A
with no checked arithmetic hypothesis. Theorem C' holds at arbitrary class
size, so nothing here is special to 3 or 4.

**Does not.** It does not close the last class, whose Theorem C hypothesis is
still checked. It does not prove that a rational chain always exists; that is
verified 17 of 17 and is a hypothesis, not a theorem. Theorem C' concludes
multiquadratic, which is what 2-elementarity requires; it does not by itself
give the literal (R1) bound `deg_Q(cos^2) <= 2` at large size, since a long
chain can enlarge the field. The measured degrees are 1 and 2 throughout, which
is measurement.

**The residue, restated.** Prove that the class-linking recursion is
well founded: that every class of size >= 3 either carries a size-2 link or is
linked to a class whose pairwise `t` are already known rational, with no
cycles. The census satisfies this, and (3,10) v103 shows the one-step version
is not enough.

## Reproducing

```sh
python scripts/shared_loop_rationality.py     # about 8 minutes, exact throughout
uv run pytest tests/test_shared_loop_rationality.py
```

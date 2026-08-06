# The polygon reformulation and Theorem C (Track B1)

Status: **REDUCED**, with a **PROVED** theorem inside it. Condition (R1) of
`docs/holonomy_two_elementary.md` at one-hop class sizes 3 and 4 is reduced to
a single exactly checkable per-class condition, the RATIONAL DIAGONAL, and
Theorem C below proves that the condition suffices. The condition holds on
17 of 17 census classes, which is a POST-HOC observation and is labelled so.

Pre-registration: `docs/prereg_polygon_diagonals.md`, committed with the
measuring script at `221fab8` before the run. Certifying artifacts:
`scripts/polygon_diagonals.py`, `results/data/polygon_diagonals.json`,
`tests/test_polygon_diagonals.py`.

## The reformulation

Lemma 2 of `docs/holonomy_rationality.md` expands `|rho_AB|^2`. Left unexpanded
the channel condition is

    rho_AB = sum_a z_a,     z_a = s_a conj(c_{u_a}) c_{v_a},   |z_a| = r_a,

so a one-hop class of size `m` is a CLOSED (m+1)-GON with prescribed side
lengths `r_1, ..., r_m, |rho_AB|`, and the within-class holonomies are its
angles. Theorem A is the law of cosines for a triangle, which is why its
radicand is a weight monomial. A class of size 3 is a quadrilateral, and a
quadrilateral splits into two triangles along a diagonal.

## Theorem C (triangulation), PROVED

> **Theorem C.** Assume (H1) rational weights and (H2) rational channel
> targets. Let `(A, B)` be a one-hop class of size 3 with pairs indexed
> `a, b, c`, and suppose `t_ab = Re(z_a conj(z_b))` is RATIONAL, equivalently
> the diagonal `|z_a + z_b|^2` is rational. Then `t_ac` and `t_bc` are roots of
>
>     D^2 t^2 - 2 p q t - (r_a^2 r_c^2 D^2 - r_a^2 p^2 - r_c^2 q^2) = 0,
>
> with `D^2 = r_a^2 + r_b^2 + 2 t_ab`, `q = r_a^2 + t_ab`, and
> `p = (|rho_AB|^2 - sum_a r_a^2)/2 - t_ab`, all rational. Hence
> `t_ac` lies in `Q(sqrt(Delta))` with
> `Delta = (r_a^2 D^2 - q^2)(r_c^2 D^2 - p^2)` rational, and
>
>     cos(Phi_ac) = t_ac / (s_a s_c r_a r_c)  lies in  Q(sqrt(Delta), sqrt(r_a^2 r_c^2)),
>
> a MULTIQUADRATIC field. In particular `deg_Q(cos^2 Phi_ac) <= 2`, so (R1)
> holds, the Galois group is contained in V4, and C4 and D4 are excluded.

*Proof.* Put `D = z_a + z_b`. The three vectors `z_a`, `z_c`, `D` lie in the
real plane `C = R^2`, so they are linearly dependent and their Gram matrix

    [ r_a^2   t_ac    q   ]
    [ t_ac    r_c^2   p   ]
    [ q       p       D^2 ]

is singular. Expanding `det = 0` gives the displayed quadratic. The entries are
rational: `D^2` and `q = Re(conj(z_a) D) = r_a^2 + t_ab` by hypothesis and
(H1); `p = Re(conj(D) z_c) = t_ac + t_bc` because Lemma 2 reads
`|rho_AB|^2 = sum_a r_a^2 + 2(t_ab + t_ac + t_bc)` and `|rho_AB|^2` is rational
by (H2) while `t_ab` is rational by hypothesis. The discriminant of the
quadratic is `4 (r_a^2 D^2 - q^2)(r_c^2 D^2 - p^2)`, so `t_ac` is at most
quadratic over Q with a rational radicand. Dividing by `s_a s_c r_a r_c`, whose
square is rational, adjoins at most one further square root of a rational. Both
radicands are rational, so the compositum is multiquadratic and its Galois
group is elementary abelian of exponent 2. ∎

**Why this is the mechanism and not more bookkeeping.** Proposition B also
lands in a quadratic extension, but its radicand `alpha^2 + beta^2 - gamma^2`
involves `sin psi` from an earlier loop, so the tower it builds is a 2-tower
and a 2-tower may close dihedrally. Theorem C's radicand is rational outright,
because the Gram relation is between the polygon's own side lengths and never
leaves the class. That is the difference between a 2-tower and a multiquadratic
field, which is exactly the open gap.

**Machine-checked.** The certificate in `polygon_diagonals.json` builds the
quadratic from rational data only and then substitutes the independently
measured `t_ac`. The Gram residual is exactly 0 on all 17 classes, every
coefficient is rational on all 17, and every discriminant is rational on all
17. Fourteen classes take the `gram_triangulation` route. The other three,
(3,10) v75 and both classes of (3,10) v103, have every pairwise `t` rational
already, so each `cos(Phi_ab)` sits in a single quadratic extension by a weight
monomial and there is nothing to triangulate; they are certified by the
`fully_rational_polygon` route.

## Scorecard

Scope 17 classes, 16 distinct states, **11 transport classes** (rule R2: the
independent content is 11, not 17). 58 proper subsets, 0 with an unresolved
degree.

| id | prediction | expectation | verdict |
|---|---|---|---|
| P-B1-1 | some proper partial sum is irrational | IRRATIONAL FOUND | **PASS**, 28 of 58 |
| P-B1-2 | every mode-coherent subset is rational | PASS | **FAIL**, 19 of 25; 6 irrational |
| P-B1-3 | mode-coherent partial sums equal the 2-RDM element | PASS | **PASS**, 25 of 25 |
| P-B1-4 | no mode-incoherent subset is rational | PASS | **FAIL**, 11 of 33 are rational |
| P-B1-5 | every irrational value has degree exactly 2 | PASS | **PASS**, 28 of 28 |
| P-B1-6 | `t_ab` rational iff `cos^2` rational | PASS | **PASS**, 0 disagreements in 54 |

Four PASS, two FAIL. The prereg named P-B1-2 as the one that mattered, and it
is one of the two that failed.

## The failure, which is the informative part

P-B1-2 predicted that a mode-coherent partial sum is rational because it is a
2-RDM element and (H2) makes 1-RDM channel targets rational. P-B1-3 confirms
the identification exactly, 25 of 25: those partial sums really are
`rho^(2)_{AC,BC}`. And six of them are IRRATIONAL, of degree 2, at (5,10) v256,
v261 and v268.

So the spectator-mode decomposition is real and arithmetically useless, which
is precisely the refutation condition the pre-registration named in advance.
The standalone finding is worth recording on its own terms:

> **(H2) does not extend from the 1-RDM to the 2-RDM.** Channel targets
> `|rho_AB|^2` are rational across the census, but the spectator-mode 2-RDM
> entries `rho^(2)_{AC,BC}` reached from them are not: 6 of 25 measured are
> irrational of degree 2.

P-B1-4 failed in the other direction: 11 of 33 subsets with no spectator mode
are rational anyway, so mode-coherence is neither necessary nor sufficient for
rationality. It is the wrong invariant, and the polygon's own geometry is the
right one.

## The post-hoc observation, labelled as such

**Every one of the 17 classes carries at least one rational diagonal.** This
was NOT pre-registered. It is what makes Theorem C bite, and it is exactly the
kind of pattern the mortality table in `docs/closing_ledger.md` says dies: ten
mined laws scored, ten failed. It is recorded as a post-hoc pattern needing an
out-of-sample test that only Stage 1 can supply.

Two things keep it from being a coincidence of counting. The condition is not
generic: 28 of the 54 within-class pairs are irrational, and 14 of the 17
classes carry EXACTLY ONE rational diagonal, so almost every class sits at the
minimum that makes the theorem apply and no more. And the condition is not
vacuous: it is exactly what a quadrilateral needs to be split.

| classes | with a rational diagonal | with exactly one |
|---|---|---|
| 17 | 17 | 14 |

## What this does and does not license

**Does.** (R1) at class size 3 follows from a rational diagonal, proved, with
C4 and D4 excluded rather than merely absent. Measured directly, `cos^2` has
degree 1 on 26 and degree 2 on 28 of the 54 within-class pairs, so (R1) holds
on every within-class holonomy in the census.

**Does not.** It does not prove (R1) census-wide at sizes 3 and 4, because the
rational-diagonal hypothesis is checked and not proved. It does not touch the
loop-basis holonomies directly: Theorem C is about within-class angles, and
carrying it to a basis of `L` needs the sines as well as the cosines, which
lie in the same multiquadratic field but have not been certified here. It does
not extend Theorem A by bookkeeping, and nothing here reads a holonomy on a
non-integer loop basis.

**The residue, stated as the new open problem.** Prove that every one-hop class
of size 3 carries a rational diagonal, or exhibit a class that does not. That
is a sharper statement than "(R1) at class sizes 3 and 4" and it is decidable
per class in exact arithmetic.

## Reproducing

```sh
python scripts/polygon_diagonals.py        # about 6 minutes, exact throughout
uv run pytest tests/test_polygon_diagonals.py
```

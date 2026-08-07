# Pre-registration: Theorem E and closing the last class

Standing rules in `docs/prereg_template.md`. Committed with
`scripts/class_closure.py` BEFORE the run (rule R1).

Target: the one large class Theorem D could not reach, (3,10) v103 channel
(1,6), and the question of whether 16 of 17 was hiding a generalisation.

## The suspicion, stated before measuring

16 of 17 is a suspicious score. Theorem D proved a rational diagonal by reading
Theorem A through a loop shared with a size-2 class, and the census then
supplied exactly `m - 2` such diagonals per class and never one more. A
mechanism that lands exactly on the minimum at every class, and misses by
exactly one class, is more likely to be a special case of something than to be
the real statement.

## Theorem E, derived before measuring

Not a prediction. Within a one-hop class, the sides lie in the real plane, so
the Gram matrix has rank at most 2 and its 3x3 minor on any triple `{a, b, c}`
vanishes:

    r_b^2 t_ac^2 - 2 t_ab t_bc t_ac + (r_a^2 t_bc^2 + r_c^2 t_ab^2 - r_a^2 r_b^2 r_c^2) = 0.

With `t_ab` and `t_bc` rational this is a rational quadratic whose discriminant
factors as

    Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2)
              = r_a^2 r_b^4 r_c^2 sin^2(Phi_ab) sin^2(Phi_bc).

So `t_ac` lies in `Q(sqrt(Delta))` with `Delta` RATIONAL, hence multiquadratic,
always. It is rational exactly when `Delta` is a rational square, in particular
whenever either linked holonomy has `cos = +/- 1`.

This supersedes the size-2 hypothesis of Theorem D: sharing a loop identifies
Gram entries regardless of the linking class's size, and size 2 was only ever
the base case that made one entry rational to begin with.

## Predictions, committed before the run

| id | prediction | expectation |
|---|---|---|
| P-E-1 | every Theorem E quadratic has rational coefficients and rational discriminant | PASS |
| P-E-2 | the factorisation `Delta/4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2)` holds identically | PASS |
| P-E-3 | the measured `t_ac` solves the quadratic built from `t_ab`, `t_bc` alone | PASS |
| P-E-4 | at (3,10) v103 channel (3,9) the chain from the two Theorem-D diagonals has `Delta = 0`, forcing `t_03` rational | PASS |
| P-E-5 | with that, all 17 large classes close | PASS |
| P-E-6 | degenerate triples (`Delta = 0`) occur ONLY at v103 | **FAIL expected** |

Reasoning:

- **P-E-1 to P-E-3** are the machine check of Theorem E. P-E-3 is the
  non-circular one: the coefficients are built from `t_ab` and `t_bc` and the
  measured `t_ac` is substituted afterwards.
- **P-E-4** is the closure. The prediction is not that `Delta` is some square
  but that it is exactly ZERO, because v103 is real up to gauge with holonomy
  `(0, pi, pi, 0, 0)`, certified in `results/data/v103_fiber_count.json`. Both
  feeding holonomies should have `cos = -1`, both sines vanish, and the
  quadratic degenerates to a perfect square. If `Delta` were a nonzero square
  the class would still close but the stated MECHANISM would be wrong, so the
  distinction is registered.
- **P-E-6 is registered as an expected FAILURE.** Other census states have
  holonomies with `cos = +/- 1` (the degree-1 population is 37 of 82 loops), so
  degenerate triples should appear elsewhere too. Predicting exclusivity and
  finding it false is more informative than not asking.

## What would REFUTE this

`Delta` nonzero and a non-square at v103 channel (3,9). Theorem E would then
give only multiquadratic there, `t_03` would be irrational, Theorem D' could
not transport rationality to channel (1,6), and that class would stay open with
its Theorem C hypothesis merely checked.

## Nonclaims

- Theorem E gives `t_ac` in a quadratic extension by a rational radicand. That
  is multiquadratic and suffices for 2-elementarity; it is NOT a claim that
  `t_ac` is rational in general, and the census contains triples where it is
  not.
- Closing all 17 census classes does not prove Holonomy Rationality at class
  sizes 3 and 4 in general. The chain still starts at Theorem A and needs a
  linking structure to exist, which is verified per class, not proved.
- All arithmetic reads `states.jsonl` closed forms, never the numerical
  natural-orbital ledger.

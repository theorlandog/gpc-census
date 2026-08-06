# Theorem E and the closure of every large class

Status: **PROVED, 17 of 17.** Every large one-hop class in the census now has
multiquadratic within-class holonomy cosines, grounded in Theorem A with no
checked arithmetic hypothesis anywhere in the chain. C4 and D4 are excluded at
class sizes 3 and 4 census-wide.

Pre-registration: `docs/prereg_class_closure.md`, committed with the script at
`7ecec03` before the run. Certifying artifacts: `scripts/class_closure.py`,
`results/data/class_closure.json`, `tests/test_class_closure.py`.

## 16 of 17 was hiding a generalisation, and it was

Theorem D proved a rational diagonal by reading Theorem A through a loop shared
with a size-2 class. The size-2 condition was never needed for the SHARING. It
was needed only to have something rational to start from. What actually
propagates rationality is the rank-2 Gram matrix, and it does not care where
the rational entries came from.

> **Theorem E (chained holonomy).** Let `a`, `b`, `c` be three pairs of a
> one-hop class. The sides lie in the real plane, so the Gram matrix has rank at
> most 2 and its 3x3 minor on `{a, b, c}` vanishes:
>
>     r_b^2 t_ac^2 - 2 t_ab t_bc t_ac + (r_a^2 t_bc^2 + r_c^2 t_ab^2 - r_a^2 r_b^2 r_c^2) = 0.
>
> If `t_ab` and `t_bc` are rational, this is a RATIONAL quadratic in `t_ac`
> whose discriminant factors:
>
>     Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2)
>               = r_a^2 r_b^4 r_c^2 sin^2(Phi_ab) sin^2(Phi_bc).
>
> So `t_ac` always lies in `Q(sqrt(Delta))`, a quadratic extension by a
> RATIONAL radicand. It is RATIONAL exactly when `Delta` is a rational square,
> in particular whenever either linked holonomy has `cos = +/- 1`.

*Proof.* Rank at most 2 kills every 3x3 minor; expanding the one on `{a,b,c}`
gives the displayed quadratic, whose coefficients are rational under (H1) once
`t_ab` and `t_bc` are. Its discriminant is `4(t_ab t_bc)^2` less
`4 r_b^2 (r_a^2 t_bc^2 + r_c^2 t_ab^2 - r_a^2 r_b^2 r_c^2)`, which collects to
the stated product. Substituting `t_ab = r_a r_b cos(Phi_ab)` turns each factor
into `-r^2 r^2 sin^2`, giving the geometric form. ∎

Theorem E is the statement Theorem D was a special case of, and it is what the
missing class needed.

## The last class, closed

**(3,10) v103 channel (1,6)** was the single class Theorem D could not reach:
its only one-hop link runs to the size-4 class (3,9), not to a size-2 class.
The chain that closes it is:

    Theorem A at (2,6), size 2   -> Theorem D -> t_02 rational   |
    Theorem A at (1,2), size 2   -> Theorem D -> t_23 rational   |  in ch(3,9)
    Theorem E on the triple (0,2,3), Delta = 0 -> t_03 RATIONAL  |
    Theorem D' shared loop ch(3,9) <-> ch(1,6) -> its diagonal rational
    Theorem C -> ch(1,6) multiquadratic

Every input is Theorem A. The artifact records
`uses_only_theorem_D_inputs: true` on that triple, so the step does not quietly
consume a measured rationality.

**Why `Delta` is zero, and it is not luck.** Both feeding holonomies have
`cos(Phi) = -1` exactly, so both sines vanish and the quadratic degenerates to
a perfect square with a single rational root. That is not a coincidence of this
computation: `results/data/v103_fiber_count.json` independently certifies that
v103's fiber is a single point with holonomy `(0, pi, pi, 0, 0)`, so the state
is REAL UP TO GAUGE. A real state has every holonomy at `0` or `pi`, hence
every sine zero, hence every chained discriminant zero. The one class the
size-2 mechanism missed is the one where the geometry is degenerate, and the
degeneracy is the property the repository already recorded about that state.

## Scorecard

Scope 17 classes, 16 distinct states, 11 transport classes (rule R2).

| id | prediction | expectation | verdict |
|---|---|---|---|
| P-E-1 | rational coefficients and rational discriminant | PASS | **PASS**, 18 of 18 |
| P-E-2 | the discriminant factorisation holds identically | PASS | **PASS**, 18 of 18 |
| P-E-3 | the measured `t_ac` solves the quadratic built without it | PASS | **PASS**, 18 of 18 |
| P-E-4 | `Delta = 0` at v103 ch(3,9), forcing `t_03` rational | PASS | **PASS**, both sines vanish |
| P-E-5 | all 17 large classes close | PASS | **PASS**, 17 of 17 |
| P-E-6 | degenerate triples occur ONLY at v103 | **FAIL expected** | **FAIL**, also at (3,10) v75 |

Five PASS, one pre-registered FAIL.

P-E-3 is the non-circular check: the quadratic's coefficients are built from
`t_ab` and `t_bc` only, and the measured `t_ac` is substituted afterwards.

## What the failure taught, which is more than the successes

P-E-6 predicted degeneracy was a v103 peculiarity. It is not. **17 of the 18
chained triples have `Delta = 0`**, spread across (3,10) v75 and both classes
of v103. Degeneracy is the rule among chained triples, not the exception, which
reframes the whole mechanism: rationality propagates mostly because linked
holonomies are pinned at `cos = +/- 1`, not because discriminants happen to be
squares.

The single non-degenerate triple is instructive and worth naming:

| triple | `Delta` | rational square |
|---|---|---|
| (3,10) v75 ch(2,8), `(0,1,2)` | `225/5903156224` | yes, `(15/76832)^2` |

So Theorem E's two modes both occur in the census: 17 triples close by a
vanishing sine, 1 closes by a genuine rational square. All 18 discriminants are
rational squares, which is why every `t_ac` in these classes is rational.

## Closure ledger

| route | classes |
|---|---|
| Theorem D supplies `m - 2` diagonals directly | 16 |
| Theorem E extends rationality along a chain | 1 |
| **total closed** | **17 of 17** |

## What this does and does not license

**Does.** Multiquadratic holonomy fields, hence 2-elementary Galois groups with
C4 and D4 excluded, at every large one-hop class in the census, PROVED from
Theorem A with no checked arithmetic hypothesis. Theorem E holds at any class
size and supersedes the size-2 hypothesis of Theorem D.

**Does not.** It does not prove Holonomy Rationality at class sizes 3 and 4 in
general. The chain still begins at Theorem A and needs a linking structure to
exist; that structure is verified per class, not proved to exist, and (3,10)
v103 already refuted the clean one-step conjecture that every class carries a
size-2 link (`docs/shared_loop_rationality.md`). It does not claim `t_ac` is
rational in general: Theorem E gives multiquadratic always and rationality only
when `Delta` is a square.

**The residue.** Two questions remain, and they are sharper than the one they
replace. First, does every large class admit a chain of shared loops back to a
size-2 class, with no cycles? Second, is the observed dominance of degenerate
triples (17 of 18) a fact about real-up-to-gauge states specifically, or does
it hold wherever chaining applies? The second is a post-hoc pattern from a
sample of 18 in 3 classes and needs an out-of-sample test.

## Reproducing

```sh
python scripts/shared_loop_rationality.py     # produces the input artifact
python scripts/class_closure.py               # seconds, exact throughout
uv run pytest tests/test_class_closure.py
```

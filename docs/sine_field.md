# Theorem F: the sine field

Status: **PROVED** (Theorem F and its propagation identity), with a **census
invariant** measured on all 54 pairs of all 17 large classes. Dissecting the
single non-degenerate triple of `docs/class_closure.md` turned up an object the
arithmetic arc had recorded but never used: the SINE field `Q(sin Phi)`, which
is not the holonomy field `Q(cos Phi)`, and which is what actually governs
Theorem E.

Pre-registration: `docs/prereg_sine_field.md`, committed with the script at
`24f30f8` before the run. Certifying artifacts: `scripts/sine_field.py`,
`results/data/sine_field.json`, `tests/test_sine_field.py`.

## The object

Write `A_ab = Im(z_a conj(z_b))`. Then the Gram entry and the area are the real
and imaginary parts of ONE product:

    z_a conj(z_b) = t_ab + i A_ab,      t_ab^2 + A_ab^2 = r_a^2 r_b^2,

and `A_ab = r_a r_b sin(Phi_ab)` is twice the signed area of the triangle
`(0, z_a, z_b)`.

> **Theorem F.** Theorem E's discriminant is built from areas, not from Gram
> entries:
>
>     Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2) = A_ab^2 A_bc^2,
>
> so `sqrt(Delta)/2 = |A_ab A_bc|`. Hence Theorem E forces `t_ac` RATIONAL if
> and only if `A_ab A_bc` is rational, which happens exactly when
> `sin^2(Phi_ab)` and `sin^2(Phi_bc)` have the same squarefree part.
>
> Moreover the plane sine-addition identity
>
>     r_b^2 A_ac = A_ab t_bc + t_ab A_bc
>
> makes the sine field PROPAGATE: with `t_ab`, `t_bc` rational, `A_ab` and
> `A_bc` in `Q(sqrt d)` put `A_ac` in `Q(sqrt d)` too.

*Proof.* The first identity is `t^2 + A^2 = |z_a conj(z_b)|^2 = r_a^2 r_b^2`
substituted into Theorem E's discriminant. The second is
`sin(x+y) = sin x cos y + cos x sin y` multiplied by `r_a r_b^2 r_c`. ∎

**Corollary.** If the sine field is constant across a class, every Theorem E
discriminant in that class is automatically a rational square, and rationality
propagates through the whole class with no case-by-case squareness check. The
checked condition is replaced by a single invariant `d` per class.

## Why 17 of 18 triples degenerate

A state that is real up to gauge has every holonomy at `0` or `pi`, so every
sine is zero, every area vanishes, and every discriminant is trivially a
square. The 17 degenerate triples are that limit. The 1 non-degenerate triple
is the general shape, which is why it was worth dissecting.

## Scorecard

Scope 17 classes, 16 states, **11 transport classes** (rule R2). The value 7
covers 4 transport classes and 47 covers 3, so the six distinct sine fields are
spread over 11 independent observations, not 17.

| id | prediction | expectation | verdict |
|---|---|---|---|
| P-F-1 | `Delta/4 = A_ab^2 A_bc^2` | PASS | **PASS**, 60 of 60 |
| P-F-2 | the sine-addition identity holds up to orientation sign | PASS | **PASS**, 60 of 60 |
| P-F-3 | the squarefree parts are recorded `sin_radicands` | PASS | **PASS**, 17 of 17 |
| P-F-4 | the sine field is a class invariant, on EVERY pair | PASS | **PASS**, 44 of 44 |
| P-F-5 | the sine field differs from the cosine field somewhere | PASS | **PASS, but weakly**, 1 of 15 |

**P-F-3 is the cross-check that makes this real rather than invented.**
`sin_radicands` is an independently generated field of `holonomy_fields.json`,
produced by a different script for a different purpose. Every squarefree part
computed here is one of those recorded radicands.

**P-F-4 is the substantive result.** It is tested on all 54 pairs, not the 26
with rational `t`. Ten pairs have vanishing area and are vacuous; the other 44
all satisfy `A_ab / sqrt(d) in Q(t_ab)` for the class-constant `d`.

**P-F-5 passed weakly and the write-up says so.** The comparison is against the
UNION of recorded cosine radicands over all loops of the state, which is a
generous set, so containment is easy and the test is not sharp. A per-loop
comparison would be the honest sharpening and has not been run.

## The sine fields, and the one that is invisible

| `d` | classes | transport classes |
|---|---|---|
| 0 (degenerate, real up to gauge) | 2 | 1 |
| 7 | 9 | 4 |
| 15 | 1 | 1 |
| 23 | 1 | 1 |
| 47 | 4 | 3 |
| 71 | 1 | 1 |

**Every nonzero area is irrational: 0 of 44.** So when the sine field is
nondegenerate it is always a genuine quadratic extension, never accidentally
rational.

**The finding worth the dissection.** The single class where the sine field is
NOT contained in the recorded cosine radicands is

    (3,10) v75 channel (2,8),  sine field Q(sqrt 15),  cosine radicands NONE,

and it is the same class that produced the single non-degenerate triple. Both
of v75's loops have degree 1, so `cos(Phi)` is rational and the holonomy census
records the state as arithmetically trivial with no radicands at all. Its sines
nevertheless carry `sqrt(15)`. A programme that reads only `Q(cos Phi)` sees
nothing there; the sine field is where its arithmetic lives.

That is the information the 1 of 18 was carrying. The 17 degenerate triples
have zero sines and tell you nothing about the field; the one non-degenerate
triple sits in the one state whose sine field is invisible to the cosine field.

## A measurement bug, caught before scoring

The first run reported `in_sine_field` holding on **0 of 44** pairs, a perfect
failure. That was `sympy.primitive_element` returning its polynomial in a
symbol of its own choosing while the degree was read against an unrelated
symbol, which returns 0 for every input. A rational-`t` pair must pass this
test trivially, so 0 of 44 was impossible as mathematics and diagnostic as
plumbing. With the degree read off the polynomial's own generator it is 44 of
44. This is the third measurement bug in this line of work caught by asking
whether a result was even possible before believing it.

## What this does and does not license

**Does.** Theorem F is proved and gives Theorem E's discriminant a geometric
meaning: it is a product of triangle areas. The rationality condition is
membership in a common sine field. The propagation identity is proved.

**Does not.** The class invariance of the sine field is MEASURED, 44 of 44 on
11 transport classes, not proved. The propagation identity shows the field is
preserved along chains with rational Gram entries, which is weaker than
invariance across a whole class. P-F-5 is a weak pass and does not establish
that the sine and cosine fields are generically independent.

**The residue.** Prove that the sine field is a class invariant, which by the
corollary would remove the squareness check from Theorem E entirely. And run
the sharp per-loop comparison of sine against cosine radicands, which P-F-5
only approximated.

## Reproducing

```sh
python scripts/sine_field.py         # a few minutes, exact throughout
uv run pytest tests/test_sine_field.py
```

# Pre-registration: the sine field (dissecting the 1 of 18)

Standing rules in `docs/prereg_template.md`. Committed with
`scripts/sine_field.py` BEFORE the run (rule R1).

Target: the single non-degenerate triple of `docs/class_closure.md`. 17 of 18
chained triples closed because a sine vanished; one closed because its
discriminant happened to be a rational square. A 17-to-1 split is a hint that
the 1 carries the general mechanism and the 17 are its degenerate limit.

## What the triage found, and it is not scored here

A triage pass before this pre-registration dissected (3,10) v75 channel (2,8),
the non-degenerate triple. Its two contributing areas are individually
IRRATIONAL, both `sqrt(15)/392`, and their product is rational only because
they carry the SAME radicand. That motivated the theorem below and the
predictions; it is not itself a scored result.

## Theorem F, derived before measuring

Write `A_ab = Im(z_a conj(z_b))`, so that

    z_a conj(z_b) = t_ab + i A_ab,    t_ab^2 + A_ab^2 = r_a^2 r_b^2,

and `A_ab = r_a r_b sin(Phi_ab)` is twice the signed area of the triangle
`(0, z_a, z_b)`. Then Theorem E's discriminant is

    Delta / 4 = (t_ab^2 - r_a^2 r_b^2)(t_bc^2 - r_b^2 r_c^2) = A_ab^2 A_bc^2,

so `sqrt(Delta)/2 = |A_ab A_bc|` and

> **Theorem E forces `t_ac` rational if and only if `A_ab A_bc` is rational,
> which happens exactly when `sin^2(Phi_ab)` and `sin^2(Phi_bc)` have the same
> squarefree part.**

That is a condition on a field the arithmetic arc has never used. The whole
holonomy programme studies `Q(cos Phi)`; this is `Q(sin Phi)`, and the two are
different. Furthermore the plane sine-addition identity gives

    r_b^2 A_ac = A_ab t_bc + t_ab A_bc,

so when `t_ab` and `t_bc` are rational the sine field PROPAGATES: if `A_ab` and
`A_bc` lie in `Q(sqrt d)` then so does `A_ac`.

**Corollary, if the sine field is a class invariant.** Every `A` in the class
lies in `Q(sqrt d)` for one `d`, every Theorem E discriminant is automatically a
rational square, and rationality propagates through the whole class with no
case-by-case squareness check. That would replace a checked condition with a
single invariant per class.

## Predictions, committed before the run

| id | prediction | expectation |
|---|---|---|
| P-F-1 | `Delta/4 = A_ab^2 A_bc^2` exactly, on every chained triple | PASS |
| P-F-2 | the sine-addition identity `r_b^2 A_ac = A_ab t_bc + t_ab A_bc` holds up to orientation sign on every triple | PASS |
| P-F-3 | the squarefree part of `A_ab^2` at rational-`t` pairs is a recorded `sin_radicand` of `holonomy_fields.json` | PASS |
| P-F-4 | the sine field is a CLASS INVARIANT: every pair of a class has `A_ab` in `sqrt(d) * Q(t_ab)` for one `d` per class | PASS |
| P-F-5 | the sine field differs from the cosine field on at least one class | PASS |

Reasoning:

- **P-F-1** is Theorem F's identity and is algebra; a failure means the exact
  arithmetic is wrong.
- **P-F-2** is the propagation mechanism. The sign is left free because `A` is a
  SIGNED area and the orientation convention of the pair ordering has not been
  fixed; only `|A|` enters the discriminant.
- **P-F-3** is the cross-check that matters most. `sin_radicands` is an
  independently generated field of `holonomy_fields.json`, produced by a
  different script for a different purpose. If my squarefree parts are those
  radicands, the new object is tied to certified data rather than invented.
- **P-F-4 is the substantive prediction.** It is tested on all 54 pairs, not
  just the 26 with rational `t`, by asking whether `A_ab / sqrt(d)` lies in
  `Q(t_ab)`. Most classes have only ONE rational-`t` pair, so a test restricted
  to those would be nearly vacuous; that weakness is why the irrational-`t`
  pairs are included.
- **P-F-5** should already be visible: (5,10) v113 loop 0 records cosine
  radicand 2 and sine radicand 94.

## What would REFUTE this

A class carrying two pairs whose `A^2` have different squarefree parts, both
nonzero. Theorem E's discriminant would then not be a square there, `t_ac`
would be irrational, and the sine field would not be a class invariant. The
17-of-18 split would then be a coincidence rather than a degenerate limit.

## Independence, per rule R2

The 17 classes span 11 transport classes. The squarefree value 7 covers 4 of
them and 47 covers 3, so the six distinct values are spread over 11 independent
observations, not 17. Any tally in the write-up reports both.

## Nonclaims

- A class invariant on the census is not a theorem. P-F-4 is a census
  measurement; the propagation identity is proved, the invariance is not.
- The sine field is not claimed to determine the holonomy field or vice versa.
- All arithmetic reads `states.jsonl` closed forms, never the numerical
  natural-orbital ledger.

# The signed Chow projection: a proved 2d-integer encoding of a GPC normal

**Status:** three theorems, proved, plus two hypotheses that are measured and
labelled as measured. Verified on 142,493 nonstructural rows across ranks 7, 8
and 9.

**Artifact:** `results/data/signed_chow_projection.json`, written by
`scripts/signed_chow_projection.py`.

This upgrades a caveat the repository has carried since the mechanism-grammar
synthesis. Singleton-incidence injectivity was recorded there as *measured on
finite populations, a canonical hash and not a theorem*. Part of it is now a
theorem, and the part that is not is separated out explicitly.

## The object

Fix `N` and `d` and let the fixed-weight slice be

```text
Omega = { chi_I : I in binom([d], N) }
```

For a primitive normal `tau`, split the slice by sign and record only the
orbital occupancy of each side:

```text
B_+ = { omega : tau . omega > 0 },    c_+ = sum over B_+ of omega
B_- = { omega : tau . omega < 0 },    c_- = sum over B_- of omega
B_0 = { omega : tau . omega = 0 }
```

`c_+` and `c_-` are vectors of `d` nonnegative integers. The pair is `2d`
integers no matter how large `binom(d, N)` is: at `(3,40)` the slice has 9,880
weights and the shadow is 80 integers, each at most `binom(39,2) = 741`.

`c_+` is exactly the `singleton_signature` of the rank-7 to rank-9 audits, so
everything below applies directly to the numbers already published there.

## T1. Determination

**Theorem.** Let `F` be any family of fixed-weight vectors with
`sum(F) = c_+(tau)`. Then `F = B_+(tau)`.

*Proof.* Every `omega` has exactly `N` ones, so `1 . sum(F) = N |F|`, and
matching shadows force `|F| = |B_+(tau)|`. From `sum(F) = c_+(tau)`,

```text
0 = tau . ( c_+(tau) - sum(F) )
  = sum over B_+(tau) \ F of tau.omega  -  sum over F \ B_+(tau) of tau.omega
```

Every term of the first sum is strictly positive, by definition of `B_+(tau)`.
Every term of the second sum is at most zero, because those weights are not in
`B_+(tau)`, so each contributes a nonnegative amount after the minus sign. A
sum of nonnegative terms vanishes only if all vanish, and the first group has
no zero terms, so `B_+(tau) \ F` is empty. Hence `B_+(tau)` is contained in
`F`, and equal cardinality gives equality. QED

Stated for two threshold families this is the fixed-weight analogue of
Chow-parameter uniqueness for linear threshold functions. The proof never uses
any structure on `F`, so the stronger statement is free: `c_+` pins `B_+`
among **all** subfamilies of the slice, not merely among threshold ones.

## T2. Monotonicity

**Theorem.** If `tau_i >= tau_j` then `c_+(tau)_i >= c_+(tau)_j`.

*Proof.* Subsets containing both `i` and `j`, or neither, contribute equally to
both coordinates. Map each `I` in `B_+(tau)` with `j in I`, `i not in I` to
`I' = I \ {j} u {i}`. Then `tau . I' = tau . I + (tau_i - tau_j) >= tau . I > 0`,
so `I'` lies in `B_+(tau)`, contains `i` and not `j`. The map is injective, so
the `i` side is at least as large. QED

So the sorted order of `c_+` recovers the chamber order of `tau` up to ties.
**The shadow already carries the placement**, which is why it is a plausible
state variable for a placement oracle and not only a deduplication key.

## T3. Reconstruction

**Theorem.** `B_0` is determined by the pair `(c_+, c_-)`, and if `span(B_0)`
has dimension `d-1` then `tau` is the primitive generator of its annihilator.

*Proof.* T1 applied to each side gives `B_+` and `B_-` from `c_+` and `c_-`, so
`B_0 = Omega \ (B_+ u B_-)` is determined. Every `omega` in `B_0` satisfies
`tau . omega = 0`, so `span(B_0)` is contained in the hyperplane `tau^perp`,
which has dimension `d-1`. If the span attains `d-1` it equals `tau^perp`, and
the annihilator of a hyperplane is the line through `tau`. QED

The admissibility condition the generator already enforces is that the zero
weights have affine rank `d-2`. Since every `omega` has coordinate sum `N`,
which is nonzero, the affine hull misses the origin and the linear span has
dimension `(d-2) + 1 = d-1`. So admissibility is exactly the hypothesis of T3.

**Corollary.** Every admissible GPC candidate normal is determined by `2d`
integers.

## What is measured rather than proved

Two things, kept separate on purpose.

**The span hypothesis.** T3 is conditional. That `span(B_0)` actually attains
dimension `d-1` is verified on every nonstructural survivor of every system
audited here, and it is not proved in general.

**`c_+` alone.** T1 says `c_+` determines `B_+`. It does **not** say `c_+`
determines `tau`: two normals could share a positive family and differ on how
the rest splits between `B_-` and `B_0`. That `c_+` alone separates normals is
a measurement, with zero collisions on 142,493 rows, and the earlier audits'
injectivity claim sits at exactly this conjectural level. The signed pair, by
contrast, determines the normal by theorem.

**Decoding.** T1 is an existence-and-uniqueness statement, not an algorithm.
Recovering `B_+` from a bare `c_+` is the Chow-parameters inversion problem and
is not solved here. This is the one place where "encoded by `2d` integers" can
mislead: the encoding is information-theoretically complete, while inversion is
a separate and nontrivial problem. T2 is the practical lever, since it fixes
the orbital order before any search begins. Certificate *verification* needs no
inversion at all: given a claimed `tau` and a shadow, one pass over the slice
recomputes the pair, and T3 turns a match into uniqueness.

## Verification

Every nonstructural trace survivor at ranks 7, 8 and 9:

| System | Nonstructural rows | Reconstructed exactly | Span is `d-1` | Monotone | Signed collisions | `c_+` collisions |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | 547 | 547 | 547 | 547 | 0 | 0 |
| `(3,8)` | 6,605 | 6,605 | 6,605 | 6,605 | 0 | 0 |
| `(3,9)` | 135,341 | 135,341 | 135,341 | 135,341 | 0 | 0 |
| total | **142,493** | **142,493** | **142,493** | **142,493** | **0** | **0** |

Reconstruction is exact rational nullspace computation, not a numerical fit:
the recovered primitive integer vector equals the original `tau` on every row.

## Why this matters for the constructor

The rank-9 audit visits 135,343 survivors and stores, per row, a normal of `d`
integers plus derived level sets. The shadow replaces the `binom(d,N)`-sized
sign pattern with `2d` integers that provably carry the same information, and
by T2 they carry the placement order too. That makes the shadow usable as:

- the primary deduplication key, with collisions now a theorem-level question
  rather than a hoped-for property;
- a compact certificate format, verifiable in one pass over the slice;
- a candidate state variable for the placement descent, since T2 shows it is
  order-aware.

It does **not** by itself give a separation oracle, and nothing here bounds the
number of realizable shadows. It removes one specific worry, that a candidate
fundamentally needs a `binom(d,N)`-bit object to be identified.

## Verification commands

```sh
uv run scripts/signed_chow_projection.py --ranks 7 8
uv run scripts/signed_chow_projection.py --ranks 9 --checkpoint-dir .cache/orbits
uv run pytest tests/test_signed_chow_projection.py
```

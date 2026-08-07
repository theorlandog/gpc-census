# The tangent determinant, classified by its matchings

Status: **AUDITED at ranks 7 and 8**, with one measured consequence already
landed in the pipeline. The determinant's permutation terms are exactly the
perfect matchings of its support graph, and classifying them says how much of
the residual algebra is genuinely needed.

Certifying artifacts: `ressayre.classify_tangent_determinant`,
`ressayre._sparse_determinant_matching_counts`,
`tests/test_determinant_matching_audit.py`.

## The classification

Every nonzero cell of the tangent matrix is a signed on-weight variable, so each
perfect matching carries a sign and a monomial. Keeping the UNSIGNED matching
multiplicity alongside the signed coefficient separates cases that the signed
sum alone conflates:

| class | meaning | certificate |
|---|---|---|
| `structural_zero` | no perfect matching | Hall violator, no algebra |
| `unique_matching` | exactly one matching | the determinant IS one signed monomial |
| `unique_exposed_monomial` | some monomial has exactly one matching | its coefficient is `+-1` and cannot cancel |
| `nonzero_after_collision` | every monomial is multiply reached, sum survives | needs the signed sum |
| `exact_cancellation` | signed sum vanishes | the genuine algebraic case |

The middle two matter because they prove nonvanishing **without summing
anything**: a monomial reached by a single matching has coefficient `+-1`, so no
cancellation is possible, and no evaluation or large-integer determinant is
required.

## Measured at rank 7

| class | count |
|---|---|
| `structural_zero` | 474 |
| `exact_cancellation` | 25 |
| `unique_matching` | 16 |
| `unique_exposed_monomial` | 24 |
| `nonzero_after_collision` | 10 |
| total | **549** |

The 474 + 25 = 499 zeros and the 16 + 24 + 10 = 50 nonzeros reproduce the
published partition exactly.

**40 of the 50 nonzero determinants, 80 percent, are certified by a unique
exposed monomial.** Only 10 need the signed sum at all. The largest matching
count seen is 30,528, so this is not a statement about tiny cases.

## The consequence that was worth more than the audit

Profiling the orbit-native constructor at rank 8 put 48.6 of its 54 seconds
inside screening, and neither enumeration (4.8 s) nor candidate reconstruction
(0.7 s) accounted for it. The cause turned out to be structural and simple.

`assess_candidate_by_evaluation` searches for a nonzero evaluation by computing
up to 32 exact Bareiss determinants at deterministic trial points. On a row
whose determinant is IDENTICALLY zero every one of them returns zero. That is
not the rare case: **6,414 of the 6,607 rank-8 trace survivors are structurally
zero**, so the loop was spending 32 exact determinants apiece to reach a
foregone conclusion, and only afterwards falling through to the Hall test that
settles it immediately.

Running the Hall test BEFORE the trial points proves the search futile instead
of performing it.

| (3,8) | seconds |
|---|---|
| materialising reference | 251.9 |
| orbit-native | 66.0 |
| orbit-native, futile search proved rather than run | **12.4** |

**20x against the reference**, with the partition unchanged: 137 certified,
6,469 determinant zeros, 1 structural, 0 unresolved.

**The recorded verdict is byte identical, including the attempt count, and that
is deliberate.** The field records that the budget yields no nonzero evaluation,
and a Hall violator establishes exactly that statement for every trial point at
once. It is the same claim by a shorter route, so the rank-7 reference manifest
still regenerates byte for byte against the stored artifact.

## What the audit does not yet license

The `unique_exposed_monomial` class is currently DETECTED by the full expansion,
which is the thing it is meant to avoid. Turning it into a cheap certifier needs
the separating-weight construction: a matching `P` and integer variable weights
whose minimum over perfect matchings is attained only at `P`, checkable with one
min-weight matching and one second-best bound. That is not implemented, and
until it is, the 80 percent figure describes what is CERTIFIABLE this way, not
what is currently certified this way.

The audit also does not explain the 25 exact cancellations. Whether they share a
common alternating-cycle mechanism is open, and it is the natural next question:
the symmetric difference of two matchings is a union of alternating cycles, and
two matchings collide exactly when their cycles have zero net variable-label
count.

## Reproducing

```sh
uv run pytest tests/test_determinant_matching_audit.py
```

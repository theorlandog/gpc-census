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

## Measured at rank 8

| class | count |
|---|---|
| `structural_zero` | 6,414 |
| `exact_cancellation` | 55 |
| `unique_matching` | 22 |
| `unique_exposed_monomial` | 94 |
| `nonzero_after_collision` | 22 |
| total | **6,607** |

193 Hall-feasible, splitting 55 zero and 138 nonzero, which is the predicted
split exactly.

**The share needing no summation RISES with rank: 80.0 percent at rank 7, 84.1
percent at rank 8.** The largest matching count seen is 9,538,560 at order 22,
so this is emphatically not a statement about tiny cases: the audit itself took
1,167 seconds at rank 8, which is the cost the certificate below avoids.

## The certificate, in polynomial time

The audit above FOUND the cheap cases with the exponential expansion they exist
to avoid, so on its own it proved only that a cheap route should exist. It does.

Give each on-weight variable an integer weight, so every perfect matching
acquires the weight of its monomial. If one matching is strictly lighter than
every other, no other matching carries its monomial, its coefficient is `+-1`,
and the determinant cannot vanish. The certificate is a weight vector, a
matching and the second-best weight; a reader replays it with one assignment
solve per matched edge.

| rank 7 | |
|---|---|
| nonzero determinants | 50 |
| have a singly-reached monomial (audit upper bound) | 40 |
| **certified by the polynomial-time certificate** | **36** |
| certificate time vs expansion time | 0.05 s vs 0.26 s |

The gap of 4 is the choice of weight vectors, not a limit of the method: those
determinants do have a singly-reached monomial, and no vector among the four
tried exposed it. Pinning both numbers keeps that distinction visible.

**Safety is the direction that matters and it is tested.** A cheap test may say
NO when the answer is yes; it may never say YES when the answer is no, because
that would certify a facet that is not there. Every certificate issued at rank 7
is checked against the exact expansion, none was issued for a vanishing
determinant, every one replays, and tampered certificates are rejected.

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

## What this does not license

The certificate is not wired into the production screening path. It is built and
verified, and the pipeline still reaches nonzero determinants by evaluation,
which already succeeds cheaply on them; the futile-search fix above removed the
cost that mattered. Wiring it in is worth doing when the determinant order grows
past what evaluation handles comfortably, not before.

Returning `None` proves nothing. The determinant may still be nonzero and the
caller must fall back, which is why the certifier is only ever allowed to say
YES.

The audit does not explain the 55 exact cancellations. Whether they share a
common alternating-cycle mechanism is open, and it is the natural next question:
the symmetric difference of two matchings is a union of alternating cycles, and
two matchings collide exactly when their cycles have zero net variable-label
count.

## Reproducing

```sh
uv run pytest tests/test_determinant_matching_audit.py
```

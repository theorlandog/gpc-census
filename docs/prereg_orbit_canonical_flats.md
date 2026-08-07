# Pre-registration: the rank-8 orbit count

Standing rules in `docs/prereg_template.md`. Committed with
`scripts/orbit_canonical_flats.py` BEFORE the rank-8 run (rule R1). The rank-6
and rank-7 orbit counts were computed before this file existed and are recorded
below as ALREADY MEASURED, not as predictions, per rule R4.

## What is already established, and what corrects the plan

**(3,7) candidate exhaustiveness is CLOSED**, and `AGENT_PLAN.md` T2 step 1 is
stale. `docs/fixed_n3_generator_methodology.md` records that the closed-flat
enumerator gives 5,341 hyperplanes in place of 1,623,160 candidate bases,
10,682 oriented rows, screening resolving all of them with none unresolved
(10,133 trace failures, 499 exact determinant zeros, 1 structural Pauli row, 49
nonstructural Ressayre rows, summing exactly to 10,682), and exact reduction
returning the four published facets. T2 step 1 currently says "candidate
exhaustiveness is bypassed, not proved", which is true of the gate-1 inner/outer
closure but NOT of the flat enumerator, and reading only T2 makes rank 7 look
open when it is not. The next gate is (3,8).

**Already measured, before this pre-registration:**

| rank | hyperplanes | `S_d` orbits | reduction |
|---|---|---|---|
| 6 | 362 | 8 | 45.2x |
| 7 | 5,341 | 19 | 281.1x |

with `images_outside_set` zero at both, so the action is closed. A first attempt
gave 15 and 35 because the permuted image was not re-normalized to the stored
sign convention; that is recorded in the script and is the trap this measurement
most easily falls into.

## Predictions, committed before the rank-8 run

| id | prediction | expectation |
|---|---|---|
| P-O-1 | the rank-8 flat lattice reproduces `[1, 56, 1540, 21420, 147630, 467082, 565208, 166420]` | PASS |
| P-O-2 | the rank-8 `S_d` action is closed, `images_outside_set == 0` | PASS |
| P-O-3 | the rank-8 orbit count lies in **[30, 70]** | PASS |
| P-O-4 | the rank-8 orbit count lies in [40, 60], the narrower band proposed externally | PASS |
| P-O-5 | the rank-8 reduction factor exceeds 2,000x | PASS |

Reasoning, so a hit is not read as a lucky guess:

- **P-O-1** is a replication. The count came from an independent numpy
  reimplementation over `GF(2^31 - 1)`; the repo's own
  `stage1_reference_series.json` already carries the same vector, so this run
  makes it three independent derivations rather than a quoted number.
- **P-O-3 and P-O-4.** Orbits went 8 to 19, a ratio of 2.4. The same ratio puts
  rank 8 near 45. The wide band is the honest one and the narrow band is the
  externally proposed one; both are registered so that hitting the wide band
  while missing the narrow one is visible rather than glossed.
- **P-O-5** follows from P-O-3: 166,420 over at most 70 is at least 2,377.

## What would REFUTE the programme

An orbit count that grows like the raw count, say above 500 at rank 8. The
canonical-augmentation strategy would then buy a constant factor rather than a
change of complexity class, and flats-only enumeration would still die between
ranks 9 and 10.

## Nonclaims

- Orbit counts at three ranks are not an asymptotic. Extrapolating to rank 11
  from 8, 19 and one more point is a projection, not a measurement, and is
  labelled as such wherever it appears.
- This script measures the orbit structure of the ALREADY ENUMERATED set. It is
  calibration for a canonical-augmentation enumerator, not that enumerator, and
  it does not by itself make rank 9 or 11 reachable.
- Orbit reduction being lossless is a theorem about Weyl invariance, not a
  claim that any particular canonicalization implementation is correct; the
  enumerator, when written, needs its own verification against these counts.

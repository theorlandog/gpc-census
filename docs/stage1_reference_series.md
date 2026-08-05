# The rank generator: gates G2 and G3

Status: **G1 and G2 are CLOSED.** The generator reproduces the published
fixed-`N=3` systems at ranks 6, 7 and 8 from first principles, blind. G3 is
attempted and its cost is measured. No new-rank claim is made here.

Certifying artifacts: `scripts/stage1_reference_series.py`,
`results/data/stage1_reference_series.json`,
`tests/test_stage1_reference_series.py`. Machinery:
`gpc_census.generation`. Methodology and labelling rules:
`docs/fixed_n3_generator_methodology.md`.

## Why gates, and why this one matters

`docs/GENERATOR_PLAN.md` fixes the ladder:

    G1  (3,6)             first-principles, sandwich-closed
    G2  (3,7), (3,8)      blind, must match the Altunbulak-Klyachko tables
    G3  (3,9), (3,10)     blind, must match the census H-data
    G4  (3,11)            the bracket
    G5  (3,12)+           frontier

A generator that has not reproduced the known tables is not evidence about an
unknown one. This note closes G2 and reports honestly on G3.

## What "blind" means here, checked rather than asserted

The generation path is: enumerate admissible affine flats, orient them into
candidate taus, screen each by the exact Ressayre criterion, reduce on the
ordered Pauli slice. It reads no stored constraints and no stored vertices.
The stored table is consulted only afterwards, in
`candidate_system._known_system_regression`, which runs on
`reduction.retained_rows` after they exist. A rank reporting status
`finalized_known_system_regression` therefore matched the published system
without having been shown it.

`tests/test_stage1_reference_series.py` re-runs rank 7 LIVE rather than
re-reading the artifact, because a test that only parses the emitted JSON
cannot establish blindness.

## Results

| rank | admissible hyperplanes | oriented taus | retained rows | published rows | gate | time |
|---|---|---|---|---|---|---|
| 7 | 5,341 | 10,682 | 4 | 4 | **PASS** | 4.5 s |
| 8 | 166,420 | 332,840 | 31 | 31 | **PASS** | 170 s |

Both match the published system as an unordered set of oriented rows, with
`generated_only` and `published_only` both empty. The screening partition is
exhaustive at both ranks: every candidate is either certified, exactly
rejected, or structural, with zero unresolved rows.

## The cost wall, measured

The flat enumeration dominates and the screening does not:

- screening runs at about **0.14 ms per tau**, so the 332,840 taus at rank 8
  cost well under a minute of the 170 second total;
- the enumeration grows by a factor of **31 in hyperplanes** from rank 7 to
  rank 8, and about 45 in wall time.

That growth factor is what gates G3 and G5, and extrapolating it puts rank 9
in the millions of hyperplanes.

The binding constraint has MOVED. It was memory: a first rank-9 attempt held
about 9 GB of 16 GB before being cut off. After the memory fix below it is
time. Both fixes are recorded separately because they are separate problems
and the second only became visible once the first was solved.

The methodology document's labelling applies: resource exhaustion is
**operational** and carries no mathematical meaning. It is not evidence that
the generator is wrong at rank 9, and the artifact records such a run as
`status: operational_failure` with `gate_passed: false`, never as a rejection.

## What this does and does not license

**Does.** The generator is now validated blind at three consecutive ranks
including two with published tables it did not read. The Ressayre screening,
the flat enumeration, the ordered-slice reduction and the finite-field
closure (whose modulus is checked to beat the Hadamard bound, so the modular
rank computations are the characteristic-zero ones) all reproduce known
mathematics exactly.

**Does not.** Passing G2 says nothing about candidate EXHAUSTIVENESS at a new
rank, which is the upstream obligation stated in the methodology document. A
new-rank output remains a facet system relative to the supplied certified
candidates until exhaustiveness is established separately. G3 is not closed,
so the ladder is not yet at the frontier.

## The speed fix, applied and measured

The profile below located the cost in two inner-loop functions, `_residual_mod`
and `_projective_key`, together 58 percent of the enumeration. Both are pure
fixed-width modular arithmetic over a whole block of candidate points, and
because the reduction is in reduced row-echelon form the residual is exactly

    residual = point - point[pivots] @ rref,

one matrix product for the entire candidate block of a parent rather than one
Python loop per point. The projective normalisation, the leading-entry inverse
(from a tabulated inverse array rather than a modular exponentiation) and the
grouping key are vectorised with it. The hyperplane incidence test at the end,
`<h, omega> == z` for every pair, is likewise one integer matmul.

| enumeration only | before | after |
|---|---|---|
| rank 7 | 5.6 s | **2.7 s** |
| rank 8 | 298 s | **121 s** |

End to end, including screening and the ordered-slice reduction, rank 8 went
from 423 s to **170 s**. Output is byte-identical at both ranks: the same
hyperplane count, the same hyperplane tuple, and the same flat count at every
rank.

Two things were deliberately NOT done. The residual grouping loop and the
Bareiss cofactor determinants are what remain, and both are exactness-critical
integer arithmetic where a vectorised rewrite trades a modest constant for real
risk. The scalar implementation is kept as `_extend_one_rank_scalar`, both as
the definition the vectorised path is checked against and as the fallback when
the modulus is too large to tabulate inverses.

## The memory fix, applied and measured

Profiling located the cost precisely. A stored flat took about 950 bytes, of
which about 870 was a materialised RREF matrix carried on every flat. That is
redundant twice over: the reduction is a function of the basis, and it is
needed only while a flat is being EXTENDED. The level dict now stores
`mask -> basis_indices` and reduces once per PARENT, which also removes one
reduction per CHILD, and children outnumber parents heavily at the wide middle
ranks.

| rank 8 | before | after |
|---|---|---|
| peak RSS | 1302 MB | **298 MB** |
| wall time | 343 s | 298 s (121 s after the speed fix above) |
| hyperplanes | 166,420 | 166,420 |
| flat counts by rank | 1, 56, 1540, 21420, 147630, 467082, 565208, 166420 | identical |

A 4.4x memory reduction with byte-identical output. The flat lattice counts
are pinned by `test_flat_lattice_counts_are_pinned` and by the artifact, so a
future change to the level representation cannot silently alter the
enumeration.

Rank validation is preserved and not weakened: every flat is reduced when it
is extended, and the final rank, which is never extended, is checked exactly
by the hyperplane closure test (`on_mask != state_mask`), which compares the
finite-field closure against an exactly recomputed integer hyperplane.

That moves rank 9 from "OOM on a 16 GB machine" to roughly 3 to 4 GB by the
same scaling, so the remaining obstacle there is TIME, not memory: rank 7 to
rank 8 cost a factor of about 53 in wall time.

## The time wall, profiled

Rank 9 needs roughly a factor of 50. A profile of the rank-7 enumeration
(12.3 s under cProfile) says where the time is, and therefore what micro
optimisation can and cannot buy:

| site | tottime | share |
|---|---|---|
| `_residual_mod` (1,192,695 calls) | 2.59 s | 21% |
| `_projective_key` (1,192,695 calls, cumulative 4.55 s) | 2.44 s | 20% |
| `_extend_one_rank` body | 1.51 s | 12% |
| `_rref_mod` (39,810 calls) | 0.87 s | 7% |
| `_bareiss_det` | 0.83 s | 7% |
| `pow` (modular inverse, 1,365,000 calls) | 0.47 s | 4% |

The two inner-loop functions together are 58 percent, which bounds what
vectorising them can buy at about 2.4x by Amdahl. That is exactly what was
obtained (2.46x on the enumeration), and it is not 50x.

The honest diagnosis is that the remaining cost is the SIZE OF THE FLAT
LATTICE, not inefficiency per flat. Rank 8 visits about 1.37 million closed
flats (1 + 56 + 1540 + 21420 + 147630 + 467082 + 565208 + 166420) and the work
is close to linear in that count. Two measurements pin that down:

- **the dedup absorbs a lot, and that is inherent.** Rank 8 generates
  15,632,106 child insertions for 1,369,356 unique flats, a redundancy of
  11.4x (8.5x at rank 7, so it grows). Each flat has many rank-`k-1` subflats
  and is therefore reached from many parents. Removing the redundancy would
  need a canonical-parent test, and the scan that finds a parent's residual
  classes has to happen anyway: the cost is (parents) times (points), which is
  59.8 million residuals at rank 8 whether or not the children are deduped.
- **pruning cannot help.** Every rank-`k` flat is on a chain to some
  rank-`d-1` flat, and the trace test that rejects 95 percent of candidates
  downstream is a property of a COMPLETED hyperplane `(h, z)`, so it has
  nothing to test on a partial flat.

So the three routes to G3 are, in order of honesty about effort:

1. **Accept the runtime, incrementally. IMPLEMENTED.** Rank 9 fits in memory
   now; it wants hours, not gigabytes. `enumerate_admissible_hyperplanes_by_flats`
   and `generate_fixed_n3_reference_rank` take a `checkpoint_dir`, and
   `scripts/stage1_reference_series.py` exposes `--checkpoint-dir`. Each
   completed level is written out atomically and reloaded on a later call, so
   a rank can be finished across several runs without a single long-lived
   process. Checkpoints are plain text (`hexmask:i1,i2,...`), never pickle,
   because a checkpoint is data and loading one must not be able to execute
   anything. `test_checkpoint_round_trip_is_identical` pins that a resumed
   enumeration reproduces the fresh one exactly, hyperplane tuple included.
2. **Compile the inner loop.** The residual and projective-key steps are now
   vectorised, which took the 2.4x that Amdahl allows. Going further means a
   compiled kernel for the whole level induction, including the residual
   grouping and the row reduction. This is the only route that plausibly gets
   the remaining order of magnitude.
3. **Change the enumeration.** A different exhaustiveness argument that does
   not traverse the whole flat lattice would be a genuine improvement, and
   there is no candidate on the table.

Streaming the hyperplanes into screening remains available but is now a minor
win, since the final level is no longer the memory peak. Parallelising the
screening buys nothing: it already runs at 0.14 ms per tau and is not the
bottleneck.

## Reproducing

```sh
python scripts/stage1_reference_series.py --ranks 7 8    # about 7 minutes
uv run pytest tests/test_stage1_reference_series.py

# a rank that wants hours, finished across several runs
python scripts/stage1_reference_series.py --ranks 9 --checkpoint-dir build/flats
```

The script accumulates across invocations, so a later rank can be appended
without recomputing the earlier ones, and `--checkpoint-dir` additionally
resumes an unfinished rank at the last completed flat level.

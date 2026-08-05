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
| 7 | 5,341 | 10,682 | 4 | 4 | **PASS** | 8 s |
| 8 | 166,420 | 332,840 | 31 | 31 | **PASS** | 423 s |

Both match the published system as an unordered set of oriented rows, with
`generated_only` and `published_only` both empty. The screening partition is
exhaustive at both ranks: every candidate is either certified, exactly
rejected, or structural, with zero unresolved rows.

## The cost wall, measured

The flat enumeration dominates and the screening does not:

- screening runs at about **0.14 ms per tau**, so the 332,840 taus at rank 8
  cost roughly a minute of the 423 second total;
- the enumeration grows by a factor of **31 in hyperplanes** from rank 7 to
  rank 8, and about 50 in wall time.

That single number is what gates G3 and G5. Extrapolating the same factor puts
rank 9 in the millions of hyperplanes, and the binding constraint is MEMORY
rather than time: the enumeration materialises the hyperplane list. A rank-9
attempt in this environment held about 9 GB of 16 GB before being cut off.

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
| wall time | 343 s | 298 s |
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

## The engineering the next step still needs

1. **Prune inside the enumeration.** About 95 percent of taus are
   trace-rejected by a cheap test (10,133 of 10,682 at rank 7). Applying an
   equivalent test during extension would prune the flat tree rather than the
   output list, which is the only change that attacks the time wall rather
   than the memory wall.
2. **Stream the hyperplanes into screening.** Screening is per-tau and
   stateless, so it can be fused with the final level instead of
   materialising both lists. This is now a minor win, since the final level is
   no longer the memory peak.
3. **Only then parallelise.** Screening is already embarrassingly parallel and
   runs at 0.14 ms per tau; parallelising it first would buy nothing.

## Reproducing

```sh
python scripts/stage1_reference_series.py --ranks 7 8    # about 7 minutes
uv run pytest tests/test_stage1_reference_series.py
```

The script accumulates across invocations, so a later rank can be appended
without recomputing the earlier ones.

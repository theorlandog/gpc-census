# Phase 1: the three-way birationality audit

**Status:** criterion comparison at `(3,7)`, `(3,8)` and `(3,9)`, on 424
determinant-nonzero Hall-feasible candidate rows with published facethood as
ground truth.

> ## ❗ CORRECTION (2026-08-09), read before the numbers
>
> The first version of this document had the probabilistic evidence **exactly
> backwards**. It claimed a `False` birationality verdict was weak and a `True`
> was the strong side. The pinned source says the reverse.
>
> `is_not_contracted` samples matrices over `Q[i]` and returns True the moment
> one has full column rank, breaking early. `Is_Ram_contracted` returns
> **`False` as soon as any divisor yields such a witness**, and returns `True`
> only after the search is exhausted without finding one. So:
>
> - **`False` = a witness was found.** Exact, and the direction that carries
>   evidence. The sampling only affects whether the witness is *found*, never
>   whether it is valid: a full-rank evaluation at a specific rational point is
>   an exact lower bound on the generic rank.
> - **`True` = no witness found.** Monte Carlo only.
>
> Worse, the first run left `random_deep` at the upstream default of **1**, a
> single trial. The old headline "every published facet is birational" therefore
> meant "one trial failed to find an obstruction for each facet".
>
> This version fixes the direction, sets `random_deep = 8` explicitly, and
> repeats the whole sweep across three independent seeds. The claims are now
> separated into the direction that proves and the direction that does not.

**Artifacts:** `results/data/bdr_birationality_candidates.json` (the candidate
population and the local column), `results/data/bdr_birationality_external.json`
(the upstream verdicts, recorded because they need a SageMath runtime) and
`results/data/bdr_birationality_audit.json` (the join and the scoring).
Generating script: `scripts/bdr_birationality_audit.py`. Guard test:
`tests/test_bdr_birationality_audit.py`.

## The question

Gate 4 measured a local proxy for linear triangularity. Gate 2 showed that proxy
is incomparable to BDR's own linear-triangular filter. Neither is what decides a
facet. BDR's **birationality** test is, and it had never been run here.

This audit puts all three verdicts on the same rows and scores each as a
predictor of facethood, rather than describing them.

The population is every Hall-feasible row with a nonzero Ressayre tangent
determinant: 48, 136 and 240 rows against published irredundant systems of 4, 31
and 52. Those rows have passed the trace test, Hall feasibility and a nonzero
determinant, and nothing else.

🟦 **Every published facet is inside the candidate population** at all three
systems, `published_facets_not_candidates` empty each time. So the scoring below
is not distorted by facets the enumeration never produced.

## Result

### The strong direction: 228 rows carry an actual obstruction

A `False` verdict is a found witness, so these rows are non-birational up to
exact arithmetic, not merely suspected:

| system | candidates | facets | witnessed NOT birational | of which facets |
| --- | ---: | ---: | ---: | ---: |
| `(3,7)` | 48 | 4 | 22 | **0** |
| `(3,8)` | 136 | 31 | 70 | **0** |
| `(3,9)` | 240 | 52 | 136 | **0** |
| **pooled** | **424** | **87** | **228** | **0** |

🟦 **No published facet ever produced a non-contraction witness**, across three
seeds at depth 8. Every one of the 228 witnessed rows is a non-facet. That is
the part of the result that carries evidence, and it is a genuine cross-check:
this repository's published systems were produced without BDR, and the two
pipelines share no code.

### The Monte Carlo direction: what survived the search

The remaining 196 rows returned `True` on every trial, meaning only that no
obstruction was found:

| criterion | tp | fp | fn | tn | recall on facets | precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `no_obstruction_found` | 87 | 109 | 0 | 228 | 1.000 | 0.444 |
| `bdr_linear_triangular` | 7 | 16 | 80 | 321 | 0.080 | 0.304 |
| `dm_fully_triangular` | 10 | 47 | 77 | 290 | 0.115 | 0.175 |

⚠️ **What this does and does not license.**

- **Necessity** (every facet is birational) is **empirically supported, not
  proved**. 0 of 87 facets produced a witness under 3 seeds at depth 8, with
  zero seed disagreement, which is meaningful evidence. It is not a proof,
  because absence of a witness is not absence of an obstruction.
- **Insufficiency** (109 non-facets are birational) is a **candidate claim
  only**. Those are `True` verdicts, which is the weak direction, so some could
  be false positives of the random search rather than genuinely birational rows.
  The earlier document stated this as established. It is not.

🟦 **Stability, which is why the Monte Carlo half is worth reporting at all.**
Raising the depth from 1 to 8 and sweeping three independent seeds changed
**nothing**: the per-system counts are identical to the single-trial run
(26/66/104 survivors), and **zero rows disagreed across seeds**. A 24-fold
increase in sampling moved no verdict. That does not convert Monte Carlo into
proof, but it does say the residual uncertainty is small and concentrated
somewhere sampling does not reach.

Precision on the surviving direction still climbs with rank, 0.154 to 0.470 to
0.500, while recall stays at 1.000.

### Both triangularity readings are poor facet tests

❌ Recall 0.080 for BDR's own linear-triangular filter and 0.115 for the DM
`1 x 1` test. Each misses roughly nine facets in ten. That is the quantitative
form of gate 4's refutation, and it extends it: the earlier statement was that
strict triangularity fails on most facets, and the sharper one is that neither
triangularity notion is usable as a facet criterion at all, in either direction.
Their precision is also poor, 0.304 and 0.175, so they are not even conservative
prefilters.

This is consistent with how the upstream treats its own filter:
`LinearTriangular` is not in the default filter list, and rows it cannot settle
go on to the birationality machinery. What this audit adds is the number
attached to that design choice.

## Caveat that travels with every BDR verdict here

`Is_Ram_contracted` ran its **probabilistic** arm, the upstream default. Its
exact arm does not run in this environment, for the reason measured and bisected
in `results/data/bdr_external_benchmark.json`: libSingular fails above roughly
48 ring variables and every system here needs 70 or more.

The asymmetry is the one stated in the correction above and it cuts **against**
the headline, not for it. The 228 witnessed rows are exact. Everything built on
`True` verdicts, which is both the 87 facets and the 109 non-facets, is Monte
Carlo. Running the exact arm is what would settle it, and this environment
cannot.

Sampling used here: seeds `20260809`, `20260810`, `20260811` at
`random_deep = 8`, against the upstream default of 1. Both are serialized into
the artifact rather than left implicit.

## Scope

Three systems. The candidate population is already filtered by the trace test,
Hall feasibility and a nonzero tangent determinant, so these are not scores
against the full Ressayre-admissible population; they are scores against the
rows a screening pipeline actually has to decide between. Nothing here is
uniform in `d`.

## Binding, and why the join fails closed

A missing external row would otherwise be scored as `unknown` and dropped from
the recall and precision denominators, so an incomplete external run would look
statistically better than a complete one. The join therefore refuses to score
anything until it has proved, in `check_binding`:

- the external artifact records the sha256 of the exact candidates file supplied;
- the tau sets are equal in both directions, per system;
- no row errored, and every row carries a verdict and a trial per seed.

Any of those failing is a hard exit, not a downgrade.

## Verification

```sh
# 1. the deterministic candidate list, repository environment
uv run scripts/bdr_birationality_audit.py \
    --candidates-out results/data/bdr_birationality_candidates.json

# 2. the upstream verdicts, in a venv with moment_cone
#    (recipe in docs/bdr_linear_triangular_gate.md)
python scripts/bdr_birationality_audit.py \
    --external-out results/data/bdr_birationality_external.json

# 3. the join, repository environment
uv run scripts/bdr_birationality_audit.py
uv run pytest tests/test_bdr_birationality_audit.py
```

Pass 3 recomputes the local column from `tau` rather than reading it back, and
binds the external artifact by sha256.

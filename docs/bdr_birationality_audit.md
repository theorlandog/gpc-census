# Phase 1: the three-way birationality audit

**Status:** criterion comparison at `(3,7)`, `(3,8)` and `(3,9)`, on 424
determinant-nonzero Hall-feasible candidate rows with published facethood as
ground truth.

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

Pooled over 424 rows, with facethood as ground truth:

| criterion | tp | fp | fn | tn | recall on facets | precision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `bdr_birational` | 87 | 109 | **0** | 228 | **1.000** | 0.444 |
| `bdr_linear_triangular` | 7 | 16 | 80 | 321 | 0.080 | 0.304 |
| `dm_fully_triangular` | 10 | 47 | 77 | 290 | 0.115 | 0.175 |

Per system:

| system | criterion | tp | fp | fn | recall | precision |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | birational | 4 | 22 | 0 | 1.000 | 0.154 |
| | linear triangular | 1 | 5 | 3 | 0.250 | 0.167 |
| | DM `1 x 1` | 0 | 14 | 4 | 0.000 | 0.000 |
| `(3,8)` | birational | 31 | 35 | 0 | 1.000 | 0.470 |
| | linear triangular | 3 | 5 | 28 | 0.097 | 0.375 |
| | DM `1 x 1` | 5 | 15 | 26 | 0.161 | 0.250 |
| `(3,9)` | birational | 52 | 52 | 0 | 1.000 | 0.500 |
| | linear triangular | 3 | 6 | 49 | 0.058 | 0.333 |
| | DM `1 x 1` | 5 | 18 | 47 | 0.096 | 0.217 |

### Birationality is necessary, and it is not sufficient

🟦 **Zero false negatives, 87 of 87 facets, across all three systems.** Every
published facet passes BDR's birationality test. The two pipelines share no code
and this repository's published systems were produced without it, so this is an
independent cross-check of both, on the criterion that actually decides.

❌ **109 non-facets also pass**, so precision is 0.444 pooled. Birationality
does not separate a facet from a valid-but-redundant row, because redundancy is
not a birationality property. Deciding which of the birational rows survive is
Farkas reduction, which is exactly the irreducible work
`docs/bdr_constructor_comparison.md` measured: the step no external generator
can remove without becoming circular.

Precision climbs with rank, 0.154 to 0.470 to 0.500, while recall stays pinned
at 1.000. So the criterion gets more useful as systems grow, and never at the
cost of a missed facet in this scope.

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

A probabilistic birationality test can in principle reject a valid row, so a
`False` verdict is weaker evidence than a `True` one. That asymmetry cuts in the
favourable direction for the headline: the zero-false-negative result is built
entirely from `True` verdicts on facets, which is the strong side. The 109 false
positives are also `True` verdicts, so they are not an artefact of randomness
either. Only the `False` column, which carries the negative results, is the
weaker one.

## Scope

Three systems. The candidate population is already filtered by the trace test,
Hall feasibility and a nonzero tangent determinant, so these are not scores
against the full Ressayre-admissible population; they are scores against the
rows a screening pipeline actually has to decide between. Nothing here is
uniform in `d`.

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

# Gate 2: the BDR backend was fetched, run, and compared row by row

**Status:** cross-implementation comparison at `(3,7)`, `(4,7)`, `(3,8)` and
`(4,8)`. Two results, one of them a correction to this repository.

**Artifacts:** `results/data/bdr_external_linear_triangular.json` (the upstream
verdicts, recorded because they need a SageMath runtime to recompute) and
`results/data/bdr_linear_triangular_gate.json` (the join). Generating script:
`scripts/bdr_linear_triangular_gate.py`. Guard test:
`tests/test_bdr_linear_triangular_gate.py`.

## The correction that comes first

`docs/bdr_constructor_comparison.md` said, and this repository then repeated
for several sessions, that the pinned Bulois-Denis-Ressayre implementation
"could not be fetched or run in this environment", on four blockers.

That conclusion is **false in this environment**, and it was never tested. The
individual blockers are mostly accurate; the inference from them is not.

| stated blocker | actually true? | load bearing? |
| --- | --- | --- |
| `plmlab.math.cnrs.fr` refused by the egress proxy | yes | **no**, the code is mirrored at `github.com/ea-icj/moment_cone` |
| `ea-icj.github.io`, `arxiv.org`, `cnrs.hal.science` refused | yes | no, none of them is needed to run the code |
| `moment-cone` is not on PyPI | yes, 404 | no, `pip install .` from the clone works |
| SageMath must never be a package runtime dependency | yes, and unchanged | no, the upstream ships a `passagemath` extra that installs into a throwaway venv |

Direct HTTPS to `github.com` does return 403 through the egress proxy, which is
probably how the wrong conclusion survived. The session's repository-attachment
path serves the mirror anyway. The clone lands on
`7ab347d9a7837d68d92bde9fac74606913f395ca`, which is the exact revision the
comparison document had already pinned, and its `LICENSE` is **MIT**, replacing
the `not_verified` record.

Reproducing the external half, outside the repository environment:

```sh
python3 -m venv /tmp/mcenv
/tmp/mcenv/bin/pip install lrcalc sympy passagemath-combinat passagemath-modules \
    passagemath-polyhedra passagemath-repl passagemath-singular
git clone https://github.com/ea-icj/moment_cone /tmp/moment_cone
git -C /tmp/moment_cone checkout 7ab347d9a7837d68d92bde9fac74606913f395ca
/tmp/mcenv/bin/pip install /tmp/moment_cone
/tmp/mcenv/bin/python scripts/bdr_linear_triangular_gate.py \
    --external-out results/data/bdr_external_linear_triangular.json
```

None of that touches `pyproject.toml`, `uv.lock` or the package's runtime
dependencies. The minimum-dependency rule is intact; what changes is that a
claim of impossibility has been replaced by a measurement.

`docs/bdr_constructor_comparison.md` stays labelled **INCOMPLETE**, because it
is: it measures wall clock and stage counts, and no BDR timing was ever taken.
Its recorded *reason* is what was wrong, and the artifact now carries the
correction alongside the original benchmark.

## Result 1: the published systems are externally confirmed

The repository's published irredundant systems have never before been checked
against an independently computed source at these ranks. They are now.

| System | published rows | upstream reference rows | in both | published only | upstream only |
| --- | ---: | ---: | ---: | ---: | --- |
| `(3,7)` | 4 | 4 | 4 | 0 | none |
| `(4,7)` | 4 | 4 | 4 | 0 | none |
| `(3,8)` | 31 | 32 | 31 | 0 | `(0,0,0,0,0,0,0,-1)` |
| `(4,8)` | 15 | 16 | 15 | 0 | `(0,0,0,0,0,0,0,-1)` |

The sets are equal. The single upstream extra at rank 8 is the structural row
`lambda_d >= 0`, which upstream keeps in the same list and this repository
stores separately from the GPC rows; the convention map in
`docs/bdr_constructor_comparison.md` already named it.

This is a genuinely independent confirmation. The two pipelines share no code,
enumerate different objects (1-PS against closed flats of the augmented
exterior weights) and certify by different criteria.

## Result 2: BDR's filter and the DM test are incomparable

`docs/levi_triangularity_gate.md` asserted, without running anything, that
BDR's recursive linear-triangular filter "is a strictly stronger criterion than
asking whether every DM block of the tangent matrix is `1 x 1`". Strictly
stronger would mean every BDR-triangular row is DM-triangular. It is not.

Pooled over the 54 rows that are in both systems and analysable on the local
side:

| | BDR triangular | BDR not triangular |
| --- | ---: | ---: |
| **DM all blocks `1 x 1`** | 4 | 6 |
| **DM has a larger block** | 6 | 38 |

Both off-diagonal cells are nonempty, so neither criterion implies the other.
The claim is **REFUTED**.

Per system:

| System | joined rows | BDR triangular | DM triangular | BDR only | DM only |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(3,7)` | 4 | 1 | 0 | 1 | 0 |
| `(4,7)` | 4 | 1 | 0 | 1 | 0 |
| `(3,8)` | 31 | 3 | 5 | 2 | 4 |
| `(4,8)` | 15 | 5 | 5 | 2 | 2 |

Named witnesses, one per direction:

```text
BDR triangular, DM single block
    tau = (-2, 1, 1, 1, 1, -2, -2)          at (3,7), order 4, blocks [4]
    tau = (2, -1, -1, -1, -1, -1, -1, 2)    at (3,8), order 6, blocks [6]

DM fully triangular, BDR rejects
    tau = (2, -1, -1, 2, -4, -1, -1, -1)    at (3,8), order 5,  blocks all 1
    tau = (1, 2, -3, -1, -2, -1, 0, 1)      at (3,8), order 14, blocks all 1
```

### Why they differ

They are not two readings of one object. They are two objects.

The local test factorizes the square Ressayre tangent matrix, whose rows are
the positive weights of `tau` and whose columns are the inversion roots, and
asks whether that matrix is triangular up to permutation.

BDR's `is_linear_triangular` works on the weight graph of the whole
representation. It builds `MG` from the inversion roots on all `dim V` weights,
forms the path-count matrix `(I - MG)^{-1} - I - MG`, and calls the equation at
a positive weight *linear* when every path count into a **non-positive** weight
is at most one. It then removes exactly the roots those linear equations
determine and recurses, requiring at each level that the number of linear
equations equals the number of roots they pin. A row is accepted only if the
recursion strips every root.

So BDR can accept a row whose tangent matrix is one irreducible block, because
its linearity test is about path multiplicity in the ambient weight graph, not
about the support pattern of the tangent matrix. And it can reject a row whose
tangent matrix is perfectly triangular, because a level of its recursion can
produce fewer linear equations than roots and stop.

### What this does and does not change

**Unaffected:** gate 4's refutation. `docs/levi_triangularity_gate.md` measures
whether the *tangent system solves by back substitution*, which is the property
the Levi Triangularity Conjecture would need in order to remove the generic
Groebner calculation. That measurement is about the tangent matrix and stands
exactly as recorded: 0 of 4, 5 of 31 and 5 of 52 published rows at ranks 7, 8
and 9.

**Corrected:** the one sentence in that document relating the DM test to the
upstream filter. The DM `1 x 1` test is not a conservative under-approximation
of BDR's filter and may not be used as one in either direction.

**New and worth naming:** BDR's own filter validates only 1 of 4, 1 of 4, 3 of
31 and 5 of 15 published rows. It is a partial validator upstream too, by
design, and the rows it cannot settle go on to the birationality machinery. So
the fraction of facets that any linear-triangularity argument can reach is
small on both readings, which is the substantive agreement underneath the
disagreement about which small set it is.

## Scope

Nothing here is uniform in `d`.

> **SCOPE CORRECTION (2026-08-09), one day later.** The first version of this
> section said four systems at `d <= 8`, because `(3,9)` "would mean running the
> full BDR pipeline rather than one filter". That is true and it is not an
> obstacle: the full pipeline runs `(3,9)` in about 20 seconds. The gate is
> being extended to every census system the upstream algorithm accepts, with
> rows for the unshipped ones computed by running `moment_cone(V)` end to end.
> The section below records the pipeline and the one system it declines; the
> comparison tables above still cover only the four shipped systems until the
> extended artifact lands.

## Running the whole pipeline, and what it costs

The upstream pipeline is seven stages in `MomentConeStep.apply`:
`GeneralStabilizerDimensionCheck`, `TauCandidatesStep` (`find_1PS`),
`SubModuleConditionStep`, `StabilizerConditionStep`,
`InequalityCandidatesStep` (`List_Inv_Ws_Mod`), `TPiPreComputationStep`, then
the chosen filters, then `ExportStep`. The filter list is a menu, not a fixed
chain: `PiDominancy`, `LinearTriangular`, `BKRCondition`, `Birationality`,
`Grobner`, defaulting to `PiDominancy` then `Birationality`. Note that
`LinearTriangular`, the subject of this gate, is **off by default**, which is
what an optional partial validator should be.

At `(3,7)` the stage-by-stage attrition is 56 tau candidates, 33 after the
submodule condition, 7 after the stabilizer condition, 27 inequality candidates,
18 after `PiDominancy`, 4 after `Birationality`.

Birationality is `Is_Ram_contracted` in `ramification.py`. It walks the strong
Bruhat covering relations of `w`, builds each boundary divisor and then `R_0`,
and asks for full column rank of a matrix sliced out of the `T_Pi` 3-tensor.
Both `--ram_schub_method` and `--ram0_method` default to **probabilistic**,
retrying `random_deep` times over `Q[i]` and breaking on full rank; `symbolic`
and `symbolic_int` are available. The recorded runs here use a fixed seed so the
computed systems are reproducible, and the set comparison against the published
system is what actually checks them.

Worth naming because it changes how the two methods compare: BDR's birationality
test is a **rank** condition on a rectangular matrix, not a determinant of a
square one. That is a different object from the square Ressayre tangent matrix
this repository factorizes, which is part of why the two triangularity criteria
came out incomparable above.

## `(3,6)` is out of scope for the upstream algorithm

Borland-Dennis is the one census system BDR declines, and it says so itself:

```text
ValueError: The general stabilizer of K in V is too big.
Namely of dimension 2. The program does not work in this case.
```

That is `GeneralStabilizerDimensionCheck` enforcing BDR Condition C, generic
isotropy of dimension one. It is a scope limit of their method, not a
disagreement between the two pipelines, and it is recorded in the artifact under
`out_of_scope_upstream` so the hole in the cross-validation is documented rather
than silent. `(3,6)` is covered here by the exhaustive reference gate in
`docs/stage1_ressayre_3_6.md`.

## Verification

```sh
uv run scripts/bdr_linear_triangular_gate.py
uv run pytest tests/test_bdr_linear_triangular_gate.py
```

The first command recomputes the local column from scratch and rejoins it to
the recorded external verdicts, so it runs on a machine with no Sage. Only
`--external-out` needs the external venv.

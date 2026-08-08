# PRE-REGISTRATION: Bulois-Denis-Ressayre constructor comparison (2026-08-08)

Written BEFORE the measuring script was run. Standing rules:
`docs/prereg_template.md`. Certifying artifacts:
`scripts/bdr_constructor_comparison.py`,
`results/data/bdr_constructor_comparison.json`,
`tests/test_bdr_constructor_comparison.py`.

## What is being tested

Whether the Bulois-Denis-Ressayre (BDR) fixed-representation moment-cone
algorithm can remove exact work from this repository's fixed-`N=3` constructor.

The question is NOT whether BDR is a good algorithm. It is whether adopting its
candidate-generation or birationality machinery reduces the work this
repository must still do, given the standing acceptance rule that every row
entering a published system carries a locally replayed exact certificate.

## External object under comparison, pinned before the run

```text
paper            An Algorithm to compute the Kronecker cone and other moment cones
authors          Michael Bulois (ICJ, AGL), Roland Denis (INSMI-CNRS, ICJ),
                 Nicolas Ressayre (ICJ, AGL)
preprint         arXiv:2505.08812 (math.AG); HAL hal-05044759, version v3
code             ea-icj/moment_cone, hosted at plmlab.math.cnrs.fr,
                 landing page https://ea-icj.github.io/
revision         7ab347d9a7837d68d92bde9fac74606913f395ca
                 (the revision already pinned by scripts/bdr_fixed_n3_backend.py)
runtime          Python-Sage (SageMath)
license          NOT VERIFIED, see the availability record below
relevant case    the fermionic cone, GL_d(C) acting on wedge^r C^d, here r = 3
```

## Availability record, established before the measurement

These probes were run before any benchmark and are feasibility facts, not
results of the comparison:

1. `plmlab.math.cnrs.fr` returns a 403 CONNECT denial from this session's
   egress proxy, recorded in the proxy status endpoint as a policy denial.
2. `ea-icj.github.io`, `arxiv.org`, `cnrs.hal.science` and
   `www.semanticscholar.org` are refused by the same egress policy.
3. `moment-cone` is not published on PyPI: `pip download moment-cone` reports
   no matching distribution.
4. SageMath is not installed in this environment, and adding it would violate
   the minimum-dependency house rule for a component that must never become a
   package runtime dependency.

Therefore the external implementation cannot be executed, and it cannot even be
fetched, in this environment. The consequences for scoring are fixed below
BEFORE the run rather than negotiated afterwards.

## Treatment of unavailable or incomplete external code

Declared in advance:

1. No BDR stage count, runtime, memory figure, or candidate count will be
   written down unless it was produced by an actual run of the pinned revision.
   An unrunnable backend yields the label INCOMPLETE, never an estimate.
2. The algorithm correspondence is still deliverable, because the repository
   already contains an audit of the pinned call graph
   (`docs/fixed_n3_generator_methodology.md`, the "Pinned BDR discovery
   backend" section) written in an earlier session that did have the code
   installed. That audit is quoted as a secondary source and labeled as such.
   It is not re-derived here and must not be presented as a fresh source read.
3. The published fermionic systems for `(3,6)`, `(3,7)` and `(3,8)` are the
   objects any correct external run must reproduce. They are available locally
   and are recertified here from scratch.
4. In place of a BDR-versus-local benchmark, which cannot be run, the
   measurement is an ORACLE CEILING defined below. This is a deliberate
   substitution and is stated before the run, not discovered afterwards.

## The oracle ceiling, and why it settles the question without running BDR

The acceptance criteria for a production integration require that every row
emitted by an external tool is independently certified locally. So the best
case for ANY external candidate generator, BDR or otherwise, is that it hands
this repository exactly the minimal facet rows for free, with zero local
candidate generation and zero local screening of rejected rows.

Define, per system `(3,d)`:

- `T_full(d)`: wall clock of the full local orbit-native construction, stage
  attributed into enumeration, screening, and composition.
- `T_oracle(d)`: wall clock of local screening plus composition when the
  published minimal rows are supplied directly as the candidate set.
- `Ceiling(d) = T_full(d) / T_oracle(d)`.

`Ceiling(d)` is an upper bound on the speedup available from replacing local
candidate generation with any external generator whose rows are still locally
certified. A perfect oracle is strictly stronger than BDR, which must at
minimum emit some rows this repository then rejects. So a measured `Ceiling`
bounds BDR from above without BDR being run.

This bound is only a bound. A large ceiling does not show BDR attains it, and a
small ceiling does show BDR cannot beat it. The asymmetry is the point:
a small ceiling is a decisive negative result, a large ceiling is
UNINFORMATIVE about BDR itself and only relocates the question.

## Overlapping systems

```text
(3,6)   published Borland-Dennis system, 1 nonstructural row
(3,7)   published AK2008 Table 2 system, 4 nonstructural rows
(3,8)   published AK2008 Table 4 system, 31 nonstructural rows
```

`(3,9)` is excluded: the local full construction is not budgeted here, and no
external path can be run to compare it against. `(4,8)` is excluded for the
same reason. Excluding them before the run prevents a post hoc choice of the
system that flatters the conclusion.

## Exact convention maps, fixed before the run

The bridge is `src/gpc_census/generation/bdr.py`, already in the repository.

```text
BDR homogeneous form        tau . lambda <= 0, tau in Z^d
local affine form           a . lambda <= b on the slice sum(lambda) = N
forward map                 tau_i = N * a_i - b, then divide by content
inverse map                 constraint_from_tau, deterministic nonnegative
                            affine representative with min coefficient zero
orientation                 positive rescaling only; a sign change is a
                            DIFFERENT row and is never used to match
primitive normalization     divide by the positive gcd of entries; content
                            zero is an error, not a silent pass
particle number             N = 3 throughout; the trace slice is sum = 3
Weyl chamber                local ordered slice lambda_1 >= ... >= lambda_d
                            with lambda_1 <= 1 and lambda_d >= 0
root orientation            negative type-A roots e_j - e_i with i < j, for
                            which H_j - H_i < 0
Ressayre H form             S = sum(tau), H_i = S - d * tau_i, z = 3 * S, so
                            H . lambda >= z equals tau . lambda <= 0
structural row              (0,...,0,-1) is lambda_d >= 0 and is separated
                            from the GPC rows, never counted with them
```

Round trips are compared after primitive positive normalization only. Any
comparison that would need a sign flip to succeed is recorded as a mismatch.

## Predictions

Scored by the artifact alone.

P1. Every published nonstructural row for `(3,6)`, `(3,7)` and `(3,8)`, when
converted to primitive `tau` and screened locally, receives final status
`certified`. None is `trace_rejected`, `symbolic_zero_rejected`,
`inadmissible`, or `evaluation_unresolved`. Expectation: PASS, because these
rows are the true facets and the local Ressayre criterion certified them
during the existing rank-7 and rank-8 gates. A failure would be a real finding
about what an externally supplied minimal system does to this pipeline.

P2. The oracle path finalizes and reproduces the published system exactly at
all three ranks: `known_system_regression` available, oriented `tau` sets
equal in both directions, `generated_only` and `published_only` both empty.
Expectation: PASS.

P3. Exact ordered-slice reduction removes zero rows from each published system
supplied alone. Expectation: PASS. A removal would mean the published table is
redundant on the ordered slice, which would be a finding in its own right.

P4. `T_oracle(3,8) < 5` seconds on this machine. Expectation: PASS, because
only 31 rows reach the determinant stage.

P5. Enumeration is the dominant stage of `T_full(3,8)`, taking more than 50
percent of its wall clock. Expectation: PASS, following the stage attribution
already recorded in `docs/orbit_native_constructor.md`.

P6. `Ceiling(3,8) >= 5`. Expectation: PASS. This is the prediction that decides
whether the BDR route is worth further work at all.

## Pass and fail criteria for integration

A production integration of any BDR component is justified only if ALL of:

1. every emitted row is independently certified locally, with a replayed
   certificate and a bound oriented `tau`;
2. the retained systems at ranks 7 and 8 are reproduced exactly, by oriented
   `tau` set equality in both directions;
3. the completeness assumption actually used is stated, and is not upgraded
   from fixed `(N,d)` to uniform in `d`;
4. candidate count or stage-attributed runtime improves materially, measured,
   not argued;
5. the external code is not required at package runtime;
6. every failure mode remains fail closed.

Integration is REFUTED for this sprint if the measured ceiling shows that
removing local candidate generation entirely cannot pay for the integration,
or if the external code cannot be pinned and reproduced. The second condition
is already met by the availability record above, so the strongest outcome
available in this sprint is a bound plus a correspondence, and the honest label
for the benchmark itself is INCOMPLETE.

## Stop conditions

Stop integration work if external conventions cannot be bound exactly, if the
code cannot be pinned and reproduced, if a hybrid loses a known row, if the
only gain comes from weakening exactness, or if the external method has no
measured advantage after local recertification.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED per prediction, decided by the
artifact alone. Timing predictions report the machine. No prediction is edited
after the run. Scope fields: these predictions are about algorithm stages and
published systems, not about sampled fiber rows, so `independent_scope` is the
number of distinct published systems behind each verdict and is reported as
such.

## SCORED (2026-08-08, after running; predictions above are unedited)

Artifact: `results/data/bdr_constructor_comparison.json`. Narrative:
`docs/bdr_constructor_comparison.md`. Machine:
`Linux-6.18.5-fc-v20-x86_64`, CPython 3.14.0rc2, single worker.

| id | verdict | scope | independent scope | measured |
|---|---|---|---|---|
| P1 | PASS | 36 rows | 3 systems | 36 of 36 `certified`, every certificate replayed and bound |
| P2 | PASS | 2 systems | 2 | oriented `tau` sets equal both ways at `(3,7)` and `(3,8)` |
| P3 | PASS | 2 systems | 2 | 0 of 4 and 0 of 31 rows removed |
| P4 | PASS | 1 | 1 | oracle total 0.98 s at `(3,8)` |
| P5 | **FAIL** | 1 | 1 | enumeration share 0.391, not above 0.5 |
| P6 | PASS | 1 | 1 | oracle ceiling 16.38x at `(3,8)` |

`(3,6)` scores under P1 only. The composition layer guards ranks 7 through 13,
so rank 6 has no composition stage to reduce or time, and it is excluded from
P2 through P6 with the reason recorded in the artifact rather than dropped.

### P5, the failed prediction, recorded as a finding

The prediction was taken from `docs/orbit_native_constructor.md`, which states
that orbit enumeration "is now the dominant stage". At `(3,8)` on this machine
it is not: enumeration 6.25 s, screening 7.64 s, composition 2.11 s, so
screening dominates and enumeration is 39.1 percent. Three runs gave 0.380,
0.394 and 0.391, so the failure is not measurement noise.

That document's measurement was taken on a different machine, where the same
construction cost 66.0 s against 16.00 s here. Multi-machine results are
reported per machine and never averaged, so this is not a refutation of it. The
status is OPEN: the "enumeration dominates" claim does not reproduce here and
needs re-measuring on the original machine before it is quoted again.

### Added after scoring, and labelled as such

The birationality-prefilter path in the artifact was added AFTER the six
predictions were scored, and carries `added_after_scoring: true`. It is a
measured quantity, not a retro-fitted prediction, and no prediction above was
edited to accommodate it.

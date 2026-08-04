# Fixed-N=3 generalized Pauli constraint generator methodology

## Scope and proof discipline

This document specifies the fixed-particle-number generator for
`N=3`. It separates four logically different questions:

1. Did the discovery backend enumerate every relevant candidate?
2. Is a supplied candidate an exact Ressayre inequality?
3. Is a valid supplied inequality redundant on the ordered Pauli slice?
4. Does an independently computed Schubert coefficient agree with the
   retained result?

The implementation answers questions 2 through 4 with exact, replayable
certificates for the rows it receives. Question 1 remains an upstream proof
obligation outside the exhaustive `(3,6)` implementation. Consequently, a
new-rank output is a facet system relative to the supplied certified
candidates unless candidate exhaustiveness has been established separately.

Artifacts form a strict claim ladder:

1. a discovery candidate collection records what one backend emitted;
2. a screened collection classifies each supplied row under the implemented
   Ressayre criterion;
3. a relative H-system is irredundant only for the supplied certified rows;
4. a known-system match reproduces one sealed reference set; and
5. a complete GPC system additionally needs an exhaustive-enumeration proof
   or an independent exact inner/outer closure.

No level implies the next without its named certificate. The serialized field
`final_system` in a new-rank relative result must be read together with
`status=finalized_relative_to_supplied_candidates` and
`system_completeness=not_established_without_candidate_exhaustiveness`.

The following labels are used throughout:

- **Proved** means that the artifact contains exact data that the local
  verifier replays using integer or rational arithmetic.
- **Conditional** means that the conclusion follows if a named external
  enumeration or theorem-applicability obligation holds.
- **Operational** describes a search result, resource limit, or diagnostic
  state with no negative mathematical implication.
- **Proposed** describes rank-scaling architecture that is not yet a theorem
  or completed implementation.

No unresolved determinant search is called a rejection. No cyclic-Schubert
coefficient, including coefficient one, is used as a general facetness or
irredundancy criterion.

## Mathematical conventions

Let

\[
\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_d,
\qquad 0 \leq \lambda_i \leq 1,
\qquad \sum_i \lambda_i = 3.
\]

The representation is \(\bigwedge^3 \mathbb C^d\). Its weights are

\[
\omega_I = \sum_{i\in I} e_i,
\qquad I\subseteq\{1,\ldots,d\},\quad |I|=3.
\]

The Bulois-Denis-Ressayre backend writes an oriented homogeneous inequality
as

\[
\tau\mathbin{\cdot}\lambda \leq 0,
\qquad \tau\in\mathbb Z^d.
\]

Positive common factors are removed without changing orientation. An affine
row \(a\mathbin{\cdot}\lambda\leq b\) and its homogeneous row are related on
the trace-three slice by

\[
\tau_i = 3a_i-b.
\]

Conversely, `constraint_from_tau` chooses a deterministic nonnegative affine
representative. Round trips are compared after primitive positive
normalization, never after a sign change.

For Ressayre certification, put \(S=\sum_i\tau_i\) and define

\[
H_i=S-d\tau_i,
\qquad z=3S.
\]

After division by the positive common content, the inequality
\(H\mathbin{\cdot}\lambda\geq z\) is exactly
\(\tau\mathbin{\cdot}\lambda\leq0\) on the trace slice. The code uses the
negative type-A roots \(e_j-e_i\), with \(i<j\), for which
\(H_j-H_i<0\).

## Implemented layers

| Layer | Implementation | Exact conclusion | Deliberate nonclaim |
|---|---|---|---|
| Exhaustive small-rank enumeration | `generation.ressayre.enumerate_admissible_hyperplanes` | Every affine exterior-weight hyperplane is enumerated by the stated brute-force argument | Scalable beyond the small-rank gate |
| Exact closed-flat enumeration | `generation.admissible_flats` | Every affine exterior-weight hyperplane is enumerated through rank-preserving projective residual closure; proved complete at rank 7 | Practical scaling to rank 13 |
| External candidate discovery | `scripts/bdr_fixed_n3_backend.py` | A revision-checked candidate export and provenance sidecar are produced; downstream artifacts bind the export hash | Local proof that the external prefilter is exhaustive |
| Parse and normalization | `generation.bdr` | Safe literal parsing, primitive orientation-preserving normalization, deduplication, and structural-row separation | Candidate validity or completeness |
| Candidate screening | `generation.candidate_pipeline` | Every supplied row receives an exact positive, exact negative, or explicit unresolved status | Facetness, redundancy, or exhaustiveness |
| System composition | `generation.candidate_system` | Only determinant-certified rows enter exact reduction; unresolved rows block finalization | Completeness beyond the supplied candidates |
| Redundancy and facet proof | `generation.redundancy` | Exact implications and relative facet witnesses on the ordered trace-three Pauli slice | Representation-theoretic validity of its input rows |
| Retained-row cross-check | `generation.cyclic_schubert` | Exact characteristic-zero coefficient, with optional modular replay | A general coefficient-one validity or irredundancy gate |
| Published-system regression | `generation.candidate_system` | Exact oriented homogeneous set equality for ranks with a lookup oracle | Proof that the discovery algorithm is complete at every rank |

## Candidate enumeration

### Exhaustive admissible-hyperplane reference algorithm

The local reference enumerator uses the fact that the exterior weights lie in
the trace affine space of dimension \(d-1\). An admissible codimension-one
hyperplane contains \(d-1\) affinely independent exterior weights. The
algorithm enumerates every \((d-1)\)-subset of weights, computes the centered
integer normal by exact cofactors, primitive-normalizes the pair \((H,z)\),
and deduplicates it.

This is a completeness proof for the reference enumeration because any
admissible hyperplane contains at least one such affine basis. It is not a
scaling strategy. With

\[
W=\binom d3,
\]

the raw subset count is \(\binom{W}{d-1}\).

### Exact closed-flat reference algorithm

The rank-7 reference replaces repeated affine bases by the lattice of closed
flats of the augmented exterior weights. From each closed rank-\(k\) flat it
groups exterior weights outside the flat by projective residual class. Each
class is exactly one rank-\(k+1\) closure. Induction generates every flat.

Rank calculations use a prime strictly above the Hadamard bound on every
augmented-weight minor, so reduction modulo that prime preserves all rational
ranks and closures. At rank 7 this enumerates 5,341 hyperplanes in place of
1,623,160 candidate bases. Both orientations give 10,682 rows. Exact
screening classifies all of them with no unresolved determinant: 10,133 trace
failures, 499 exact determinant zeros, one structural Pauli row, and 49
nonstructural Ressayre rows. Exact reduction returns the four published
facets. This closes candidate exhaustiveness and system completeness at rank
7, while remaining a coverage oracle rather than the proposed rank-13 search.

### Pinned BDR discovery backend

Ranks above the brute-force regime use the optional upstream backend
`ea-icj/moment_cone`, pinned to revision

`7ab347d9a7837d68d92bde9fac74606913f395ca`.

The wrapper performs the following reproducibility checks and modifications:

1. It reads the installed PEP 610 revision, or the revision of a local
   checkout, and rejects a mismatch by default.
2. It rejects tracked changes in a local checkout by default.
3. It replaces the randomized stabilizer filter by a pass-through.
4. With `--candidate-superset`, it requests an empty upstream inequality
   filter list.
5. With `--candidate-output`, it captures pending and validated homogeneous
   rows immediately after the upstream candidate step.
6. It records expected and observed revisions, dirty status, delegated
   arguments, filter mode, output path, and run completion in a JSON sidecar.

The override flags `--allow-unpinned-backend` and
`--allow-dirty-backend` make a run explicitly provisional. Such a run cannot
support an exhaustiveness claim.

The current sidecar is not yet a complete proof manifest. In particular, the
wrapper checks tracked changes but deliberately ignores untracked files, does
not hash the installed source tree or wrapper file, and does not bind the raw
export hash in the discovery sidecar. The downstream screening commands do
hash the exact export bytes. A proof-eligible discovery artifact must close
all three provenance gaps.

**Current proof status: Conditional.** The wrapper removes two known
probabilistic loss points and preserves a prefilter candidate export, but the
local repository does not yet contain an independent proof that every branch
of the pinned upstream `TauCandidates` and submodule-condition computation is
complete for this specialization. Until that audit is closed, the export is
a supplied candidate collection, not a proved exhaustive list. The word
`superset` in the wrapper flag names the requested prefilter mode; it is not a
local theorem that the resulting file contains every required candidate.

An audit of the pinned path gives the more precise call graph. `TauCandidatesStep`
calls `tau.find_1PS`; `SubModuleConditionStep` then filters those rows by
`Tau.is_sub_module`; the wrapper bypasses `StabilizerConditionStep`; and
`InequalityCandidatesStep` calls `List_Inv_Ws_Mod` before the wrapper captures
pending and validated inequalities. With an empty inequality-filter list,
capture precedes the explicitly probabilistic `PiDominancy` and birationality
filters. `List_Inv_Ws_Mod` itself is a deterministic partition recursion over
compatible inversion sets.

Three coverage gaps remain. The pinned source does not supply a proof that
the recursive pruning in `find_1PS` emits every required one-parameter
subgroup, does not justify in this repository that `is_sub_module` is a
completeness-preserving rejection, and does not prove that the symmetry and
partition bounds in `List_Inv_Ws_Mod` cover every required Weyl representative.
In addition, the separate `GeneralStabilizerDimensionCheck` still calls the
randomized stabilizer routine unless delegated arguments include
`--no_dim_check`; bypassing `StabilizerConditionStep` does not bypass that
earlier check. Therefore the exact proof status is: the wrapper faithfully
captures the pinned prefilter output under recorded arguments, but
candidate-exhaustiveness is unproved and every derived new-rank H-system is
relative to that captured collection.

## Exact candidate screening

### Normalization and structural rows

`screen_tau_candidates` first validates integral row widths, removes positive
common factors, sorts and deduplicates the primitive rows, and separates

\[
(0,\ldots,0,-1)\mathbin{\cdot}\lambda\leq0,
\]

which is the structural inequality \(\lambda_d\geq0\). Count-conservation
checks bind raw rows, normalized rows, structural rows, and nonstructural
candidates.

For the ranks finalized by `candidate_system`, a deterministic finite-field
zero-stabilizer certificate proves the ambient moment cone is full
dimensional. That certificate is replayed independently before reduction.

### Admissibility

For a primitive oriented \(\tau\), the code constructs \((H,z)\), selects the
weights satisfying \(H\mathbin{\cdot}\omega=z\), and proves that these
weights have affine rank \(d-1\). A modular nonzero minor is already an exact
rank witness over \(\mathbb Q\); symbolic rational rank is the fallback. A
row failing this condition is `inadmissible` for the Ressayre criterion. This
is not, by itself, a witness that the associated halfspace is violated by the
moment polytope.

### Trace condition

Let

\[
\Phi_{<} = \{\omega:H\mathbin{\cdot}\omega<z\}
\]

and let \(R_{<}\) be the selected negative roots. The exact Ressayre trace
condition is

\[
|\Phi_{<}|=|R_{<}|.
\]

A mismatch is `trace_rejected`. It is an exact negative result for Ressayre
certification, not a general proof that the written inequality is false. A
redundant valid inequality need not itself be a Ressayre facet element.

### Determinant nonvanishing

The tangent matrix has rows indexed by \(\Phi_<\), columns indexed by
\(R_<\), and entries that are integral linear forms in variables attached to
weights on the hyperplane. Exterior-root actions determine every matrix
entry, including its fermionic sign.

At a deterministic integer assignment, the code constructs the integer
matrix and computes its determinant by fraction-free Bareiss elimination. A
nonzero value is a mathematical certificate that the determinant polynomial
is nonzero. The serialized certificate records the oriented \((H,z)\), the
admissibility witness, on and below weight indices, negative roots, integer
evaluation, and determinant value. Verification rebuilds the matrix, replays
the determinant, and binds \((H,z)\) back to the exact oriented \(\tau\).

Trying finitely many assignments and finding only zero does not prove that
the polynomial is zero. That outcome is `evaluation_unresolved`.

### Symbolic fallback

An optional fallback constructs the determinant polynomial exactly. It is
attempted only after admissibility and trace equality have been replayed. A
matrix-order limit can gate it, but the automatic rank-at-most-13 default has
no order cap unless the caller supplies one.

- A nonzero polynomial is converted to an exact nonzero integer evaluation,
  replayed, and reported as `certified`.
- Exact polynomial zero is `symbolic_zero_rejected` for the Ressayre
  nonvanishing condition.
- A configured determinant-order limit is
  `resource_limit_unresolved` inside the row's fallback payload.
- An operational exception is `exception_unresolved` inside the fallback
  payload.

The last two cases leave the row's final status
`evaluation_unresolved`. They are not rejection evidence. The screening API
can process independent rows in worker processes, but ordered mapping keeps
the serialized result deterministic.

### Complete status partition

| Final row status | Mathematical meaning | May enter reduction? |
|---|---|---|
| `certified` | Replayed nonzero Ressayre determinant bound to oriented `tau` | Yes |
| `trace_rejected` | Exact failure of the Ressayre trace condition | No |
| `symbolic_zero_rejected` | Exact failure of Ressayre determinant nonvanishing after trace equality | No |
| `inadmissible` | Exact failure of Ressayre exterior-weight admissibility | No |
| `evaluation_unresolved` | No positive witness and no exact negative proof | No, and it blocks finalization |

The `--require-resolved` screen-only command writes its artifact first and
then exits with status 2 when unresolved rows remain.

## Exact reduction and facet certificates

`compose_candidate_system` replays the screening object before using it. It
requires ranks supported by the full-dimension and ordered-slice certificate,
requires candidate-wide cyclic checks to have been deferred, verifies every
positive determinant certificate, and confirms that every negative or
unresolved row carries no positive certificate.

If any row remains unresolved, composition returns
`blocked_unresolved_candidates`. No reduction is run and no final system is
exposed.

Otherwise, only certified rows are passed to the ordered slice

\[
\sum_i\lambda_i=3,
\quad \lambda_1\leq1,
\quad \lambda_1\geq\cdots\geq\lambda_d,
\quad \lambda_d\geq0.
\]

For a removed row \(a\mathbin{\cdot}x\leq b\), the certificate gives
nonnegative rational multipliers \(y_j\), a free trace multiplier \(z\), and
nonnegative slack such that

\[
a=\sum_j y_j A_j+z\mathbf1,
\qquad
b=\sum_j y_j b_j+3z+\operatorname{slack}.
\]

This is an exact Farkas implication. Removals are replayed in order, so the
proof has no circular dependency.

For every retained row, a rational point satisfies all other final rows and
strictly violates the target. A shared rational point satisfies every final
row strictly. The segment between these points intersects the target
hyperplane while remaining strict for all other rows. Therefore the target is
a relative facet of the final polyhedron in the trace slice and is
irredundant there.

Floating-point linear programming discovers supports and approximate points
only. Every returned multiplier and point is reconstructed over
\(\mathbb Q\) and replayed exactly. Reconstruction failure raises
`CertificateSearchError`; it is not converted into a mathematical verdict.

## Cyclic-Schubert cross-validation

Exact cyclic-Schubert evaluation is deferred until after redundancy
elimination. This changes its workload from every certified candidate to only
the retained rows. The characteristic-zero divided-difference result is
stored in full. An optional prime-field computation independently checks its
residue. This is one selected cyclic coefficient, not a complete general
Schubert-calculus engine.

The cross-check never changes the Ressayre status. Coefficient zero, a
positive nonunit coefficient, coefficient one, a degree mismatch, and an
unavailable test are preserved as exact diagnostic outcomes. In particular:

- coefficient one is not assumed to be a general validity theorem;
- a positive nonunit coefficient is not rejected;
- coefficient zero concerns only the selected cyclic test;
- degree mismatch or a non-weight threshold makes that test unavailable;
- a nonzero modular residue proves the exact integer is nonzero, while a zero
  residue alone is inconclusive;
- facetness comes from the ordered-slice certificate, not from the
  coefficient value;
- candidate completeness does not follow from a retained-row coefficient.

A genuinely independent retained-row audit should use a second
implementation or a general Schubert route rather than treating exact and
modular evaluation of the same recurrence as complete independence.

`generation.series` and `scripts/fixed_n3_series.py` predate the candidate
screen and deliberately require coefficient one for their known final-system
regression payloads. They are not the general candidate-validity API. New
candidate or frontier work must use `candidate_pipeline` followed by
`candidate_system`, which preserves nonunit coefficients and keeps cyclic
evaluation separate from Ressayre certification.

## Regression and conditional completeness

For published ranks, the composition layer compares primitive oriented
homogeneous rows with the lookup system by exact set equality. A mismatch
returns `blocked_known_system_mismatch` and exposes no final system. Passing
this gate is an exact software and data regression for that rank. It is not a
general proof of the upstream enumeration algorithm.

For a rank without a lookup oracle, successful composition is labeled
`finalized_relative_to_supplied_candidates`.

The complete-system implication is the following conditional proposition.

### Conditional completeness proposition

Assume:

1. for every nonstructural Ressayre element required by the cited
   full-dimensional completeness theorem, the export contains a positive
   multiple of its oriented primitive `tau`, modulo only explicitly proved
   symmetry reductions;
2. the ambient full-dimension hypothesis and the applicability of the stated
   Ressayre characterization are proved;
3. every supplied row is resolved exactly as certified or rejected;
4. every positive certificate and every reduction certificate passes its
   verifier; and
5. the structural Pauli and dominance rows describe the intended ordered
   slice.

Then the retained rows and structural rows give a complete irredundant facet
description of that ordered moment polytope.

**Proof.** Assumption 1 reduces validity and completeness to the exhaustive
candidate list. Assumptions 2 through 4 classify each candidate without
discarding an unresolved case. The exact Farkas certificates show that every
removed certified row is implied by the rows active at its removal. The
facet witnesses prove every retained row is necessary and facet defining on
the slice in assumption 5. Therefore reduction preserves the intersection
defined by all valid candidates and returns its irredundant facets.

The proposition is proved as a composition statement. Its first assumption
is currently conditional above the exhaustive small-rank reference gate.

### Independent inner/outer closure route

Candidate-enumeration completeness can also be bypassed for a proposed
H-description when all of the following are exact:

1. every H-row has a validity proof;
2. the H-polyhedron is vertex-enumerated without omission; and
3. every enumerated H-vertex has an exact attainable-state certificate.

Validity gives \(P\subseteq H\), where \(P\) is the moment polytope. Exact
attainability of every H-vertex gives \(H\subseteq P\) by convexity. Hence
\(P=H\). This is the inner/outer sandwich used by the exhaustive `(3,6)`
reference gate. Numerical sampling or an incomplete vertex list cannot replace
either exact inclusion.

### Blind known-system regression protocol

The strongest reproduction test keeps lookup data inaccessible during
discovery, screening, and reduction, seals and hashes the generated artifact,
and only then loads the reference system. It compares primitive oriented
`tau` sets in both directions and preserves exact differences. The current
local composition performs the exact bidirectional set comparison; campaign
automation should additionally enforce the lookup-isolation and precomparison
hash steps. A pass proves reproduction of the sealed reference, not the
general completeness of the discovery algorithm.

## Complexity

Let

- \(r\) be the number of supplied rows;
- \(c\) be the number of determinant-certified rows;
- \(f\) be the number retained after reduction;
- \(W=\binom d3\) be the number of exterior weights;
- \(R=d(d-1)/2\) be the number of negative type-A roots before filtering;
- \(q=|\Phi_<|=|R_<|\) for a trace-matched candidate;
- \(m\) be the number of weights on its hyperplane; and
- \(a\) be the deterministic evaluation-attempt limit.

Primitive normalization and sorting cost
\(O(rd\log r)\) integer comparisons plus gcd arithmetic. Weight-level scans
cost \(O(Wd)\) per candidate. Modular affine-rank discovery is polynomial in
the number of on-hyperplane rows and dimension, with exact rational fallback.

One evaluated tangent matrix costs \(O(qmd)\) operations in the current dense
weight representation because each exterior-root action scans and rebuilds a
length-\(d\) weight. Bareiss elimination costs
\(O(q^3)\) integer operations, excluding coefficient bit growth. Thus the
simple upper bound for evaluated screening is

\[
O\bigl(r(Wd+R)+a\sum(qmd+q^3)\bigr).
\]

Exact symbolic determinant expansion can have exponential expression growth.
It is optional, resource-gated, and never required to turn finite evaluation
failure into rejection.

The ordered-slice reducer solves \(O(c)\) discovery linear programs with
\(O(c+d)\) constraints or dual variables. Its exact replay costs
\(O(c^2d)\) rational operations in the stated worst-case bound, and the
serialized certificate is \(O(cd)\) for basic sparse Farkas supports.

Cyclic divided differences can have exponential state growth in their word
length. Deferring them changes the number of calls from \(c\) to \(f\).
Candidate screening is embarrassingly parallel across rows, although memory
use and symbolic expression growth limit useful worker counts.

## Architecture toward rank 40

At \(d=40\),

\[
W=\binom{40}{3}=9880,
\qquad R=\binom{40}{2}=780.
\]

The brute-force count \(\binom{9880}{39}\) rules out subset enumeration.
Dense determinant and global redundancy steps also become dominant. The
following architecture is proposed; only the exact local certificate kernels
and row-parallel screening are implemented today.

1. **Proof-carrying candidate pruning.** Generate candidates by Weyl-chamber,
   Bruhat, trace, and submodule conditions, but associate every discarded
   class with a proved pruning lemma. Empirical small-edge or cyclic-window
   patterns may prioritize work but may not delete candidates.
2. **Canonical orbit representatives.** Quotient only by symmetries proved to
   preserve the dominant problem. Store the orbit action and orientation so
   every representative expands back to the complete candidate set.
3. **Streaming normalization.** Primitive-normalize and hash rows as they are
   generated. Shard by canonical signatures and preserve deterministic merge
   order, source revision, command line, and content hashes.
4. **Scalable ambient certificates.** Replace dense action-matrix storage by
   sparse modular zero-stabilizer witnesses and verify them rank by rank. The
   current composition and ordered-slice APIs intentionally stop at rank 13;
   remove those guards only after the full-dimension and verifier contracts
   are proved for the larger range.
5. **Cheapest exact filters first.** Apply admissibility and trace counts
   before determinant work. Cache exterior weights, level signatures, root
   sets, affine-rank witnesses, and repeated tangent sparsity patterns.
6. **Sparse nonvanishing certificates.** Build tangent matrices sparsely. A
   nonzero determinant residue modulo a recorded prime at an exact integer
   evaluation is already a characteristic-zero nonvanishing proof. Replay
   selected positive rows over integers or independent primes as an audit.
7. **Fail-closed zero handling.** Use symbolic factorization or structured
   matching only where it proves zero or constructs a witness. Time, memory,
   or expression limits remain unresolved and block finalization.
8. **Incremental exact reduction.** Use sparse active-set discovery, parallel
   witness searches, and checkpointed rational Farkas reconstruction. The
   final merged certificate must still replay sequentially without a cycle.
9. **Cross-check only retained orbits.** Run exact cyclic-Schubert evaluation
   after reduction, with modular replay and a second implementation on audit
   samples.
10. **Stabilization as a theorem target.** Record primitive normal signatures,
   padding behavior, coefficient patterns, and orbit data across ranks.
   Promote a stable family to a pruning or generation rule only after proving
   both forward validity and converse coverage.
11. **Hierarchical verification.** Verify each row and shard independently,
    then verify merge conservation, reduction implications, facet witnesses,
    and the final artifact hash. A failed shard leaves the global rank
    unresolved rather than partially finalized.
12. **End-to-end provenance.** Hash the candidate export in the discovery
    sidecar as well as downstream artifacts. Record the backend revision,
    package metadata, dirty state, delegated arguments, Python and dependency
    versions, deterministic seeds, worker count, shard hashes, and merge
    manifest. Overrides must mark every descendant artifact provisional.

The rank-40 objective therefore depends less on faster symbolic algebra than
on proof-preserving candidate compression, sparse certificates, and resumable
artifacts. Any shortcut that cannot explain which candidates it removed is a
heuristic, not part of the definitive generator.

## Reproducible command layers

The discovery wrapper is optional and external:

```text
python scripts/bdr_fixed_n3_backend.py \
  --candidate-superset \
  --candidate-output CANDIDATES.py \
  --provenance-output CANDIDATES.provenance.json \
  [pinned upstream arguments]
```

Screening alone preserves every exact and unresolved outcome:

```text
PYTHONPATH=src python scripts/screen_fixed_n3_candidates.py \
  --input CANDIDATES.py --rank D --output SCREEN.json \
  --workers P --require-resolved
```

Screening, exact reduction, retained-row cross-validation, and published-rank
regression are composed by:

```text
PYTHONPATH=src python scripts/generate_fixed_n3_candidate_system.py \
  --input CANDIDATES.py --rank D --output SYSTEM.json --workers P
```

The two local screening commands retain source hashes and proof-status fields.
The discovery wrapper records revision and command provenance, while direct
export hashing at that layer remains a rank-scaling requirement above. A
nonzero exit caused by unresolved candidates or a known-system mismatch is a
successful fail-closed diagnostic, not permission to discard the blocking
rows.

## Verification map

- `tests/test_bdr_bridge.py` checks literal parsing, trace-gauge conversion,
  normalization, and lookup round trips.
- `tests/test_bdr_fixed_n3_backend_script.py` checks revision, dirty-tree,
  candidate capture, and provenance behavior of the wrapper.
- `tests/test_ressayre_tau_bridge.py` checks admissibility, exact evaluation,
  orientation binding, unresolved semantics, and tamper rejection.
- `tests/test_candidate_pipeline.py` checks the full status partition,
  symbolic zero, resource failure, cyclic deferral, count conservation, and
  full-dimension separation.
- `tests/test_redundancy.py` checks exact Farkas replay, rational facet
  witnesses, structural rows, and tamper rejection.
- `tests/test_candidate_system.py` checks fail-closed composition, known-rank
  regression gates, reduction ordering, and retained cyclic cross-checks.
- `tests/test_cyclic_schubert.py` checks exact and modular divided-difference
  certificates without treating coefficient one as a general gate.

The exhaustive `(3,6)` reference proof and its independent inner/outer
closure are documented in `docs/stage1_ressayre_3_6.md`. The complete
closed-flat `(3,7)` proof is documented in
`docs/stage1_ressayre_3_7.md`. The broader
representation-theoretic specification and its still-unproved pruning ideas
are recorded in `docs/stage1_klyachko_spec.md`.

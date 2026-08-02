# Execution roadmap

## Phase 0: Consolidate the branch

- Apply all research patches to a single branch based on current `main`.
- Remove duplicate scripts and merge overlapping theorem notes.
- Regenerate provenance and checksums.
- Run the complete repository test suite.
- Produce one auditable diff.

**Exit criterion:** clean branch, all tests passing, no result described only in
chat.

## Phase 1: Theory freeze

- Formalize T1--T8.
- Add assumptions and nonclaims to every theorem.
- Produce second-path sign and rank verifiers.
- Search small carriers exhaustively for counterexamples.
- Determinacy closure on the unresolved 28 is done and resolved none
  (`results/data/rdm_determinacy_closure.json`). What remains is deciding
  whether to attempt the nonlinear rank-one closure the artifact's nonclaim
  field leaves open, or to report 771/799 as final for the order-two
  certificate.
- Add a claims ledger mapping each sentence to proof, exact artifact, numerical
  evidence, or conjecture.

**Exit criterion:** referee can reproduce every theorem and exact count.

## Phase 2: Prior-art and adversarial review

- Build the closest-result table from `docs/prior_art_roadmap.md`, which
  already carries the outside-GPC map with primary citations, rather than
  starting a second survey.
- Write an internal hostile referee report.
- Commission or solicit independent mathematical review.
- Downgrade or remove any claim with ambiguous novelty.
- Verify bibliography and quotations.

**Exit criterion:** novelty statement survives comparison with geometric
quantum mechanics, marginal uniqueness, phase retrieval, TD-2RDM, matrix
completion, active-space, and tensor-network literature.

## Phase 3: Paper

- Draft manuscript in `paper/observable_manifolds/`.
- Separate main theorem proofs from computational supplement.
- Include one finite transport model and one infinite family.
- Put the 799-record atlas in the supplement, not at the narrative center.
- Create reproducibility bundle and standalone verifier.
- Run journal-specific style and submission checks.

**Exit criterion:** complete manuscript, supplement, cover letter, referee
checklist, and reproducibility archive.

## Phase 4: Reference package

- Create `obdyn` package.
- Stabilize schemas and APIs.
- Implement automatic carrier certification and atlas construction.
- Implement minimal regular charts plus Cartesian boundary fallback.
- Implement compiled sparse RHS.
- Add direct-wavefunction verification backend.

**Exit criterion:** installable package reproduces all paper figures and tables.

## Phase 5: Integrations

- PySCF adapter.
- OpenFermion adapter.
- FCIDUMP adapter.
- Import measured 2-RDM data with uncertainty.
- Export reconstructed higher RDMs and diagnostics.

**Exit criterion:** an external electronic-structure calculation can enter and
leave the observable workflow without custom glue code.

## Phase 6: Performance and visualization

- Numba/JAX compiled kernels.
- Symplectic midpoint/variational integrator.
- Atlas optimizer based on conditioning and Hamiltonian sparsity.
- Interactive carrier/coherence graph viewer.
- Benchmark suite across FCI, carrier CI, dense TD-2RDM, sparse observable,
  and tensor-network baselines where appropriate.

**Exit criterion:** benchmark claims are reproducible on at least two machines
and include equally optimized controls.

Partial status. The four existing benchmarks have been run on a second
machine: variable counts, evaluation and step counts, and accuracies
reproduced, while every runtime factor moved (compiled kernel 1.44 to 1.17,
exact propagator versus carrier 2.29 to 1.76). Method ordering held
throughout. The equally-optimized-control half is not met yet: the dense
2-RDM baseline wraps the sparse implementation rather than being an
independent one, and the generalized family's honest control is the
factorized amplitude update, not the dense 2^q carrier.

## Phase 7: Release strategy

- Preprint theory after novelty and proof audit.
- Tag the exact reproduction release.
- Keep performance kernels and integration roadmap moving beyond the paper.
- Decide licensing only after identifying external users and a defensible
  support/integration value proposition.

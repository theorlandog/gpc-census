# Agent Track: Entanglement Polytopes, Quantum Marginals, and Spin Systems

## Mission
Test whether the current moment-polytope / fiber / observable-range machinery generalizes into a useful state-certification engine outside fermions.

These topics belong together because entanglement polytopes and spin/multipartite moment polytopes are instances of the broader quantum marginal / moment-map problem.

## Phase 1: entanglement-polytopes reproduction
Select low-dimensional canonical multipartite systems where entanglement polytopes are known.
- reproduce published local-spectra polytopes from first principles or trusted reference data;
- identify active faces and associated state restrictions;
- verify whether local spectra can exclude entanglement classes with exact certificates;
- quantify how many additional measurements are needed to separate residual classes.

## Phase 2: minimal certification measurements
For a target entanglement class or property:
- formulate compatible-state fibers from local spectra / reduced states;
- find minimal extra observables needed to certify membership or nonmembership;
- compare with standard entanglement witnesses, tomography, randomized measurements, and shadow-based methods;
- prefer problems where local spectra alone already collapse the candidate set sharply.

## Phase 3: general quantum marginal backend
Abstract the pipeline:
`representation + subgroup action + moment data -> polytope face -> compatible variety/fiber -> target range`.
Determine exactly which components of gpc-census are fermion-specific and which can be generalized.
Build the smallest reusable interface, not a giant framework.

## Phase 4: spin-system pilot
Use a small spin / qubit many-body representation:
- fixed magnetization / total spin / symmetry sectors;
- local reduced information;
- target correlators or Hamiltonian energy;
- test whether moment/polytope constraints reduce needed measurements in a way conventional grouping would not know.

Favor trapped-ion / neutral-atom / spin-simulation style observables if a clean benchmark exists.

## Phase 5: commercial verdict
Score use cases:
- entanglement certification without tomography;
- many-body simulator validation;
- quantum marginal consistency checking;
- hardware experiment preflight.

Require a concrete expensive workflow and a credible buyer before recommending productization.

## Deliverables
- `docs/qvalidate_entanglement_marginal_spin.md`
- `results/data/qvalidate_entanglement_marginal_spin.json`
- exact reproduction tests for known polytopes
- at least one minimal-measurement benchmark
- `HANDOFF_ENTANGLEMENT_MARGINAL_SPIN.md`

## Guardrails
- Browse literature aggressively before novelty claims.
- Do not confuse class exclusion from local spectra with full state certification.
- Do not claim scalability from tiny multipartite examples without evidence.

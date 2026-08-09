# Agent Track: QEC, Stabilizer, and Magic-State Validation

## Mission
Test whether the qvalidate idea is substantially more valuable in quantum error correction and fault-tolerant resource-state validation than in fermionic chemistry.

These topics are related and should be handled by one agent serially because stabilizer geometry, logical-observable identifiability, and magic-resource certification share algebraic structure and baselines.

## Phase 1: stabilizer/QEC feasibility
Choose small canonical codes first: repetition code, 5-qubit code, Steane [[7,1,3]], surface-code patches if tractable.
For each:
- represent the code space and stabilizer constraints exactly;
- define target logical observables and figures of merit: `<X_L>`, `<Z_L>`, logical fidelity, syndrome-conditioned quantities;
- determine the minimum physical observables needed to identify/certify each target under the code-space assumption;
- compare against physical-state tomography and obvious stabilizer-measurement baselines;
- distinguish state certification from syndrome extraction already required operationally.

## Phase 2: noise / near-code-space robustness
Extend from exact code-space membership to approximate membership:
- bounded stabilizer violation;
- leakage outside code space;
- noisy measurement intervals;
- derive rigorous target-error bounds where possible.

The economically interesting question is not "can stabilizers characterize a code?" but:
`Given measurements already available or cheap to obtain, what additional information is minimally sufficient to certify a logical claim?`

## Phase 3: magic-state / non-Clifford resource validation
Study small stabilizer polytopes / resource theories of magic.
Test problems such as:
- certify non-stabilizerness from sparse observables;
- certify a lower bound on a chosen magic monotone or threshold-relevant quantity;
- determine minimal measurements needed to distinguish a target magic state from the stabilizer polytope under noise;
- compare against known witnesses / tomography / randomized protocols.

Do not invent a novel polytope if the standard convex body already answers the question.

## Phase 4: fault-tolerant commercial relevance
Estimate whether any demonstrated reduction affects a genuinely expensive workflow:
- logical-qubit characterization;
- code-space verification;
- magic-state factory acceptance testing;
- logical observable validation.

Require at least one benchmark with either >=10x reduction in measured settings/shots or an exact certification capability unavailable to the obvious baseline.

## Deliverables
- `docs/qvalidate_qec_magic.md`
- `results/data/qvalidate_qec_magic.json`
- executable benchmark(s)
- standalone verifier(s)
- comparison table of codes / targets / baseline / optimized measurements / robustness
- `HANDOFF_QEC_MAGIC.md` with a go/no-go verdict

## Guardrails
- Search current literature before novelty claims.
- Do not count stabilizer measurements the hardware must already perform as "savings" unless the comparison is explicit.
- Keep fault-tolerance resource estimation separate from measurement identifiability.
- A cute small-code identity is not a business result.

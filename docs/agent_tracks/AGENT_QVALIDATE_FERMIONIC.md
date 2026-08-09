# Agent Track: Fermionic qvalidate / GPC Commercial Validation

## Mission
Determine whether GPC-aware state-space constraints can produce acquisition-grade reductions in real quantum measurement cost, and if so, identify the narrowest defensible product wedge.

## Work serially in this order

### Phase 1: rebase and inventory
- Pull current `main` from https://github.com/theorlandog/gpc-census.
- Record HEAD SHA.
- Inventory the current observability atlas, the 28 unresolved carriers, fixed-1RDM fiber machinery, rank-11 frontier artifacts, and any new qvalidate-related work.
- Do not trust prose summaries when an artifact can be regenerated.

### Phase 2: hard-carrier identifiability
Focus on the 28 carriers where the collision-free 2-RDM graph certificate is silent.
For each carrier:
- construct the exact polynomial map from carrier amplitudes/coherences to selected 2-RDM entries;
- quotient global phase and any obvious gauge;
- test generic injectivity, finite-to-one behavior, degree, and birationality where tractable;
- distinguish exact proofs from Monte Carlo evidence;
- search for minimal subsets of 2-RDM entries that determine the state or a specified target observable;
- emit exact counterexamples when injectivity fails.

Prefer symbolic elimination, Jacobian rank, Groebner/resultant methods, exact finite-field/rational witnesses, and independent replay over floating-point-only conclusions.

### Phase 3: target-observable optimization
Do not optimize full-state reconstruction unless it is genuinely cheaper.
For each promising carrier and target family:
- compare direct target estimation with reconstruction-based estimation;
- optimize the actual number of real measured observables / Pauli groups / shots under a stated model;
- test one-target and many-target amortization separately;
- include higher-RDM targets and dynamics where exact lower-RDM closure is available;
- use a strict kill criterion: no acquisition-grade claim unless the method gives >=10x reduction in a realistic cost metric or enables an exact capability a baseline cannot certify.

### Phase 4: chemistry mapping
Map the best structural examples to real or realistic electronic Hamiltonians:
- H2, LiH, BeH2 or other small molecules first;
- then at least one larger toy active space where GPC structure matters;
- compare against conventional Pauli grouping, shot allocation, classical shadows / fermionic shadows, and low-rank RDM completion where relevant;
- quantify chemical-accuracy target and cost reduction.

### Phase 5: quasipinning
Investigate whether approximate facet saturation yields rigorous state-space or observable-error bounds.
Target theorem form:
`D(lambda) <= eps => ||(I-P_S) psi||^2 <= C eps`
or a directly useful observable bound.
Do not claim a theorem without proof.

### Phase 6: commercial verdict
Produce a short decision memo:
- strongest surviving use case;
- exact technical moat;
- realistic buyers;
- likely pilot price;
- whether this justifies a proprietary Rust engine;
- whether to continue, narrow, or kill the company thesis.

## Deliverables
- `docs/qvalidate_fermionic_campaign.md`
- `results/data/qvalidate_fermionic_campaign.json`
- generator script(s)
- standalone verifier(s)
- tests that fail closed when artifacts are absent/stale
- a ranked table of best commercial benchmarks
- a concise `HANDOFF_QVALIDATE_FERMIONIC.md`

## Guardrails
- Never use Hilbert-space compression as a proxy for measurement savings.
- Never compare only to a straw-man baseline.
- Separate exact proof, rigorous bound, exhaustive finite search, and empirical evidence.
- Keep commercially general optimization ideas out of GPL code if they are not required for the research result; document proposed proprietary boundary instead.

# Coordination and Parallelization Plan

Run these tracks in parallel:

1. `AGENT_GENERATOR_FRONTIER.md`
2. `AGENT_QVALIDATE_FERMIONIC.md`
3. `AGENT_QEC_MAGIC.md`
4. `AGENT_ENTANGLEMENT_MARGINAL_SPIN.md`
5. `AGENT_ANYON_NEWTON_SCOUT.md` (lower priority)

## Base state

These tracks were refreshed against HEAD `476bf0b` ("Prove Conjecture S: the
pattern set determines the values", #127). One correction was applied at that
rebase and is recorded here so no agent reintroduces it: the level-spread
question named "Conjecture S" in the original generator track is PROVED at that
HEAD and is now Theorem S in `docs/profile_native_generation.md`. The generator
track's Phase 2 has been rewritten against the open problem that replaced it.
Agents should re-derive the current state rather than trusting this paragraph
if their HEAD differs.

## Why these are decoupled
- Generator work changes candidate/polytope production and proof completeness.
- Fermionic qvalidate consumes existing certified carriers and can proceed without new generator results.
- QEC/magic uses stabilizer/resource-state structure rather than fermionic moment polytopes.
- Entanglement/marginal/spin develops a generalized moment-map backend.
- Anyon/Newton is exploratory and should not block anything.

## Serial dependencies inside tracks
- QEC before magic: stabilizer geometry and logical-observable machinery are reusable.
- Entanglement before general marginal before spin: reproducing known moment polytopes validates the backend abstraction first.
- Fermionic hard-carrier identifiability before chemistry: do not spend time integrating chemistry until a real structural advantage survives.
- Generator spread-bound and determinant-stage work may proceed partly in parallel internally, but rank-11 completeness depends on both.

## Shared handoff protocol
Every agent should finish with:
- HEAD SHA used;
- exact files changed/created;
- claims table with PROVED / EXACT COMPUTATION / EMPIRICAL / OPEN;
- strongest positive result;
- strongest falsification/negative result;
- next three experiments ranked by information gain;
- any cross-track result another agent should consume.

Agents must not silently promote another track's empirical result into an assumption.

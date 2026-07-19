# PROVENANCE
Constraint systems ranks 9-10: parsed from M. Altunbulak, PhD thesis (Bilkent, 2008), Appendix A;
counts validated against Altunbulak-Klyachko CMP 282, 287 (2008): 52/60/93/125/161.
(3,9): one inequality lost at thesis page boundary; repaired (recovered facet:
3l1+2l2-l3+l5+l7+2l8-l9 <= 4); repaired system: 52 inequalities, 58 vertices.
Vertices: exact enumeration (lrslib, rational arithmetic). Validation: embedding coherence
(rank r -> r+1), frozen-core lifts, particle-hole self-duality of (4,8) and (5,10) at
inequality/vertex/verdict level, agreement with published rank-8 vertex tables and (4,9)=103.
Census verdicts: two-stage solver certificates (CP-SAT integer stage, ortools 9.15.6755;
CBC real stage). (3,10) vertex 89 resolved with extended real-stage budget (7200 s).
Rank-10 constraint lists are conjecturally complete per AK2008; rank-10 rows conditional.
Verdicts transport to dual systems (4,7),(5,8),(5,9),(6,9),(6,10),(7,10) by particle-hole
complement bijection on determinants.

## states.jsonl (constructed extremal states)

Per-vertex extremal states, produced by the routed state pipeline
(scripts/solve_all.py --all): DESIGN-INT vertices built directly from their
classification witness (exact by construction), DESIGN-REAL and INTERFERENCE
vertices through the weights-first block/clique solver (max_clique=3), certified
by exact characteristic-polynomial identity where the closed form is recognized
(status EXACT), else labeled TIER-C (exact numeric state, closed form not yet
recognized) or FAIL (outside the current ansatz family). Each record carries the
classification verdict, so the dataset cross-checks census_master.csv.

This file is PARTIAL and resumable: it currently covers (3,9) and (4,9), still in
progress on the slow k>=3 interference vertices; scripts/solve_all.py skips
already-recorded (system, index) pairs, so rerunning resumes and extends it to
(3,10), (4,10), (5,10). Regenerated fresh it is deterministic up to the
interference phase gauge, which does not affect the certified spectrum.

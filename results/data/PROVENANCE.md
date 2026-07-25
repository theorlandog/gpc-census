# PROVENANCE
Constraint systems ranks 9-10: parsed from M. Altunbulak, PhD thesis (Bilkent, 2008), Appendix A;
counts validated against Altunbulak-Klyachko CMP 282, 287 (2008): 52/60/93/125/161.
(3,9): 52 inequalities, 58 vertices.
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

This file is COMPLETE: all 799 vertices of the nine determinate systems
(ranks 6 to 10) carry a certified closed-form extremal state. Regenerated fresh
it is deterministic up to the interference phase gauge, which does not affect
the certified spectrum.

Correction (v0.11): three DESIGN-REAL records, (3,9) index 53, (3,10) index 106,
and (4,10) index 53, shipped valid extremal states on one-hop-connected supports
rather than genuine design witnesses. They were replaced by design witnesses
transported (padding / frozen-core lift) from the (3,8) parent design at
denominator 12, each verified by the exact characteristic-polynomial identity;
their feature_table.csv rows (one-hop-class and max-class-size counts) were
regenerated to match.

Update (v0.12): the (4,9) index 65 record (v_B = (20,14,14,14,14,4,4,4,4)/23)
now ships the REAL extremal state: the theta=0 endpoint of the wall family
k(t) = (1+t, 8-t, 4, 3, 2-t, 2+t, 1, 2)/23 on the same support, at
t0 = 7/17 - 18*sqrt(35)/85 (reality wall 85 t^2 - 70 t - 119), squared
weights in Q(sqrt(35)). Generated and certified by scripts/vb_real_state.py
(exact characteristic-polynomial identity; independently re-verified by
scripts/verify_states_standalone.py). closed_form.weights for this record
holds exact radical strings instead of integers; its feature_table.csv
all_real flag is now 1. The previously shipped phase-carrying fiber point is
the t=0 rational anchor of the same family, recoverable from the family
metadata in the record.

# CONSTRAINT GENERATOR -- architecture v0 and rig session spec (Session 4)

## Why now
This week accidentally built the generator's decision layer: the contraction
attack is a fast attainability + DISTANCE oracle (Hoffman-Wielandt: the residual
floor equals squared distance from a spectrum to the moment polytope), and the
SOS pipeline (scripts/sos_nonattain.py, symmetry reduction validated) is a
rigorous non-attainability certifier. Facet enumeration therefore no longer
needs Schubert calculus as a PREREQUISITE: the polytope can be reconstructed
from oracles, with Schubert (upstream's character-table work) demoted to the
final validity-proof layer. This is the architecture that makes rank 12+
tables producible.

## Architecture (oracle-driven reconstruction)

Layer 0 -- relation detection (WORKS, validated at (3,6)): sample spectra,
  SVD null space of [lambda | 1], lattice-reduce to integer relations. Blind
  run at (3,6) recovered exactly the 3 particle-hole dualities (span verified).
Layer 1 -- support oracle (THE core build): h(c) = max over states of
  <c, sorted-spec(gamma(Psi))>, computed by Wirtinger L-BFGS with multi-start
  (nonsmooth at spectrum degeneracies -- use soft-sort smoothing or subgradient
  handling; validate against known values: chamber directions give Slater
  subset sums exactly).
Layer 2 -- hull-by-oracle: beneath-beyond / incremental convex hull driven by
  support queries (NOT sampling: gate-1 showed cloud hulls never reach the
  boundary -- random states are interior-biased). Output: candidate facet
  normals + offsets.
Layer 3 -- rationalization: facet normals via PSLQ/limit_denominator with
  cross-validation (re-query support in the rationalized direction; the gap
  must be < 1e-12).
Layer 4 -- certification, per object, labeled honestly:
  - V-description: every claimed vertex gets an exact certified state (the
    census pipeline, unchanged) -> CERTIFIED.
  - H-description validity: for each candidate facet, SOS-certify
    non-attainability of exterior probe points along it (evidence), and where
    feasible prove validity via the Schubert/character-table verifier ->
    CONJECTURED+SOS-TESTED, upgraded to PROVED per facet as the verifier lands.
  - Completeness: cross-check V against H (vertex enumeration of the
    candidate H-description must reproduce the certified V exactly).

## Gates (validation-first; no new-rank claims until all pass)
G1: (3,6) -- relations DONE; facets: blind-reconstruct the known table via
    Layers 1-3, diff against src/gpc_census constraint tables. EXACT match.
G2: (3,7) and (3,8) blind -- exact match against A&K published tables.
G3: (3,9)/(3,10) blind -- match the census's own H-data (larger, known).
G4: (3,11) -- reproduce the 19-vertex bracket INDEPENDENTLY: the generator
    must find the level-<=4 facets + confirm the four excluded candidates sit
    outside (cross-check with the SOS deltas ~ dist^2).
G5: FRONTIER -- (3,12) then rank 12+: new tables, labeled per Layer 4.

## Notes
- Reuse: scripts/contraction_attack.py (oracle core), sos_nonattain.py /
  sos_symmetry.py (certifier), upstream Murnaghan-Nakayama tables (validity).
- The support oracle at degenerate spectra is the main numerical-analysis
  risk; budget real effort there (soft-sort epsilon-scheduling, many starts).
- Deliverable framing for the papers/EV: "certified V-description +
  SOS-tested H-description" is the honest product; chemists can use
  conjectured facets with stated status -- that is exactly how the rank-10
  list was used by the field for 15 years, minus the honesty labels.

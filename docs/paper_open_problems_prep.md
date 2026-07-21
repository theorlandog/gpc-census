# Prep note: open problems the census closes (material for the final paper edit)

NOT part of the paper. Staging ground for the single consolidated edit of
results/report/main.md: exact source quotes, verified data points, and the tier
structure, so the framing is ready to drop in. Final residual count is 785
certified / 14 open (rig sweep at max_card 44, max_clique 4, max_cliques 0
cracked none of the 14; the off-clique cancellation gap is not a walltime issue).

## Verified data anchors (checked against results/data at this HEAD)

- Design/interference onset: (3,6) and (3,7) are pure design (0 interference);
  (3,8) is the first system with interference (11 of 38). Matches the thesis
  claim that the Dadok-Kac construction first fails at wedge^3 H_8.
- Per-system design / interference: (3,6) 4/0, (3,7) 10/0, (3,8) 27/11,
  (4,8) 22/0, (3,9) 38/20, (4,9) 87/16, (3,10) 71/42, (4,10) 134/25,
  (5,10) 250/42.
- Certified: 785 of 799. Uncertified 14, all Tier-A SOLVE-FAIL: (4,9) 2,
  (3,10) 8, (4,10) 2, (5,10) 2 (12 of the 14 are rank-10).

## Tier 1: explicitly posed, closed by shipped artifacts

1.1 v_A / v_B (the founding problem). Altunbulak-Klyachko, Comm. Math. Phys.
    282, 287 (2008), Sec. 6.2.2, verbatim: "The remaining two vertices
    [16,16,16,6,6,6,6,6,6]/21 and [20,14,14,14,14,4,4,4,4]/23 were checked only
    numerically." Census: v_A DESIGN-INT (explicit integer design), v_B
    INTERFERENCE with cos gamma = 3/(4 sqrt(14)), both exact. This is the
    headline; frame it as closing that sentence. Cite [@AK2008].

1.2 The exact domain of the Borland-Dennis / Dadok-Kac construction (PROMOTE to
    a named result). Altunbulak thesis, Sec. 3.6, verbatim: "for a fixed support
    T the set of unordered spectra ... form a convex polytope. It is not known
    when this approach gives the whole moment polytope. The smallest fermionic
    system where it fails is wedge^3 H_8." The disconnected-support / diagonal
    1-RDM construction IS the design mechanism (one-hop-free support => diagonal
    1-RDM), so the trichotomy answers this with per-vertex certificates: the
    construction reaches exactly DESIGN-INT and DESIGN-REAL, provably misses
    every INTERFERENCE vertex, first fails at (3,8), failure fraction stabilizes
    at 34-37% (N=3), 16% (N=4), 14% (N=5) across ranks 9-10. A 46-year-old
    heuristic (Borland-Dennis 1972 [@BD1972]) gets its exact validity domain.
    Refs: thesis [@AltunbulakThesis], Mueller J. Phys. A 32, 4139 (1999) (add
    if used).

1.3 The lost rank-9/10 extremal-state tables. A&K print states only for ranks
    <= 8; rank-9 data was electronic supplementary (now dead), rank-10 never
    released (see RESEARCH.md). Census: results/data/states.jsonl is the first
    public, certified, machine-readable extremal-state library for ranks 6-10,
    a strict superset of everything A&K distributed. Argues for a Zenodo DOI on
    results/data at submission.

## Tier 2: closed modulo the 14 (state precisely)

2.1 First rigorous completeness for rank 9/10. A&K Sec. 6.2.2: "Only for
    smallest system wedge^3 H_9 we have a rigorous justification of
    completeness." Thesis Ch. 6 (rank 10): "the constraints are complete only
    beyond a reasonable doubt ... some [vertices] still evaded all the efforts.
    For the latter we resort to the numerical optimization." Convexity closure:
    certified attainment of every outer-polytope vertex forces P_out = P.
    - wedge^4 H_9: census certifies 101/103; the 2 uncertified lie inside A&K's
      rigorously-proved 101, so census + A&K = complete. Say it that way.
    - rank 10: census gives the first published, independently checkable
      attainability certificates for 552/564 vertices, reducing an 18-year
      "beyond reasonable doubt" gap to an explicit 12-vertex list, each
      compute-bound (Tier-A SOLVE-FAIL), not open mathematics.

## Tier 3: one-script queries of the census (immediate-applications section)

3.1 Maciazek-Tsanov doubly-excited coverage. J. Phys. A 50, 465304 (2017)
    [@MT17]: inner bound exact for (3, d<=7); where it is exact was left open.
    Census: exactly on the interference-free systems; the exceptional-vertex
    list per system reads off the trichotomy. Collaboration hook.
3.2 Slater-sparsity at the boundary. Chakraborty-Mazziotti, J. Chem. Phys. 148,
    054106 (2018) [@CM2018]. Census: states.jsonl has exact support, weights,
    phase type per vertex; the real-suffices vs complex-forced dichotomy IS the
    trichotomy; two-exchange-channel bound is a new structural law.
3.3 Experimental target states. Hackl-Li-Akopian-Christandl, Phys. Rev. A 108,
    012208 (2023) [@Hackl2023]. Census: certified closed-form states are the
    prep targets; Slater ranks bound the entanglement depth. Theory input they
    lacked, now supplied (experiments still to be done).
3.4 Exact polytope invariants. Reuvers, J. Math. Phys. 62, 032204 (2021)
    [@Reuvers2021]. Census: exact V-representations make volumes, f-vectors,
    facet incidence a finite computation; see docs/polytope_invariants.json and
    scripts/polytope_invariants.py (exact dim / vertices / facets / simple
    vertices per system).

## Tier 4: enabled, not solved (outlook only)

RDMFT exact functionals (Schilling-Schilling, Phys. Rev. Lett. 122, 013001
(2019)); MCSCF face hierarchy (NJP 22, 023001/023002 (2020)); fermionic
entanglement polytopes / SLOCC [@WDGC2013]. Each needs its own derivation.
EXCLUDE: spin-adapted magnetism questions (different lattice; note as the
natural H_nu extension).

## Placement plan

1.1 headline (done). 1.2 a named result with the thesis Sec. 3.6 quote as
epigraph. 1.3 and 2.1 in the introduction contribution list with the exact A&K
quotes. Tier 3 a compact "immediate applications" section (one paragraph each,
pointer to the emitting script). Tier 4 outlook. Pose the 12 remaining rank-10
vertices as the paper's own explicit open problem.

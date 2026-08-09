# Agent Track: Anyon and Newton-Polytope Scout

## Mission
Perform cheap falsification-oriented scouting on two more distant directions. This is intentionally a lower-priority track and should not consume major implementation time unless a strong wedge appears.

The topics are unrelated scientifically but both are exploratory enough to share one serial scout agent rather than occupy separate agents.

## Part A: anyons / fusion constraints
Investigate whether fusion rules, topological charge sectors, and braid/fusion measurement structure produce a qvalidate-style problem with a tiny but high-value buyer set.

Questions:
- Can partial fusion measurements certify a logical/topological state more cheaply than standard protocols?
- Are there convex/polyhedral structures genuinely analogous to the moment-polytope machinery, or would this require a completely different engine?
- Is there a commercially expensive validation problem at current topological-QC companies/labs?

Produce one small exact toy benchmark if feasible. Stop if the geometry does not align or the buyer/use case is too speculative.

## Part B: Newton polytopes / sparse polynomial systems
Investigate whether any of the new generator machinery (profiles, cocircuits, Chow-like shadows, symmetry-native enumeration) transfers to a commercially meaningful sparse-polynomial problem.

Candidate use cases:
- mixed-volume / homotopy path reduction;
- sparse elimination;
- polynomial optimization preprocessing;
- symbolic computation.

Compare against mature tools and algorithms. The bar is high: only continue if there is a clear capability gap and buyer class.

## Deliverables
- `docs/anyon_newton_scout.md`
- `results/data/anyon_newton_scout.json`
- one-page go/no-go verdict for each domain

## Guardrails
- No broad platform implementation.
- No novelty claims without current literature review.
- Kill quickly if the connection is merely aesthetic.

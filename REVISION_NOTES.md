# Revision notes, main paper

## v0.11 → v0.12 (July 2026): J. Phys. A submission pass

Scope: journal-style conversion, bibliography audit against publisher
records, one scientific correction, and consistency wiring. No data files
changed; SHA256SUMS untouched and re-verified.

**Style (IOP / J. Phys. A: Math. Theor.):**
- Citation style switched from the vendored APS CSL to the vendored IOP
  numeric CSL (`institute-of-physics-numeric.csl`); Makefile updated; the
  APS CSL removed. Rendering validated with pinned pandoc 3.6.4 +
  pandoc-crossref v0.3.18.2: zero unresolved citations or references.
- Keywords line added to the abstract; date line marks the target journal.
- `nocite: "@*"` removed: every bibliography entry is now cited in the
  text, so the reference list contains exactly the cited works.
- "Comm. Math. Phys." corrected to "Commun. Math. Phys." (abstract and
  bib); "Amer." to "Am."; the single UK spelling ("neighbouring")
  normalized to the manuscript's US convention.
- Remark numbering restored to document order: the census-section remark
  is now Remark 3 and the realification remark (rem:orbit) is Remark 4;
  all cross-references updated.

**Scientific corrections:**
- The abstract no longer implies the v_B phase is forced at the vertex:
  it now states the phase is forced only within the integer-weight
  single-block class, and that the vertex admits real extremal states at
  quadratic-irrational weights (matching Theorem 4 and Remark 4).
- The denominator-34 closure vertex is no longer called "first-order
  rigid" without qualification. Per the silence-rank lemma
  (docs/silence_rank_lemma.md), full rank of the five silence conditions
  on the incidence kernel is a fixed-1-RDM statement; the fixed-spectrum
  fiber is positive-dimensional because cross-block channels impose no
  first-order constraint. The text now says exactly that.

**Bibliography audit (every entry verified against publisher records):**
- Hackl2023: wrong title and initials corrected to L. Hackl, D. Li,
  N. Akopian, M. Christandl, "Experimental proposal to probe the extended
  Pauli principle," Phys. Rev. A 108, 012208 (2023).
- Reuvers2021: subtitle corrected to "the Pauli principle dominates"
  (was "the low-density limit"); arXiv:1911.00471 added.
- vdBerg25: updated to the published STOC 2025 version (pp. 756-765,
  proceedings title), arXiv:2510.08336 retained.
- Castillo2021: updated from arXiv preprint to the published version,
  Ann. Henri Poincare 24, 2241-2321 (2023).
- LPW2020: updated from arXiv preprint to the published version,
  Int. Math. Res. Not. 2023 (19), 16778-16836.
- MS20: volume 22, article 023002 added; title capitalization fixed.
- Klyachko2006: page range completed (72-86). BCMW17, W92, S98,
  Kirwan1984, Hardt1980, SL91, Ressayre, WDGC2013, BD1972: page ranges
  completed. Klyachko2009 retyped as @misc (arXiv-only, never published).
- Added (both now cited in the text): Dadok and Kac, "Polar
  representations," J. Algebra 92, 504-524 (1985), cited at contribution
  5 for the Borland-Dennis/Dadok-Kac construction alongside BD1972; and
  Cox, Little, Schenck, "Toric Varieties" (AMS GSM 124, 2011), cited in
  Proposition 2 for the simplicial toric quotient (closing the v0.11
  "remaining actions" item 2).

**Consistency wiring (closing v0.11 "remaining actions" item 6):**
- New tests/test_manuscript_consistency.py runs
  scripts/check_manuscript_counts.py and re-verifies SHA256SUMS inside
  the default pytest suite, so `make test` and CI now fail on any
  manuscript/data/checksum drift. Full suite green (78 passed); the
  standalone checker re-passes all 799 certificates.

# Revision notes — main paper v0.10 → v0.11

Commit reviewed: 7543ca9 (same tip the referee reviewed). Every claim below was
re-verified computationally against the released artifacts; scripts are included.

## Independent verification performed

**Confirmed correct (exact recomputation):**

- Theorem 1 (ψ_A): 1-RDM rebuilt symbolically with independent fermionic-sign
  bookkeeping; off-diagonals vanish identically; spectrum exactly v_A.
- Theorem 3 (ψ_B): only nonzero off-diagonal is ρ₈₉ = (11+i√215)/92; exact
  spectrum v_B; block eigenvalues 14/23, 4/23.
- Library state for v_B (states.jsonl): exact characteristic-polynomial
  identity passes after `expand_complex` hardening.
- Census table: all nine rows of Table 1 match states.jsonl exactly
  (631 INT + 12 REAL + 156 INTERFERENCE = 799; rank-10 total 564; interference
  fractions 34→37% / 16% / 14%; the full (4,9) interference vertex list;
  both (4,9) DESIGN-REAL vertices).
- README table (27 designs at (3,8), 643/156/799 totals) matches the dataset.
- Toric claims at v_B: 8 active facets, simplicial cone, primitive-ray SNF
  diag(1,23,…,23) over L={u∈Z⁹: Σu=0}, index 23⁷; dual facet-normal index
  exactly 23, cyclic; every nonzero pairing ⟨W_j,r_k⟩ = 23. All verified in
  exact arithmetic (now formalized as Proposition 2 per referee item 7).
- Remark 3 real state: t₀ = 55/109 − (42/109)√15 gives cos θ = 1 exactly,
  interior to positivity. The second root (55+42√15)/109 ≈ 1.9969 is also
  interior and also θ=0 (added to the text).
- Conjecture 2 norm identities: 11²−21=10², 817²−1633=816²,
  2651²−48²·2769=805² all check.
- Minimal polynomial 32x⁴+55x²+32 and field Q(√2,√−119) for cos = 3√2/16;
  derived the analogous 56x⁴+103x²+56 and Q(√14,√−215) for ψ_B's holonomy.
- Loop census: 93 of 156 interference supports are loop-free; 63 carry 82
  independent loops; closure vertices have kernel dims 1 (den-26) and 5
  (den-34), matching the text.
- All 799 closed forms: standalone checker (see below).

**Errors found beyond the referee report:**

1. §galois attributed the kernel vector (1,−1,0,0,−1,1,0,0) and the holonomy
   arctan(√119/3) (cos 3√2/16) to "the eight support determinants of
   Theorem 3." Wrong: Theorem 3's support has kernel (0,1,−1,0,−1,1,0,0) and
   holonomy cosine 3/(4√14) = cos γ. The quoted kernel and holonomy belong to
   the *library* state for v_B, a different fiber point (support
   {1235,1239,1246,1347,1458,1489,2479,3459}, weights (1,8,4,3,2,2,1,2)/23).
   Both verified by direct computation. The revision presents holonomy
   explicitly as a per-state invariant that varies along the fiber and gives
   both states' data.
2. §posgeo's "wall polynomial 85t²−70t−119" is the Cayley–Menger locus of the
   *library-state* family, not of Remark 3's Theorem-3 family, whose wall
   polynomial is 109t²−110t−215 (roots (55±42√15)/109, matching Remark 3's
   Q(√15) weights). Both families and both polynomials verified; the revision
   states both and labels which family each belongs to.
3. Proposition 1's "exotic vertex (5,5,5,5,2,2,2,2)/28" is wrong in the
   paper's own normalization (sums to 1, not 3). The dataset's integer form is
   (15,15,15,15,6,6,6,6)/28. Corrected.
4. The (3,8) DESIGN-REAL vertex is (3,3,3,3,3,1,1,1)/6 with a rational design
   at denominator 12; now stated explicitly in Prop. 1.

## Changes mapped to the referee's blocking items

1. **Census contradictions** — all counts regenerated from states.jsonl:
   abstract and Prop. 1 now say 27/38 with 26 integer; contribution 2 says
   27/38. Provenance note added explaining the earlier 22 (shallower search).
   New CI test `scripts/check_manuscript_counts.py` (added) regenerates every
   number in the paper from the dataset and fails on drift. It passes.
2. **Interference accounting** — replaced "the remaining 142" with the
   referee's recommended unambiguous decomposition (643 = 631+12 design;
   156 = 142 engine + 12 audit + 2 fiber-sparsified), displayed as a block
   quote in §census910 and enforced by the consistency test.
3. **Rank-10 status** — new Theorem 5 with the validity hypothesis explicit
   and the four-step convexity proof; all "conjecturally complete /
   conditional" language removed elsewhere; table caption points to Thm 5.
4. **Interference misdescription** — every "forced off-diagonal 1-RDM"
   formulation replaced with the termwise-vs-coherent-cancellation definition
   (Remark 4 rewritten; §framework language removed with the section).
5. **Fiber stratification** — the "Inverse geometry" section is deleted as a
   section; its rigorous content (Hardt existence, semialgebraically definable
   invariant, explicit non-claims about orbit-type/symplectic/Hardt strata)
   moved to §outlook with the referee's recommended wording.
6. **Certificate types** — [E]/[CA]/[N] classes defined in the introduction,
   attached to every theorem; certificate matrix added as Table 2 in a new
   §Certificate architecture; Thm 2's proof now states explicitly that
   two-solver agreement is not a distributable proof object and that an
   exportable Farkas/branch certificate is planned. Abstract no longer says
   "solver-independent certificates."
7. **Toric claim** — formal Proposition 2 with ambient lattice, simpliciality,
   ray data, SNF, dual index, pairing values, canonicity caveat (lattice
   equivalence; pseudoreflection question left open).
8. **Positive geometry** — "one positive geometry" removed; explicitly framed
   as spectrahedral feasibility geometry with a sentence disclaiming any
   ABL-canonical-form construction; adjoint degree F−dim−1 stated as an
   observation on the computed polytopes with no general claim.
9. **Universal claims** — "for every such family, not only v_B" and the
   universal density–density pinning claim removed; both restated as verified
   for the two computed v_B families with the general statements labeled
   conjectures and the failure modes the referee listed named explicitly.
10. **Scaling resources** — vCPU-month and wall-clock claims removed; deferred
    to the companion note with the required benchmark metadata; abstract's
    final sentence replaced per the referee's wording.

## Softened claims (referee's second list)

- Title: "a complete census" → "a certified census."
- "Exact domain of validity" of Borland–Dennis/Dadok–Kac → "vertex-level
  reach within the census through rank 10" (contribution 5).
- "The vertices that resisted AK are precisely of this interference type"
  removed (contradicted by v_A being design); replaced with the correct
  statement.
- "First genuinely complex extremal state known" → scoped novelty claim with
  an explicit statement of what was and was not surveyed.
- Numerical non-realifiability removed from the abstract; Remark 3 now splits
  [E] and [N] content and uses "extensive numerical optimization found no
  realifying orbital rotation."
- Theorem 4 restated with its three restrictions inseparable from the claim.
- Terminology block added (rank, published system, natural denominator);
  DESIGN-INT/REAL/INTERFERENCE naming caveat added at Definition 1 (labels
  kept to match the released dataset; renaming would desynchronize paper and
  artifact).
- Abstract cut to ~45% of its previous length, restructured to the referee's
  five points.
- Speculative material (Hamiltonian pilot, toric-fiber conjecture, rank-11/12,
  scaling program) consolidated in §outlook with per-item status labels and a
  citation to the planned companion research-agenda note, per the layered-
  release strategy.

## New files

- `scripts/check_manuscript_counts.py` — manuscript–dataset consistency test
  (passes at HEAD; add to CI).
- `scripts/verify_states_standalone.py` — independent exact certificate
  checker for all 799 closed forms (no shared code with the pipeline).
- `results/report/main_revised.md` — the v0.11 draft.

## Remaining before arXiv (author actions)

1. Mint the companion-note DOI (or strike the two forward references to it).
2. Add a toric-geometry reference (e.g. Cox–Little–Schenck 2011) to
   references.bib for Proposition 2's standard quotient description.
3. Decide the Thm 5 hypothesis question: if AK2008's derivation establishes
   validity of each listed rank-10 inequality (as opposed to completeness),
   say so with a pinpoint citation and the theorem becomes unconditional.
4. Produce the exportable MILP proof object for Theorem 2, or keep the
   planned-work language.
5. Release the holonomy minimal-polynomial table (24 polynomials) in
   machine-readable form and extend the Galois tally to the 12 audit-recovered
   states, or keep the stated scope.
6. Wire both new scripts into CI (Makefile target + GitHub Actions).
7. Regenerate the Zenodo deposit after the text changes so the DOI-stamped
   PDF matches the dataset.

## Post-verification addendum: dataset finding (new, found by the checker)

Full standalone verification of all 799 closed forms: **796 pass all checks;
3 pass the spectral certificate but fail the design-support check**:
(3,9) idx 53, (3,10) idx 106, (4,10) idx 53 — the padding/frozen-core family
of the (3,8) DESIGN-REAL vertex (3,3,3,3,3,1,1,1)/6. Their shipped states are
valid extremal states (exact char-poly identity holds) but sit on one-hop-
connected supports with nonzero off-diagonals — i.e. they are not design
witnesses, contradicting the manuscript sentence "the 12 real designs from a
rational weighted design on a one-hop-free support." The classification
labels are correct (design witnesses exist); only the shipped states are the
wrong kind. Correct replacements are one transport away from the (3,8)
parent's genuine design (weights (2,1,3,3,1,2)/12): padding for (3,9)/(3,10),
frozen-core lift for (4,10). All three replacements are constructed and
verified exactly in `patches/design_real_transport_fix.jsonl`; regenerate
them through the pipeline (for provenance + SHA256SUMS) before release, or
apply the patch and update checksums.

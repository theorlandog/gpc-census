# Research state and continuity notes

This is the living research document. Read it before touching the science.
House rules live in AGENTS.md. The full evidence log of the census era
(pre-registration scorecards, campaign session notes, retractions with
mechanism) is frozen in `docs/RESEARCH_ARCHIVE.md`; status tags below point
into it by section name. The archive is the record; this file is the map.

House rule inherited from the archive (non-negotiable): every 🟦/❌ status
tag carries a pointer to its certifying artifact (script, jsonl, doc, or
archive section). A claim with no in-repo artifact cannot carry an evidence
tag.

Status legend:

- ✅ Theorem (proved in project or standard reference)
- 🟨 Theorem candidate (appears to follow from standard mathematics but needs a written proof)
- 🟦 Confirmed by certified project computations
- ❌ Disproved by certified project evidence
- 🔬 Open
- 💡 Speculative

# Part I. The census program (established base)

## What this project is

Complete classification program for fermionic moment polytopes (generalized
Pauli constraints, GPCs). The pure-state N-representability problem: which
one-body occupation spectra arise from N-fermion pure states in d orbitals.
Altunbulak and Klyachko (CMP 282, 287 (2008)) computed constraint systems
through rank 10 and left two vertices of the wedge^4 H_9 polytope resolved
only numerically. This program closed both (paper in results/report), then
classified every vertex of every known system, and now extends the frontier.

## The trichotomy (introduced by this program)

Every polytope vertex is exactly one of:

- DESIGN-INT: an integer weighted design exists at the natural denominator
  (nonnegative integer determinant weights, prescribed mode sums, support
  one-hop free, meaning no two support determinants share N-1 modes).
- DESIGN-REAL: a real nonnegative design exists but no integer one.
- INTERFERENCE: neither exists; complex phase cancellation is forced.

Feasible design verdicts carry exact witnesses. The 11 negative verdicts in
Proposition 1 and the v_B verdict in Theorem 2 now carry exact global
disjunctive Farkas/branch certificates: 192 primitive integer vectors, 5141
symmetry-orbit hitting clauses, and 12 branch DAGs verified without a solver.
The other 144 INTERFERENCE verdicts still rest on the pinned two-solver census
and structural transport checks, so their verdict-level status remains
computer-assisted. The one-hop exclusion is essential; a degree-sums-only
model is strictly weaker and misclassifies (this bug was made and caught once;
see tests pinning v_A and v_B).

## Key results encoded in results/data

- Census complete for every determinate system, ranks 6 through 10
  (799 vertices; duals covered by particle-hole verdict transport).
- Certified closed-form extremal states for ALL 799 vertices (ledger
  regenerated, CI-green): every design vertex from its witness, every
  interference vertex through the routed pipeline (per-phase recognition,
  the constructive off-diagonal-target exactifier of stage 3b, PSLQ),
  14 closed by state transport (scripts/transport_states.py), and the two
  rank-10 cancellation-regime closures v89 and v103, found off the rational
  search grid after the attainability audit self-corrected
  (RESEARCH_ARCHIVE "RANK-10 CLOSED").
- v_A (16,16,16,6,6,6,6,6,6)/21 is DESIGN-INT; v_B (20,14,14,14,14,4,4,4,4)/23
  is INTERFERENCE with cos(gamma) = 3/(4*sqrt(14)). The v_B phase is not
  special: it is the k=2 instance of the off-diagonal-target mechanism that
  fixes every interference phase.
- Interference onset: absent at rank 8 in the N=4 series, first appears at
  rank 9. Fractions stabilize across ranks 9 to 10 within each series
  (34-37 percent at N=3, 16 at N=4, 14 at N=5).
- Functorial structure: trailing-zero padding and frozen-core lifts
  (lambda -> (1, lambda)) preserve verdicts and generate most interference
  from lower-rank originals; at (4,10) they generate all of it.
- LEDGER NOTE (2026-07): all 12 DESIGN-REAL states ship at
  state-den/spectrum-den ratio 2 (the archive's mined "9 of 12" figure is a
  stale snapshot; the ratio is representative-dependent, exactly as the
  P2/P3/P4 caveat predicts). The census min-support column (support of the
  certified state) is an UPPER bound on s_Q^NO for the 645 records that pass
  the diagonality gate, and on s_Q^free only for the other 154; see the
  support taxonomy below before quoting it, and use
  results/data/states_natural_orbital.jsonl for an s_Q^NO bound at an
  interference vertex. It is NOT known strict at v103.

## Cross-realization invariants (2026-08-07)

Every census object recomputed in three presentations other than the
combinatorial one it shipped in, all exact, all 799 rows for the first two.
[docs/symplectic_invariants.md; docs/prereg_symplectic_invariants.md;
results/data/symplectic_invariants.json;
results/data/quantization_multiplicities.json]

- 🟦 PROVED on the census, 799 of 799, zero violations: the SUPPORT DETERMINES
  THE TORAL SYMMETRY. With `A` the 0/1 support-incidence matrix,
  `dim(Stab_{U(d)}([psi]) intersect u(1)^d) = d - rank([A_S - A_S0]_S)`, a
  number with no amplitude input. This is the one clean identification between
  a combinatorial and a symplectic invariant that the campaign found.
- ❗ AND IT STOPS THERE. The support does NOT determine the full stabilizer: 67
  of 411 support-invariant classes carry more than one `dim Stab`, and
  `dim Stab` exceeds its toral part at 445 of 799 vertices, by up to 16. The
  boundary between the two is the sharpest open target this campaign produced.
- 🟦 The moment polytopes are TORIC-SINGULAR above rank 6. Exactly 183 vertices
  have a simplicial tangent cone, matching the simple-vertex column of
  docs/polytope_invariants.json system by system, and only 3 have a UNIMODULAR
  one, all three in (3,6). Lattice indices reach `37^8` at (5,10)#144.
  `index = q^(dim-1)` holds at only 29 of 183, so it is not a law.
- 💡 The census minimum `dim Stab = 2` is attained at five vertices, two of
  which are v89 and v103, the rank-10 cancellation-regime closures that carry
  the fiber program's hardest open questions. An invariant that knows nothing
  about how they were found picks them out.
- ❗ QUANTIZATION PERIOD, a new invariant, and a scored FAIL that is the
  finding. The multiplicity of `k*lambda` in `Sym^k(wedge^n C^d)` is 1 on
  multiples of a period `p` and 0 elsewhere. `p = q` at ten of the fourteen
  (3,6) and (3,7) vertices and `p = 2q` at the other four, never anything else,
  with zero exceptions to the periodic model. Every nonzero multiplicity is
  exactly 1, so the reduced space over a census vertex is a point as far as
  multiplicity growth sees; what the polytope fails to record is at which `k`
  the vertex is realized. (3,6)#3 is the classical case: `lambda = (1/2)^6`
  makes `k*lambda` a determinant power, and `C[wedge^3 C^6]^{SL(6)}` is
  generated by a quartic, so `p = 4`. The mechanism at (3,7)#5 and (3,7)#8 is
  OPEN.
- 🟦 THE PERIOD IS EXPLAINED (2026-08-07, follow-up campaign). The quantization
  period is exactly the order of a finite isotropy character of the certified
  attainer: with `omega(h) = c(h) / prod_i det(h_i)^{lambda_i}` for a
  `lambda`-preserving `h` fixing the state's line, a nonzero multiplicity at
  `k` forces `omega(h)^k = 1`. Predicted period equals measured period at
  **14 of 14** vertices, with all four doubling cases carried by one explicit
  permutation each. `omega` is trivial on the identity component, so the whole
  obstruction lives on `pi_0`, which is why the continuous stabilizer
  dimensions of the first campaign could not see it. (3,7)#7's 18 obstructing
  elements are exactly (3,6)#3's 18 padded with a fixed orbital, so the honest
  count is 3 independent doubling instances. The mechanism at (3,7)#5 and
  (3,7)#8, previously OPEN, is now explicit.
  [docs/quantization_character.md; docs/prereg_quantization_character.md;
  results/data/quantization_characters.json]
- 🟦 AND THE OBSTRUCTION IS EXACTLY A Z/2, census wide. The character side
  needs only the attainer, so unlike the multiplicity side it is not capped at
  rank 7: all 799 predicted periods compute in 100 seconds with none
  unresolved, and `predicted_period / q` is 1 at 745 vertices and 2 at 54,
  never anything else, across all nine systems. ⚠️ Those 785 unverified rows
  are PREDICTED: divisibility is derived and holds at 799 under the stated
  semi-invariant hypothesis, equality is verified only at the 14 with a
  measured multiplicity sequence. Never quote them as measured periods; the
  artifact flags each row with `verified`.
- 🟦 THE SLICE FIXES THE VACUOUS KERNEL (2026-08-07). The symplectic slice
  `V_x = ker d(mu)_x / T_x(K_mu . x)` is the honest transverse object. At the
  uniform Borland-Dennis vertex (3,6)#3 the raw kernel is 19 but `dim V_x = 0`,
  which CERTIFIES the reduced germ is a point and the fixed-rho fiber near psi
  is exactly the `K_mu` orbit. ❗ And a nonzero slice does NOT mean a
  positive-dimensional germ: the Slater vertex has `dim V_x = 20` over a
  reduced space that is a single point, so the slice dimension is a local model
  INPUT, not a fiber dimension. ⚠️ REFUTED METHOD, recorded so it is not
  retried: `dim V_x - 2 dim(K_x . v)` at generic `v` in `V_x` returns -10 at
  the Slater vertex; the formula needs a generic point of `mu_V^{-1}(0)`.
  PARTIAL: v_B, v89 and v103 exceeded the exact-symbolic budget and are named
  in the artifact's `unfinished` field.
  [docs/symplectic_slice.md; docs/prereg_symplectic_slice.md;
  results/data/symplectic_slices.json]
- ⚠️ NOT CLAIMED: `dim M_lambda` from symplectic reduction. First-order tangent
  rank is vacuous at these ranks (`ker d(mu)` has dimension at least
  `2*binom(d,n) - d^2` for trivial reasons; the (3,6) Slater vertex, whose
  fiber is a single point, returns 20). Only multiplicity growth is used, and
  only at finitely many `k`. Every stabilizer number is UNREDUCED and
  FIXED-RHO.
- CONVENTION TRAP, now pinned by a test. `d(rho) = [rho, X^T]`, so the
  subalgebra preserving `rho` is the commutant of `rho^T`, not of `rho`. The
  two coincide for the 645 real-diagonal rows and differ for the other 154.
  The forced identity `stab_dim == commutant_stab_dim` failed at one row out of
  505 and is what exposed it; it is now clean on all 799, as is the independent
  check `commutant_dim == sum m_i^2` against the vertex spectrum.

## The support taxonomy (adopted 2026-07; do not collapse these)

Four different minima, with different lower bounds. Conflating them produced a
wrong reinstatement of a retracted bound, so every support claim in this
project must name which one it is about.

- s_M(lambda): least |S| with A_S q majorized by lambda. The airtight
  a-priori baseline.
- s_I(lambda): least |S| with the diagonal incidence system A_S p = lambda,
  p >= 0, sum p = 1 feasible. A DIAGONAL notion. s_M <= s_I.
- s_Q^NO(lambda): least support of a state whose 1-RDM is EXACTLY
  diag(lambda). s_I <= s_Q^NO. CORRECTION 2026-07: this is NOT what the
  library column bounds for every record. Only the 645 records that pass the
  diagonality gate (all 643 design states plus v89 and v103) have an exactly
  diagonal 1-RDM; the other 154 INTERFERENCE records attain the spectrum in
  another basis, so their support column bounds s_Q^free instead. Check
  results/data/attainer_gate_audit.json before quoting a support against
  s_Q^NO or s_I. [scripts/attainer_gate_audit.py; docs/rigidity_rationality.md]
- s_Q^free(lambda): least support of any state with spec(rho) = lambda.
  s_M <= s_Q^free <= s_Q^NO. **s_I does NOT bound s_Q^free**, so the two must
  never appear in one inequality.

Per-STATE (not per-vertex) there is a fifth: the ORBIT-MINIMAL support of psi,
the least support over the full U(d) orbit of that one state. It is a
basis-free complexity measure of the state, bounded below by s_M, and it upper
bounds s_Q^free for the vertex psi attains.

## The validation law (non-negotiable)

Every real error in this program was caught by structural invariants, never
by inspection: a sign-parse bug (caught by particle-hole self-duality of
(5,10)), a wrong classify port (caught by a Farkas feasibility
contradiction), a wrong solver gradient (caught by finite differences), and
the levi.py conjugation and fiber-conflation bugs (caught by the invariance
gate; RESEARCH_ARCHIVE "BUGS FIXED IN src/gpc_census/levi.py").
Consequences:

- No result ships without passing gpc_census.validate and the test suite.
- THE ATTAINER GATE (added 2026-07, enforced by
  gpc_census.validate.check_vertex_attainer). A state offered as attaining a
  vertex must have (1) a diagonal 1-RDM, (2) diagonal exactly lambda, (3) unit
  norm, and, if it is claimed to be a NEW attainer, (4) a 2-RDM spectrum
  DIFFERENT from the library representative's. Conditions 1 and 2 are what make
  its support comparable with s_I at all; a state passing only
  spec(rho) = lambda attains the vertex in some other basis and bounds
  s_Q^free instead. Condition 4 is basis-free and cheap: the 2-RDM spectrum is
  invariant under every one-body unitary, so a match means the state is the
  library one in rotated orbitals, which is a gauge statement and must be
  reported as one. This gate would have caught the mislabelled v103 support-10
  bound immediately.
- The preflight gate: any state-solving change must reconstruct v_B and
  certify its closed form end to end (scripts/solve_all.py --preflight)
  before campaigns run. The gate uses the weights-first solver
  (solve_vertex_exact_first with certify_tier_b): Tier A must attain the v_B
  spectrum exactly and Tier B must exactify to a certified closed form. It
  completes in minutes. The historical attain-based path
  (--legacy-preflight) is a Tier-A-only regression check, not the gate.
- Pin solvers. Census verdicts were produced under ortools 9.15.6755
  (a required runtime dependency); the CBC fallback labels its backend in
  every verdict.

## The falsification protocol (retained as method, not history)

Every new fiber law gets a pre-registered scored prediction before the
narrative is written, and a FAIL is a publishable finding. Mined-law
mortality ran about one in two, and both rank-10 closures falsified
pre-registered laws with mechanism. Protocol text: RESEARCH_ARCHIVE
"SKEPTICAL AUDIT + PRE-REGISTRATION". Final tally over the closure vertices
(RESEARCH_ARCHIVE "PRE-REGISTRATION SCORED: v103"): P1, P2, P5 pass; P3 and
P4 falsified with mechanism; P6 never in scope. The archive's P4 verdict
for v103 ("first-order rigid") is a fixed-rho statement and is annotated in
place there; the fixed-spectrum fiber at v103 is positive-dimensional (see
the silence-rank lemma below). P4's v89 verdict stands in both senses.

## Data-integrity findings (2026-08-05, sync-invariant cleanup)

Infrastructure, not science, but two of these were shipped defects.

1. [E] `census_master.csv` had NO generating script, was MALFORMED CSV, and was
   INCOMPLETE. The system field is written `(3,6)` and its comma was never
   quoted, so every row parsed as seven fields against a six-field header and
   any standard reader mis-split the file. It also carried 761 rows against 799
   vertices, missing the entire `(3,8)` block, with the pre-v0.11 DESIGN-REAL
   count of 11. Now generated by `scripts/build_census_master.py` from
   `states.jsonl`: 799 rows, 631 DESIGN-INT / 12 DESIGN-REAL / 156
   INTERFERENCE, reconciling with the README coverage table at 643 design and
   156 interference. Guarded by `tests/test_orphan_artifacts.py`.
2. [E] Completeness audit (`scripts/audit_data_completeness.py`,
   `tests/test_data_completeness.py`): every file under `results/data/` must
   carry a SHA256SUMS entry, a PROVENANCE row, a named generator and a guard
   test. 52 orphans on first run, 0 now.
3. OPEN `results/data/cyclic_window_benchmark.json` and
   `entangling_window_hamiltonian.json` have no surviving generator, no doc and
   no consumer anywhere in the repository. Declared in PROVENANCE.md, waived by
   name in the audit's `FAMILY_RULES`, and pinned by digest so they cannot
   drift. They should be regenerated from a recovered method or deleted; the
   decision is deliberately left open rather than taken incidentally.
4. [E] Doc tables that restate artifact data are now machine-emitted between
   `sync:NAME` markers by `scripts/emit_doc_tables.py`, checked by
   `tests/test_doc_sync.py`. Coverage: README coverage table, (3,11) bracket
   tally, and the four benchmark tables, whose hand-typed metrics are deleted.
   Benchmark footers carry the per-machine rule and the still-unmet
   equally-optimized-control caveat. Run `make emit`; never hand-edit a block.
5. [E] `docs/CLOSURE_PLAN.md`: the 797/799 status header is marked HISTORICAL
   (the ledger is 799/799) and the parallel SOS cand-44 item is marked
   SUPERSEDED, do not run, since the level-5 certificates settle that branch.

## Generator ladder gate 1: (3,7) is EXACT (2026-08-05)

[E] PROVED. (3,7) moves from RELATIVE to EXACT on the claim ladder, by the
independent inner/outer closure route in docs/fixed_n3_generator_methodology.md,
not by candidate enumeration. All seven pre-registered predictions PASS
(prereg frozen at 9c55f1b before the run).
Validity of every H-row gives P inside H; exact attainability of every
enumerated H-vertex gives H inside P by convexity; hence P = H.
Rechecked rather than read: all four Ressayre certificates replay from (h,z)
with weights, level sets, roots and tangent matrix rebuilt from scratch and
the determinant recomputed by Fraction elimination rather than Bareiss, all
four reproducing 48841, 1928920, 1455387, 4817604; each of the ten attaining
spectra recomputed from its own support and amplitudes rather than from
integer_form, the field under test.
Controls: the ordered trace-three slice alone has 13 vertices, so the four
rows cut it to 10, and dropping any single row changes the vertex set.
CAVEAT, and it is the whole point of the route: candidate exhaustiveness at
(3,7) is BYPASSED, not proved. The manifest's system_completeness field is
unchanged. Nothing here bears on (3,8) or higher; gate 2 is a separate prereg.
[docs/generator_gate1_3_7.md; docs/prereg_generator_gate1_3_7.md;
scripts/generator_gate1_3_7.py; results/data/generator_gate1_3_7.json;
tests/test_generator_gate1_3_7.py]

## Current campaigns: the (3,11)/(3,12) constraint-generation frontier

NOT superseded by the fiber program; this is where the polytope frontier
stands.

1. Stage 1, the constraint generator for d >= 11: docs/stage1_klyachko_spec.md
   holds the extracted algorithm (Theorem 3.2.1, cubicle extremal edges,
   Schubert coefficient test via lrcalc). Test-first: it must reproduce the
   five known N=3 systems before any new output counts. Until then
   constraints ship as a lookup (src/gpc_census/data/constraints.json).
2. The (3,11) bracket (docs/bracket_3_11.json) is CLOSED: all 50 outer
   candidates are settled WITHOUT Stage 1, in exact arithmetic, at
   19 TRUE + 31 REFUTED + 0 OPEN.
   scripts/settle_bracket_3_11.py regenerates docs/bracket_3_11_settlement.json
   at 19/31/0 using five exact methods. The first four settle the same 46 they
   always did: 19 TRUE (17 by state transport from the rank-9/10 census, the
   uniform (3/11)^11 by an explicit Z11 difference design, and a frozen-core
   paired state); 27 REFUTED (frozen-core pinning to N=2 plus the
   even-degeneracy pairing theorem, and zero-restriction against known
   lower-rank GPCs). The fifth, LEVEL5-RESSAYRE, is applied last and refutes
   exactly the four that were OPEN (cands 23, 26, 34, 44), each record naming
   the certified row it violates and the excess; see item 3.
   CAVEAT unchanged: the refutations do NOT finish (3,11) as a polytope;
   cutting them creates new true vertices absent from the outer list, so
   Stage 1 (or a tighter bracket) is still required. The settled points are a
   mandatory acceptance oracle for any (3,11) generator.
3. The level-5 series inequalities are PROVED, which closed the 4 OPENs
   (docs/level5_ressayre_3_11.md, scripts/level5_ressayre_3_11.py,
   tests/test_level5_ressayre_3_11.py). The four rows that kept 23, 26, 34
   and 44 open, AK-RMK-4.2.1 and KLY09-EXT-1/2/3, all of the form
   lambda_A + lambda_11 <= 2 for a second-level quadruple A, carry exact
   Ressayre certificates from the repository's own generation package. That
   is branch (a) of the old fork, and it lands: each candidate violates at
   least one certified row, so each is outside the polytope. Branch (b),
   attaining idx 44 to falsify KLY09-EXT-1, is dead; its design route was
   already closed by exact counting in
   docs/rank11_open_candidate_designs.md.
   scripts/partial_families.py tiers these four RESSAYRE at rank 11 ONLY. The
   rank-12 rows share index sets but are wedge^3 C^12, a different
   representation, and stay CLAIMED: restriction carries a valid (3,12) row
   down to (3,11), never the reverse.
   The certificate proves VALIDITY of each row, not facetness or
   irredundancy, and the screen records both as not established. Validity is
   what refutes a candidate, so that is enough here.
   The four certificates are replayed by
   scripts/verify_level5_ressayre_3_11_standalone.py, which imports nothing
   from gpc_census, derives the fermionic sign by inversion counting rather
   than position bookkeeping, and recomputes the determinant by exact
   Fraction elimination rather than Bareiss. 64 checks, and the test suite
   asserts it fails on a tampered determinant or hyperplane. This closes the
   gap the rank-7 checkpoint records, that the Ressayre determinants sit
   outside its independent verifier's trust boundary.
   The screen is calibrated at ranks 9 and 10, the largest published fixed-N=3
   systems, by a COMPLETE test: 2,898 subset rows screened, 50 certified, and
   every certified row maximized by LP over the published polytope with
   0 invalid. Negative control: 1,857 demonstrably false rows, 0 certified.
   Caveat recorded in the artifact and the doc: the trace count did nearly all
   the rejecting, only 6 false rows reached the tangent determinant, and a
   targeted hunt over perturbed published facets found no more, so the
   determinant step is the least stressed part of the calibration.
   Superseded by this, retained as history: the normalized contraction audit
   in results/data/rank11_contraction.json, which left positive numerical
   floors on all four and identified cand 44 at 11/4320 as the sharpest
   target. It pointed the right way and was never a certificate.
4. The (3,12) footholds: 19 certified true vertices by face embedding of the
   (3,11) settlement (docs/bracket_3_12_true_vertices.json), and the exact
   plethysm Stage-0 inner cloud (docs/bracket_3_12_stage0_inner.json,
   scripts/plethysm_inner_hull.py), whose cloud-extremality is NOT
   vertexhood until the cloud converges.

Other census-era leads (vertex-cone arithmetic, holonomy Galois structure,
state transport, the v96 campaign, core/hole isomorphisms, Levi symmetry
theorem and its audit) live in the archive under their section names; the
Levi identity rank_deficiency == stab - 1 (799/799) and the hard/soft
plethysm criterion feed Problem 7 below.

## Source documents in docs/

- facet_index_stability.md: the facet-index probe. Refutes the AK2008
  coefficient-one facetness criterion with 59 counterexamples, refutes
  trailing-zero padding as a map on facets, lands the normal-shape compression
  census over all twelve published systems, and prices the surviving sound
  prefilter at 1.66x and 1.79x on the reduction stage. Pre-registration in
  prereg_facet_index_stability.md, all six predictions PASS; artifact
  facet_index_stability.json; guard test tests/test_facet_index_stability.py.
- bdr_constructor_comparison.md: the external-constructor comparison. Bounds
  what any external candidate generator can remove under the local
  certification rule (1.82x, 6.27x, 16.38x at `(3,8)`), recertifies all 36
  published rows of `(3,6)`, `(3,7)` and `(3,8)`, and labels the benchmark
  INCOMPLETE because the pinned backend cannot be fetched or run.
  Pre-registration in prereg_bdr_constructor_comparison.md, P5 FAIL recorded;
  artifact bdr_constructor_comparison.json; guard test
  tests/test_bdr_constructor_comparison.py.
- quantization_character.md: the follow-up that closed the period question.
  The quantization period is the order of a finite isotropy character of the
  certified attainer, matched at 14 of 14 vertices. Pre-registration in
  prereg_quantization_character.md; artifact quantization_characters.json;
  standalone replay through wedge minors in
  scripts/verify_quantization_character_standalone.py.
- symplectic_invariants.md: the cross-realization campaign. Lattice tangent
  cones and exact `U(d)` stabilizers for all 799 rows, plus quantization
  multiplicities at the 14 vertices of (3,6) and (3,7). Scorecard P1-P5 PASS,
  P6 and P7 FAIL, both informatively. Pre-registration in
  prereg_symplectic_invariants.md; artifacts symplectic_invariants.json and
  quantization_multiplicities.json; standalone replay in
  scripts/verify_symplectic_invariants_standalone.py.
- sine_field.md: THEOREM F, and the SINE FIELD, an object the arithmetic arc
  recorded but never used. 🟦 PROVED: with A_ab = Im(z_a conj(z_b)) the signed
  area, the Gram entry and the area are the real and imaginary parts of ONE
  product, t_ab^2 + A_ab^2 = r_a^2 r_b^2, and Theorem E's discriminant is
  Delta/4 = A_ab^2 A_bc^2. So rationality turns on whether sin^2(Phi_ab) and
  sin^2(Phi_bc) share a squarefree part: a condition on Q(sin Phi), NOT on the
  holonomy field Q(cos Phi) that the whole arithmetic arc studies. The
  sine-addition identity r_b^2 A_ac = A_ab t_bc + t_ab A_bc makes that field
  propagate. ❗ CROSS-CHECK, 17 of 17: every squarefree part computed here is a
  recorded sin_radicand of holonomy_fields.json, generated by a different script
  for a different purpose, so the object is tied to certified data rather than
  invented. 🟦 The sine field is a CLASS INVARIANT, 44 of 44 non-vacuous pairs,
  tested on all 54 pairs including the 28 with irrational t so the measurement
  is not near vacuous. Values 0, 7, 15, 23, 47, 71 over 11 transport classes.
  Every nonzero area is IRRATIONAL, 0 of 44, so a nondegenerate sine field is
  always a genuine quadratic extension. 💡 THE FINDING FROM DISSECTING THE 1 OF
  18: (3,10) v75 ch(2,8) has sine field Q(sqrt 15) and NO cosine radicands at
  all, because both its loops have degree 1. The holonomy census records that
  state as arithmetically trivial; its arithmetic lives in the sine field. It is
  also the class that produced the single non-degenerate Theorem E triple, so
  the 17 degenerate triples say nothing about the field and the one
  non-degenerate triple sits exactly where the field is invisible. ❗ P-F-5
  passed only WEAKLY, 1 of 15, because the comparison is against the union of
  cosine radicands over all loops of a state rather than per loop; the sharp
  version has not been run. ❗ A measurement bug was caught before scoring:
  in_sine_field read 0 of 44 because sympy.primitive_element returns its
  polynomial in its own symbol and the degree was read against an unrelated one.
  A rational-t pair passes trivially, so 0 of 44 was impossible as mathematics.
  🔬 The class invariance is MEASURED, not proved; proving it would remove the
  squareness check from Theorem E entirely.
  [docs/prereg_sine_field.md; scripts/sine_field.py;
  results/data/sine_field.json; tests/test_sine_field.py]
- class_closure.md: THEOREM E, and EVERY LARGE CLASS CLOSED, 17 of 17. 🟦
  THEOREM E PROVED, the statement Theorem D was a special case of: the sides lie
  in the real plane so the Gram matrix has rank at most 2, and its 3x3 minor on
  any triple {a,b,c} is a RATIONAL quadratic in t_ac as soon as t_ab and t_bc
  are rational, with discriminant factoring as 4 (t_ab^2 - r_a^2 r_b^2)(t_bc^2 -
  r_b^2 r_c^2) = 4 r_a^2 r_b^4 r_c^2 sin^2(Phi_ab) sin^2(Phi_bc). So t_ac always
  lies in a quadratic extension by a RATIONAL radicand, and is rational exactly
  when that is a square, in particular whenever a linked holonomy has cos = +/-
  1. The size-2 hypothesis of Theorem D was only ever the base case. 🟦 THE LAST
  CLASS CLOSED: (3,10) v103 ch(1,6) could not be reached by Theorem D because
  its only one-hop link runs to the size-4 class (3,9). Two entries of (3,9) are
  Theorem-D-proved from size-2 links at (2,6) and (1,2); Theorem E on the triple
  (0,2,3) has discriminant EXACTLY ZERO, forcing t_03 rational; the shared loop
  transports it and Theorem C finishes. The triple carries
  uses_only_theorem_D_inputs, so nothing measured is consumed and the chain is
  grounded in Theorem A throughout. ❗ The vanishing discriminant is not luck:
  v103_fiber_count.json independently certifies v103 as real up to gauge with
  holonomy (0, pi, pi, 0, 0), so every sine vanishes. The one class the size-2
  mechanism missed is the one whose geometry is degenerate, and the degeneracy
  was already on record. ❌ P-E-6 FAILED as pre-registered, and taught more than
  the successes: degeneracy is NOT a v103 peculiarity, 17 of the 18 chained
  triples have discriminant 0, spread across (3,10) v75 and both v103 classes.
  Rationality propagates mostly because linked holonomies are pinned at cos =
  +/- 1, not because discriminants happen to be squares. The single
  non-degenerate triple, v75 ch(2,8) (0,1,2), has discriminant 225/5903156224 =
  (15/76832)^2, so both modes occur. 🔬 POST-HOC and needing an out-of-sample
  test: whether degenerate dominance is a fact about real-up-to-gauge states
  specifically. Scope 17 classes, 16 states, 11 transport classes (R2).
  [docs/prereg_class_closure.md; scripts/class_closure.py;
  results/data/class_closure.json; tests/test_class_closure.py]
- shared_loop_rationality.md: THEOREM D and THEOREM C', the induction that
  turns Theorem C's checked hypothesis into a proof and lifts its conclusion off
  class size 3. 🟦 THEOREM D PROVED: if two pairs of a class (A,B) have A-side
  determinants one hop apart, u_b = u_a - C + E, then v_b = v_a - C + E too,
  both determinant pairs lie in the class (C,E), and the two Lemma 1 loops are
  the SAME integer vector; when (C,E) has size 2 Theorem A applies there and the
  weight monomial M = (r_a r_b)^2 cancels identically, so t_ab is RATIONAL. The
  rational diagonal is Theorem A read through a shared loop, not a coincidence.
  🟦 THEOREM C' PROVED: the sides lie in the real plane so the Gram matrix has
  rank at most 2, a class of size m needs m-2 rational diagonals, and a rational
  chain makes every adjoined radicand rational, hence multiquadratic at ANY
  size. That is what Proposition B could not give, since B's radicand carries
  sin psi from an earlier loop and only builds a 2-tower, which may close
  dihedrally. 🟦 Together: 16 of 17 large classes PROVED end to end from Theorem
  A with NO checked arithmetic hypothesis, and all 16 sit exactly at the minimum
  m-2 with never a spare diagonal. Cross-check 17 of 17: the predicted t_ab is
  built from the LINKING class and matches the one computed inside the original
  class, which is where a fermionic sign error would surface. ❌ REFUTED, and
  pre-registered as an expected failure: "every class of size >= 3 contains a
  pair linked by a class of size at most 2" is FALSE, counterexample (3,10) v103
  channel (1,6), whose only one-hop link sits in the size-4 class of the same
  state. It survives only because its diagonals are rational for another reason,
  which names the honest generalisation: replace "size 2" with "all pairwise t
  already rational", making the induction a mutual recursion whose
  well-foundedness is the new open problem. ❗ A measurement bug was caught
  before scoring: the chain search skipped orderings with order[0] > order[-1]
  as reverses, but reversing changes which partial sums appear, and the bug
  would have manufactured a refutation of the programme's own hypothesis.
  [docs/prereg_shared_loop_rationality.md; scripts/shared_loop_rationality.py;
  results/data/shared_loop_rationality.json;
  tests/test_shared_loop_rationality.py]
- polygon_diagonals.md: THEOREM C, the mechanism that forbids D4 at one-hop
  class size 3, and the one to read before any further attack on (R1) at large
  classes. Lemma 2 unexpanded says a class of size m is a closed (m+1)-gon with
  side lengths fixed by the weights, so Theorem A is the law of cosines for a
  triangle and a size-3 class is a quadrilateral. 🟦 THEOREM C PROVED: if a
  size-3 class carries a RATIONAL DIAGONAL then z_a, z_c and D = z_a + z_b are
  three vectors in the real plane, their Gram matrix is singular, and the
  vanishing determinant is a quadratic in t_ac with RATIONAL coefficients and
  RATIONAL discriminant; so every within-class cosine lies in a multiquadratic
  field, (R1) holds and C4 and D4 are excluded. ❗ This is what Proposition B
  could not give: B's radicand carries sin psi from an earlier loop and builds a
  2-tower, which may close dihedrally, while C's radicand is rational outright
  because the Gram relation never leaves the class. 🟦 Machine-checked, Gram
  residual exactly 0 on 17 of 17 classes (14 by triangulation, 3 already fully
  rational). ❌ P-B1-2 FAILED and is the informative part: mode-coherent partial
  sums ARE the 2-RDM entries rho^(2)_{AC,BC}, verified exactly 25 of 25, and 6
  of them are IRRATIONAL of degree 2, so hypothesis (H2) does NOT extend from
  the 1-RDM to the 2-RDM. P-B1-4 failed the other way, 11 of 33 mode-incoherent
  subsets are rational, so mode-coherence is the wrong invariant entirely. 🔬
  POST-HOC, labelled: all 17 classes carry a rational diagonal, 14 of them
  exactly one, which is what makes Theorem C bite and which needs an
  out-of-sample test. The residue, and the new open problem, is sharper than
  what it replaces: prove every size-3 one-hop class carries a rational
  diagonal, or exhibit one that does not. Scope 17 classes, 16 states, 11
  transport classes (R2).
  [docs/prereg_polygon_diagonals.md; scripts/polygon_diagonals.py;
  results/data/polygon_diagonals.json; tests/test_polygon_diagonals.py]
- screening_orbit_structure.md: SCREENING IS NOT A WEYL-ORBIT INVARIANT, a
  refutation of my own next step, with the decomposition that survives. ❌
  REFUTED: screening one representative per orbit and transporting the verdict.
  At rank 7, 25 of 35 oriented-tau orbits carry MIXED verdicts, several holding
  both certified and trace_rejected. The reason is structural: the dominant
  chamber is a fundamental domain and so is exactly what the Weyl group does not
  preserve, and the moment cone's Weyl invariance is a statement about the
  AMBIENT cone, not about a chamber-facing inequality list. 🟦 WHAT SURVIVES,
  and it is the generator's new cost model: the trace test compares negative
  roots with h_j - h_i < 0 against weights strictly below, the second count is
  S_d INVARIANT (verified on every row at ranks 6, 7, 8, zero failures) and the
  first is EXACTLY the inversion number of h (verified against
  _negative_root_set), which spreads over a Mahonian range inside one orbit. So
  the test reads inv(sigma h) == B with B invariant: B costs one evaluation per
  ORBIT rather than per row, and each row then costs only an inversion count. 🟦
  Measured survivor fractions 8.84%, 5.14%, 1.98% at ranks 6, 7, 8, FALLING, so
  the expensive determinant stage shrinks as a share of the work. The 549
  survivors at rank 7 cross-check exactly against the methodology doc's 10,133
  trace failures out of 10,682. ❗ orbit_canonical_flats.md is CORRECTED: it said
  orbit invariance "is true for validity" and justified the reduction by
  claiming orbit members are facets together or not at all. Both wrong. The
  enumeration reduction stands instead on the S_d STABILITY of the candidate
  set, which is separately verified with images_outside_set zero at every rank.
  [scripts/screening_orbit_structure.py;
  results/data/screening_orbit_structure.json;
  tests/test_screening_orbit_structure.py]
- parabolic_survivor_traversal.md: THE PARABOLIC REFORMULATION IS RIGHT AND THE
  SPEEDUP IT MOTIVATES IS NOT. 🟦 CONFIRMED: arrangements of h are the parabolic
  quotient S_d/(S_m1 x ... x S_mr), the trace inversion number is Coxeter length
  on minimal coset representatives, and the Gaussian multinomial already used by
  trace_survivor_count is that quotient's Poincare polynomial. ❗ CONVENTION IS
  LOAD BEARING: minimal coset length equals the arrangement's inversion number
  only with h sorted INCREASINGLY; sorted decreasingly the base arrangement is
  the LONGEST representative and neither the min nor the max coset length
  matches. ❌ REFUTED, first obstruction: a fixed-length Gray code cannot exist,
  because an adjacent transposition changes Coxeter length by exactly one, so
  two distinct length-B representatives are never adjacent. Measured, 562
  adjacent swaps between same-multiset arrangements, ZERO preserve the inversion
  number. ❌ REFUTED, second obstruction: the DFS-prefix increment that DOES
  exist is exact (partitions bit identical to a full rebuild at ranks 7 and 8)
  and saves almost nothing, ratio 1.03 and 1.07 against a predicted d/3 of 2.33
  and 2.67, and it is SLOWER in wall time (0.11 s to 0.17 s at rank 8). The
  prediction assumed prefix sharing and there is almost none: 52.2 triple sums
  per survivor against 56 for a full rebuild, so the DFS tree is close to a star.
  🔬 The stage was under one percent of the run anyway, 0.11 s of 12.4 s, after
  the futile-evaluation fix took rank 8 from 66 s to 12.4 s. What survives is the
  identification itself as the right language, with its convention pinned.
  [scripts/parabolic_survivor_traversal.py;
  tests/test_parabolic_survivor_traversal.py]
- determinant_matching_audit.md: THE TANGENT DETERMINANT CLASSIFIED BY ITS
  MATCHINGS, and the profiling result that came out of it. 🟦 Each perfect
  matching of the support graph carries a sign and a monomial, so keeping the
  UNSIGNED multiplicity beside the signed coefficient separates cases the signed
  sum conflates. A monomial reached by exactly ONE matching has coefficient +-1
  and cannot cancel, which proves nonvanishing with no summation, no evaluation
  and no large-integer work. Rank 7: 474 structural_zero, 25 exact_cancellation,
  16 unique_matching, 24 unique_exposed_monomial, 10 nonzero_after_collision,
  summing to 549 and reproducing the published 499/50 split. 40 of the 50
  nonzero determinants, 80%, are certifiable without summing anything; largest
  matching count seen is 30,528, so this is not a small-case artifact. Rank 8
  reproduces the predicted split exactly: 6,414 structural, 55 exact
  cancellation, 193 Hall-feasible splitting 55 zero and 138 nonzero, and the
  share needing no summation RISES, 80.0% to 84.1%, at a largest matching count
  of 9,538,560 and order 22. 🟦 THE CERTIFICATE IS NOW POLYNOMIAL: weight each
  variable, and a matching strictly lighter than every other carries a monomial
  no other matching reaches, so its coefficient is +-1 and the determinant
  cannot vanish; the certificate is a weight vector, a matching and the
  second-best weight, replayed with one assignment solve per matched edge. It
  issues on 36 of the 50 rank-7 nonzeros against the audit's upper bound of 40,
  the gap being weight-vector choice rather than a limit of the method. Safety
  is the tested direction: no certificate was ever issued for a vanishing
  determinant, every one replays, and tampered ones are rejected.
  ❗ THE BIGGER FINDING was in the profile, not the audit: at rank 8, 48.6 of
  the orbit-native constructor's 54 seconds sat inside screening, with
  enumeration 4.8 s and candidate reconstruction 0.7 s. The cause is that
  assess_candidate_by_evaluation computes up to 32 exact Bareiss determinants
  looking for a nonzero evaluation, and on an IDENTICALLY ZERO determinant every
  one returns zero; 6,414 of 6,607 rank-8 survivors are structurally zero, so
  the budget was being spent to reach a foregone conclusion before falling
  through to the Hall test that settles it at once. Running Hall FIRST proves
  the search futile instead of performing it: (3,8) goes 251.9 s reference to
  66.0 s orbit-native to 12.4 s, 20x, partition unchanged at 137/6,469/1/0.
  The verdict including the attempt count is byte identical by design, because
  the field records that the budget yields no nonzero evaluation and a Hall
  violator establishes exactly that for every trial point at once; the rank-7
  manifest still regenerates byte for byte. 🔬 unique_exposed_monomial is
  currently DETECTED by the full expansion it is meant to avoid; the cheap
  separating-weight certifier is not implemented, so 80% is what is CERTIFIABLE
  this way, not what is certified this way. The 25 exact cancellations are
  unexplained.
  [docs/determinant_matching_audit.md;
  tests/test_determinant_matching_audit.py]
- orbit_native_constructor.md: THE FAST COMPONENTS ARE NOW ON THE ACTUAL PATH,
  which they were not. ❗ FINDING FIRST: orbit enumeration, direct survivor
  generation and the Hall filter were each proved and measured, and the
  production constructor called NONE of them. It materialised every hyperplane,
  oriented each into two taus and screened the whole list, and the screening
  fallback ran the exponential determinant DP with no matching pre-test. Every
  measured speedup sat beside the pipeline rather than in it. 🟦 HALL is now
  inside the symbolic fallback with a DELIBERATELY INDISTINGUISHABLE outcome,
  same status and same reason, so stored manifests stay byte identical and it is
  a pure accelerator: the rank-7 reference manifest regenerates byte for byte
  and rank-8 screening counts match the stored artifact. Of 6,469 rank-8 rows
  reaching the fallback, 6,414 are decided by matching and 55 reach the DP.
  🟦 ORBIT-NATIVE CONSTRUCTOR builds a complete system from orbit
  representatives alone and never constructs the hyperplane list, the oriented
  tau list, or a trace-rejected row: 326,233 of 332,840 rank-8 rows are never
  built. Coverage is CHECKED by count conservation, orbit-stabilizer plus the
  Gaussian multinomial against what the run did, so a lost orbit cannot balance
  it. ❗ SELF-PAIRING TRAP: the oriented orbit of (h,z) already contains (-h,-z)
  when z = 0 and the multiset of -h equals that of h, so processing both
  orientations double counts; there are 1, 3 and 3 such orbits at ranks 6, 7, 8,
  which is exactly why the oriented counts are 15, 35, 109 and not 16, 38, 112.
  🟦 GATED, not merely run: at rank 7 the emitted taus equal the reference
  survivor set exactly and the partition reproduces the published 1/49/499/0; at
  rank 8 retained rows, reduction certificate, regression, certified tau set and
  symbolic-zero tau set are all identical, 66 s against 252 s. 🔬 The speedup is
  3.8x, not the 719x survivor generation shows alone, because the orbit
  ENUMERATION is now the dominant stage. The row scan is gone and enumeration
  replaced it as the bottleneck.
  [src/gpc_census/generation/orbit_native.py; tests/test_orbit_native.py]
- orbit_canonical_flats.md: WEYL-ORBIT CANONICAL ENUMERATION, and the first
  exhaustive hyperplane count at rank 9. ❗ CORRECTION FIRST: (3,7) candidate
  exhaustiveness is CLOSED by the closed-flat enumerator (5,341 hyperplanes for
  1,623,160 candidate bases, 10,682 rows all resolved, exact reduction giving
  the four published facets). AGENT_PLAN T2 step 1 said "bypassed, not proved",
  which describes only the gate-1 inner/outer closure and made rank 7 look open;
  the plan is corrected and the next gate is (3,8). 🟦 Orbit reduction is the ONE
  pruning the standing rule permits, because the moment polytope is Weyl
  invariant and the ordered slice is a fundamental domain, so nothing is
  discarded. Measured: 362 hyperplanes in 8 orbits at rank 6, 5,341 in 19 at
  rank 7, 166,420 in 56 at rank 8, reductions 45x, 281x, 2,972x, action closed
  at all three. 🟦 THE AUGMENTATION ENUMERATOR holds only representatives and
  never materialises the lattice; expanding orbits by d!/|Aut| rebuilds the
  materialised counts EXACTLY at EVERY level, not just the top. Rank 8 in 2.9 s
  against 125 s materialised. 🟦 RANK 9 REACHED, which the materialising
  enumerator could not finish: it died at level 6 having checkpointed levels 1
  to 5, and the orbit expansion reproduces those five stored counts (84, 3486,
  78274, 968121, 6383706) then continues, giving (3,9) = 10,004,154 ADMISSIBLE
  HYPERPLANES through 231 orbits, peak 494 at any level, in 407 s. 🟦 Rank-8
  level vector now has THREE independent derivations. ❗ CONVENTION TRAP, pinned
  by a test: hyperplanes are stored with first nonzero entry positive, so a
  permuted image must be re-normalized; without it the action sends 201 of 362
  images outside the set and the orbit count reads 15 instead of 8, a valid
  orbit count of a valid action on the wrong set. 🔬 This enumerates
  HYPERPLANES, not facets: each representative still needs Ressayre screening
  and exact reduction, and no (3,9) system is produced here. Orbit ratios 2.4,
  2.9, 4.1 over four points are a projection, not an asymptotic.
  🟦 COST, measured and then fixed, twice. Canonicalization: the worst per-flat
  cost was exactly d! at EVERY rank, not a large-cell effect, because plain
  colour refinement splits nothing on the most symmetric flats. Individualization
  (the nauty step, with an isomorphism-invariant target cell so the canonical
  form stays well defined) drops the (3,8) peak worst case 40,320 to 1,440 and
  the (3,8) enumeration 401 s to 16.7 s, reaching the floor leaves >= |Aut| on
  300 of 313 flats. That floor is then BROKEN where it was the real obstacle: a
  mode in no present weight is swapped with any other such mode by a
  transposition fixing the flat, so refinement can never split them and |Aut|
  carries a k! that is pure waste. At rank 1 of (3,11) that is 8 of 11 modes,
  |Aut| = 3!*8! = 241,920, and the run emitted nothing in 43 s of CPU; handling
  them analytically gives 6 leaves for the same flat with orbit size 165 =
  C(11,3), checked against a count known independently. Memory: children were
  materialised
  before being collapsed, so peak scaled with the CHILD count, which is the
  shape that killed the (3,11) run at 6.2 GB; streaming them makes peak the
  ORBIT count of the next level. (3,9) reruns to the same ladder and the same
  10,004,154 in 308 s at 112 MB. A test pins the leaves-to-|Aut| RATIO rather
  than a timing, so the factorial tail cannot return unnoticed.
  [docs/prereg_orbit_canonical_flats.md; scripts/orbit_canonical_flats.py;
  src/gpc_census/generation/orbit_canonical.py;
  results/data/orbit_canonical_flats.json; tests/test_orbit_canonical.py]
- direct_survivor_generation.md: THE ROW SCAN IS GONE, which is the limit
  screening_orbit_structure.md named for itself ("every oriented row still has
  to be touched, which is what sets the reachable rank"). 🟦 The trace test is
  inv(sigma h) == B with B an orbit invariant, so an orbit's survivors are the
  arrangements of a multiset with a PRESCRIBED inversion number; those are
  generated directly rather than filtered out. Output linearity is the claim and
  it is exact, not heuristic: pruning to residual budget in [0, maxinv(rest)]
  with maxinv = C(m,2) - sum C(m_i,2) is exactly the attainable range, because
  the Gaussian multinomial has strictly positive coefficients across its whole
  degree range, so every surviving branch completes and no node is wasted.
  🟦 Verified SET against set, not merely by count, on all 35 rank-7 orbits: a
  generator emitting the right NUMBER of wrong arrangements would silently
  replace the candidate list. Measured 549 of 10,682 rows at rank 7 (133x) and
  6,607 of 332,840 at rank 8 (719x), both reproducing
  screening_orbit_structure.json exactly, with 45 of 109 rank-8 orbits killed
  outright before any arrangement is built and 326,233 rows never visited. The
  speedup RISES with rank because it is the reciprocal of the falling survivor
  fraction, so this is a change of shape rather than a constant. 🔬 It does not
  reduce the ORBIT count, does not touch the determinant stage, and does not
  make screening an S_d invariant, which stays refuted.
  [docs/prereg_direct_survivor_generation.md;
  src/gpc_census/generation/orbit_canonical.py; tests/test_orbit_canonical.py]
- census_inner_sandwich.md: COMPLETENESS WITHOUT CANDIDATE EXHAUSTIVENESS, and
  the one to read before spending any more effort on exhaustive enumeration.
  The claim ladder allows "an exhaustive-enumeration proof OR an independent
  exact inner/outer closure", and the second route does not care where the
  candidate rows came from. 🟦 P = O PROVED at ALL NINE census systems, (3,6)
  through (3,10) and (4,8), (4,9), (4,10), (5,10), with exhaustiveness unused;
  ranks 9 and 10 had no in-repo completeness proof before, and no N > 3 system
  did. Outer half: every stored row screens as an exact Ressayre certificate
  (4/4, 31/31, 52/52, 93/93, 15/15, 60/60, 125/125, 161/161), because
  certifying VALIDITY of a supplied row costs 0.14 ms and only EXHAUSTIVENESS
  was ever expensive. Inner half: every vertex of O is the spectrum of a
  certified census state, verified by an exact characteristic-polynomial
  identity. 🟦 The screening pipeline is NO LONGER FIXED-N=3: everything under
  it (exterior weights, roots, tangent determinant, certify_candidate, the
  full-dimension witness) was already written for wedge^n C^d, so the
  restriction was seven literal 3s in candidate_pipeline.py plus one
  combinations(range(d), 3) in full_dimension.py; screen_tau_candidates now
  takes particle_number and the n=3 artifacts are byte-identical. ❗ The
  general-N screening is guarded against being vacuous by two tests: perturbing
  each (4,8) row rejects 15/15 and certifies 0, and screening the (4,8) rows at
  particle_number=3 does not certify them all. 💡 The inner_half_only label is
  what made the fix obvious: the first version of the script called the N > 3
  systems closed when only their inner half was certified, the label recorded
  the gap honestly, and the gap turned out to be plumbing. The label is kept and
  still tested, now unused. ❗ NON-CIRCULARITY is the
  whole game: the census VERTEX lists came from the constraint tables and would
  be circular, the attaining STATES did not and reference no table, which is
  why the script re-runs the state certificates rather than trusting the vertex
  list. ❌ The predicted wall was wrong: exact vertex enumeration of the
  degenerate polytope was expected to block this, and gpc_census.exact_polytope
  does every census system reproducing the certified counts, so lrs-style
  pivoting, symmetry-aware orbit enumeration and containment checking are all
  still unspent. 🟦 It did bite at the wide end, and the fix was inside the
  existing enumerator: exact_polytope.cut recomputed a vertex's tight set once
  per PAIR and an exact rank once per PAIR, so caching the tight set per vertex,
  rejecting pairs whose common tight rows cannot reach rank n-1 on an integer
  count, and memoising the rank on the common tight set give 12->1.5 s at
  (3,9), 247->11 s at (3,10), 211->22 s at (4,10) and, at (5,10), over 40
  minutes abandoned -> 166 s. All three are identities, not heuristics, and a
  test pins that by re-running the whole (3,8) enumeration against the
  unmemoised _adjacent and comparing the cut sets at every step. 🔬 BOUNDED NULL: the self-improving loop (a resisting
  vertex points at a missing facet, contraction attack + SOS + level-5 Ressayre
  turn it into a row) is implemented but never fired, zero vertices resisted at
  all nine systems, so it is exercised only on its trivial branch. Its first
  real test is (3,11), where the inner side is 19 bracket vertices rather than
  a census. 💡 Orbit reduction under S_d is recorded as available and NOT
  implemented; it is the one provably safe pruning, since the moment polytope
  is Weyl-invariant and the ordered slice is a fundamental domain.
  [scripts/census_inner_sandwich.py;
  results/data/census_inner_sandwich.json;
  tests/test_census_inner_sandwich.py]
- stage1_reference_series.md: the RANK GENERATOR's validation ladder. 🟦 GATE
  G2 CLOSED: the generator reproduces the published fixed-N=3 systems at rank
  7 (4 rows) and rank 8 (31 rows) BLIND, matching as oriented sets with
  generated_only and published_only both empty. Blindness is structural, not
  claimed: generation enumerates admissible affine flats, screens by the exact
  Ressayre criterion and reduces on the ordered Pauli slice, and the stored
  table is read only afterwards in _known_system_regression; the rank-7 gate
  is re-run LIVE in the test rather than re-read from the artifact. 🟦 The
  screening partition is exhaustive at both ranks, zero unresolved rows. 🟦
  MEMORY FIX, 4.4x with byte-identical output: a stored flat carried a
  materialised RREF (about 870 of its 950 bytes) that is a function of its own
  basis and is needed only while extending, so the level dict now stores
  mask -> basis_indices and reduces once per parent instead of once per child.
  Rank 8 peak RSS 1302 MB -> 298 MB; flat counts by rank pinned by test.
  🔬 G3 (ranks 9, 10) is NOT closed. The wall is now TIME, not memory: rank 7
  to rank 8 costs a factor of about 53 in wall time and 31 in hyperplanes, and
  the only fix that attacks it is pruning inside the enumeration (about 95
  percent of taus are trace-rejected downstream by a cheap test). Per the
  methodology labels, resource exhaustion is OPERATIONAL and is recorded as
  status=operational_failure with gate_passed=false, never as a rejection.
  ❗ Passing G2 says nothing about candidate EXHAUSTIVENESS at a new rank;
  a new-rank output stays relative to the supplied certified candidates.
  [scripts/stage1_reference_series.py;
  results/data/stage1_reference_series.json;
  tests/test_stage1_reference_series.py; gpc_census.generation]
- holonomy_two_elementary.md: the PROOF ATTACK on 2-elementarity, and the
  current state of the classes-3-and-4 gap. ✅ REDUCTION THEOREM: the
  multiquadratic conjecture follows from two conditions on one quantity per
  loop, (R1) deg_Q(cos^2 Phi) <= 2 and (R2) N(cos^2 Phi) a rational square at
  degree 4, the second being exactly the classical V4 criterion for an even
  quartic (D4 and C4 are the alternatives it excludes). ✅ THEOREM A',
  strengthening the shipped Theorem A: at one-hop class size <= 2 the cosine
  is n/sqrt(M) with n, M rational, so cos^2 is RATIONAL and the minimal
  polynomial is EVEN, which is the part that generalises. 🟦 (R1) holds 82/82
  and (R2) 13/13, and the reduction DERIVES every certified Galois group
  (C1 37, C2 32, V4 13) with no group computation. 🟦 Theorem A' has 0
  violations, and every cos^2-irrational loop sits in a state with a class of
  size 3. ❌ "The minimal polynomial is always even" is FALSE: two (5,10)
  loops have deg cos = deg cos^2 = 2; 2-elementarity survives trivially there,
  but any proof of evenness must explain them. 🔬 (R1) at class sizes 3 and 4
  is the remaining open step, and Proposition B does not supply it (squaring
  leaves a cross term). 💡 The conjugation-involution/Galois-commutant
  mechanism is stated precisely: time reversal and the fiber sign flips are
  commuting INVOLUTIONS, so the group they generate is elementary abelian and
  has no element of order 4, which is what D4 and C4 both need; the
  obstruction is that Galois permutes COMPLEX solutions while those symmetries
  act on the REAL ones, the same gap Rigidity-Rationality's descent step hits.
  [scripts/holonomy_two_elementary.py;
  results/data/holonomy_two_elementary.json;
  tests/test_holonomy_two_elementary.py]
- conditional_2rdm_program.md: steps 4-7 built on the order-three theorem.
  🟦 HIDDEN PAIR CURRENT: states with exactly the same complete 1-RDM carry
  different two-body pair currents, over the exact flux-INDEPENDENT range
  [-sqrt(15)/10, sqrt(15)/10]; the energy/current joint body is convex and is
  recovered by sweeping the support function. 🟦 The unit-edge discriminant
  gives interior branch-exchange fluxes at Phi ~ 0.0593 and ~ 3.0823, the
  genuine non-analyticities of the constrained-search functional. 🟦 The
  Slater-Condon carrier projection reproduces the full CI submatrix with max
  error 0.0, and its fermionic signs are recomputed from occupation-string
  algebra rather than assumed. 🟦 QUASIPINNING: eps_H = 2||H||(eps+eps^2) +
  ||w-w0||_1 ||y||_inf, where the second Lipschitz constant is the dual
  certificate itself; no quasipinning constant is assumed. Pinning plus the
  1-RDM narrows the achievable range by 1-2 orders of magnitude. 🟦 ATLAS: all
  799 supports tiered, 75 exactly spectrahedral (m<=3) and 724 outer-bounded,
  and only 72 have a complete coherence graph, so MOST census coherences are
  invisible to any two-body observable. ❌ Two model choices failed
  instructively and are recorded: spinful Hubbard (degenerate ground state, so
  the 1-RDM is undefined) and the symmetric ring (degenerate occupations, so
  the carrier is ambiguous). 🟦 The physical complex-flux gate is now met by
  the asymmetric spinless t-V ring with exact boundary phase `(4+3i)/5`:
  10 of 12 scan points pass, five have genuine carrier flux, and the selected
  V=2 point has ground gap 1.435, minimum occupation gap 0.00167,
  off-carrier weight 0.00466, and Im(chi)=-0.00578. An independent verifier
  rebuilds the occupation-space Hamiltonian and exterior basis change and
  passes 15/15 checks. The complex run also corrected two real-only blind
  spots: passive CI coefficients use `U.T`, not `U.conj().T`, and the
  quasipinning error uses the centered norm of the full Hamiltonian, not PHP.
  [gpc_census.pinned_physics; scripts/conditional_2rdm_program.py;
  scripts/verify_complex_flux_physical_benchmark_standalone.py;
  results/data/conditional_2rdm_program.json;
  tests/test_conditional_2rdm_program.py]
- conditional_2rdm_body.md: the REFRAMING of the observable-range work, and
  the one to read first. The endpoints of a two-body observable's range are the
  SUPPORT FUNCTION of K(gamma) = conv{Gamma^(2)_Psi : Psi -> gamma}, so the
  solver is an exact support oracle for the conditional 2-RDM body, and the
  headline is an exact conditional N-representability statement rather than two
  solved triangle Hamiltonians. 🟦 ORDER-THREE THEOREM: the variable part of
  K(gamma) is an affine image of the scaled complex elliptope E_w, whose
  extreme points obey r^2 <= m (Li and Tam), so at m = 3 every extreme point is
  a phase matrix and the support function is an EXACT SDP, not a relaxation.
  🟦 The sextic of complex_fixed_1rdm_observable_ranges is exactly the KKT
  elimination ideal of that SDP (verified, up to primitive content), which also
  explains structurally why only |h_ij|^2 and Re(chi) enter. 🟦 Exact dual
  certificates for four carriers: optimality is witnessed POINTWISE by PSD-ness
  of Diag(y)-H from principal-minor signs, so no exhaustive-root argument is
  needed, and the zero/pi-flux cases the sextic engine excludes as degenerate
  are covered uniformly. ❌ The theorem does NOT extend past order three: an
  explicit rank-two extreme point at m = 4 gives a relaxation gap of 6.9e-2, so
  m >= 4 supports ship as outer bounds until proved otherwise. Steps 4-7 are
  now implemented in conditional_2rdm_program.md; the principal open problem
  here is structured exactness beyond carrier order three.
  [gpc_census.elliptope; scripts/conditional_2rdm_body.py;
  results/data/conditional_2rdm_body.json; tests/test_conditional_2rdm_body.py]
- complex_fixed_1rdm_observable_ranges.md: exact constrained observable ranges
  on the pinned Borland-Dennis fiber for rational Gaussian Hermitian carriers.
  🟦 The gauge-invariant cycle flux chi = h01*h12*conj(h02) is the new
  invariant: edge phases are gauge-dependent but chi is not, so two carriers
  with equal diagonals and edge magnitudes can have different ranges. Phase
  elimination gives a quartic branch curve whose critical values are a
  primitive sextic in Z[z]; endpoints are its extreme real roots by exact
  Sturm isolation, with a subresultant linear lift proving every real root has
  a real (not spuriously complex) critical preimage. ✅ Verified independently
  here: brute-force optimisation of the original two-phase objective
  reproduces both ranges to 1e-15, the fiber's 1-RDM is the stated target and
  is diagonal (the support is one-hop free), A*D-C^2 = delta^2 holds, and the
  three fermionic channel signs of the at-most-two-body realisation check out
  by hand. 🟦 The unit-edge flux family contains BOTH earlier real-symmetric
  regressions as exact u = +-1 boundary points, agreeing with
  fixed_1rdm_observable_ranges.json, which a different engine produced; pinned
  by test_flux_family_reproduces_the_real_symmetric_regressions. Scope limits
  worth respecting: rational Gaussian entries only, generic (nonzero flux,
  squarefree sextic, coprime lift) only, and one fixed three-determinant
  fiber. Complex flux is NOT generic for ordinary real-orbital molecular
  Hamiltonians; it wants magnetic fields, twisted boundaries, complex Bloch
  orbitals, or spin-orbit structure.
  [scripts/complex_fixed_1rdm_observable_ranges.py;
  scripts/verify_complex_fixed_1rdm_observable_ranges_standalone.py;
  results/data/complex_fixed_1rdm_observable_ranges.json;
  tests/test_complex_fixed_1rdm_observable_ranges.py]
- lambda3_c7_polytopes.md: the ORBIT-CLOSURE layer (entanglement polytopes),
  one polytope per SLOCC class rather than one per (N,d). 🟦 CERTIFIED EXACT
  for 11 of 13 class polytopes: all four of Lambda^3 C^6 (reproducing Walter,
  Doran, Gross, Christandl) and seven of the nine of Lambda^3 C^7, including the
  decisive gate that the OPEN ORBIT reproduces the ambient (3,7) polytope
  exactly. Classes VI and VII ship as certified inner/outer brackets. Class
  VIII is closed by an exact six-term interference state with spectrum
  `(5/9)^3,(1/3)^4` and symbolic orbit tangent rank 34. Both
  bounds are exact: the inner from torus orbits of one-hop-free supports (a
  COMPLETE enumeration, because one-hop-free supports have vanishing
  coefficient moduli) plus certified degenerations; the outer from a
  Schur-Horn/Newton dominant-weight family that generalises Ky Fan. 🔬 The
  outer family is SATURATED (sweeping the weight bound 3 to 7 changes
  nothing), so the last two classes need Ressayre inequalities FOR THE ORBIT
  CLOSURE, or covariant nonvanishing. Do NOT read stage1_ressayre_3_7 as
  closing them: that manifest derives the AMBIENT cone, which is the open
  orbit's polytope (class IX), already EXACT here. The agreement is a
  cross-check on the ambient input and is pinned by a test; the bracketed
  classes are untouched by it. 🔬 Class VII has a numerically attained missing
  inner vertex and is the next exactification target. Whether class V (orbit
  dim 26) lies in the closure of class VI (orbit dim 28) is open and recorded,
  not guessed.
  [gpc_census.orbit; scripts/lambda3_c7_polytopes.py;
  scripts/lambda3_c7_class_viii_witness.py;
  scripts/verify_lambda3_c7_class_viii_witness_standalone.py;
  results/data/lambda3_c7_polytopes.json;
  results/data/lambda3_c7_class_viii_witness.json;
  tests/test_lambda3_c7_polytopes.py]
- HANDOFF.md: the note to the next agent from the 2026-07 gauge session, which
  supersedes two instructions of the handoff that preceded it (the v103
  endpoint is not diagonal, and 8 of the 71 cascade endpoints are new states).
  Read it before acting on any support claim.
- RESEARCH_ARCHIVE.md: the frozen census-era evidence log (this file's
  status tags point into it).
- THEOREMS.md: numbered theorem statements (T1 Jacobian program, T4
  endpoint theorem, T7 stratification reduction, v_B anomaly resolution).
- holonomy_rationality.md: the arithmetic of loop holonomies. Theorem A
  (PROVED) settles the Holonomy Rationality conjecture for one-hop channel
  classes of size at most 2, with the adjoined radicand a weight monomial;
  Proposition B (PROVED) is the reduction step at larger classes and yields
  only a 2-tower, which is weaker than multiquadratic, so the conjecture
  stays open there on 82/82 evidence. Also carries the trisection
  dissolution and the analysis of the one non-rational-weight ledger record.
- prereg_denominator_arithmetic.md: the pre-registered denominator-arithmetic
  predictions, committed before the census run, scored FAIL/FAIL/PASS in
  results/data/denominator_arithmetic.json with mechanisms for both FAILs.
- silence_rank_lemma.md: the silence-rank lemma, PROVED (first-order
  statements, degenerate perturbation theory), with verified census
  corollaries: v89 fixed-spectrum rigid, v103 fixed-rho rigid but
  fixed-spectrum deformable.
- unified_fiber_dimension.md: the unified first-order fiber-dimension
  formula (both regimes), its corpus verification over all 156 interference
  states, the v103 sparsification cascade, and merge notes. Its s_Q(v103)
  <= 12 is superseded by s_Q^free(v103) <= 10, certified exactly; see
  cascade_census.md.
- gauge_minimization.md: minimization of every certified state's support over
  its natural-orbital gauge family, the exact profile-class/exterior-contraction
  lower bound (670 of 799 certified gauge-minimal), the null result on the other
  129, the power calibration of the search, and the scored pre-registration.
- prereg_gauge_minimization.md: that pre-registration, committed before the run.
- cascade_census.md: the census-wide sparsification cascade, the exact
  certification of all 71 improved support bounds, the v103
  support-10 state, the cancellation-regime census, and the
  representative-swap check against the v0.12 real v_B record.
- prereg_bracket_3_11_3_12.md: the scored out-of-sample pre-registration on
  the (3,11) and (3,12) bracket holdout (4 PASS, 2 FAIL), including the
  power caveat stated before scoring, AND the 2026-07 PROVENANCE CORRECTION:
  34 of the 38 rows are transports or face embeddings of 17 rank-10 donors,
  so the effective independent sample is 2 states, not 38 rows. Read the
  corrected readings, not the original verdicts.
- prereg_template.md: the standing rules for every pre-registration in this
  project. R2 is the one added by that correction: a measurement derived by
  transport, padding, or face embedding is NOT an independent out-of-sample
  observation, and every scorecard must report independent sample size
  separately from row count.
- two_block_family.md: the four two-block interference vertices in full, and
  the census-wide result that a CROSS-BLOCK one-hop channel is NECESSARY for
  the two fibers to disagree (799/799, zero exceptions). Also kills the
  proposed cross-block counting rules and shows two-block is not the
  operative invariant.
- fiber_symmetries.md: the pair-incidence separation proposition, exact
  nonzero-minor certificates for all 799 shipped supports, the time-reversal
  involution on the exact v_B slice, and the cubic 2-RDM moment that recovers
  its family parameter modulo every one-body unitary.
- prereg_two_block_inversion.md: that study's pre-registration, committed
  before the cross-block counts were computed.
- prior_art_roadmap.md: the outside-GPC prior-art map (phase retrieval,
  rigidity theory, structured inverse eigenvalue problems, numerical
  algebraic geometry, tensor rank, ED degrees, SOS) with an ordered action
  list mapped to the open problems below.
- stage1_klyachko_spec.md: the constraint-generation algorithm, verbatim
  from the Altunbulak thesis, with the implementation plan.
- ak2008_extracted.md: inequality and vertex tables for ranks 6-8 from the
  2008 paper, including extremal states (surd coefficients need re-checking
  against the PDF before use as ground truth).
- bracket_3_11.json / bracket_3_11_settlement.json: the (3,11) Stage-0
  bracket and its 46/50 settlement (19 true, 27 refuted, 4 open);
  reproducible via scripts/settle_bracket_3_11.py.
- prereg_rank11_candidate44.md: the frozen theorem/falsification fork for the
  normalized candidate-44 floor, exact projected spectrum, and required proof
  objects. The numerical audit is results/data/rank11_contraction.json.
- partial_families_3_11_3_12.json: the published stable Grassmann families
  for N=3 at ranks 11 and 12, tagged PROVED / WEAK / CLAIMED. NECESSARY
  conditions, an outer bound only; use as a Farkas-style acceptance oracle
  for any Stage-1 output (scripts/partial_families.py).
- bracket_3_12_true_vertices.json: 19 certified true vertices of (3,12) by
  face embedding of the (3,11) settlement.
- bracket_3_12_stage0_inner.json: the (3,12) exact-plethysm inner cloud
  (scripts/plethysm_inner_hull.py); cloud-extremality is not vertexhood
  until convergence.

## Facet-index stability: the coefficient-one criterion is refuted (2026-08-08)

Follow-up to the BDR comparison, which identified redundancy elimination as the
stage with the worst long-run scaling. The question was whether an intrinsic
criterion can replace it.

❌ DISPROVED, 59 exact counterexamples. The AK2008 coefficient-one criterion is
NOT a facetness test. Exact cyclic-Schubert evaluation is deferred past
reduction in the production pipeline, so the coefficients of the rows reduction
REMOVES had never been computed at any rank. Computing them for every certified
row at `(3,7)` and `(3,8)`, 49 and 137 rows: every facet has coefficient one
(4 of 4, 31 of 31) and so do 23 of 45 and 36 of 106 removed rows. Coefficient
one is necessary and not sufficient. The standing refusal in
`docs/fixed_n3_generator_methodology.md` to use it as a facetness or
irredundancy gate now has counterexamples behind it, not just caution. Scope:
this is one selected cyclic coefficient, NOT the general Belkale-Kumar or
Ressayre deformed product, which remains untested here.
[docs/facet_index_stability.md; docs/prereg_facet_index_stability.md;
results/data/facet_index_stability.json;
tests/test_facet_index_stability.py]

🟦 CONFIRMED at two ranks: the forward direction survives, which makes
coefficient one a SOUND PREFILTER rather than a criterion. It cuts the
reduction input from 49 to 27 and from 137 to 67, gives 1.66x and 1.79x on the
reduction stage, and returns an identical retained set. NOT gate-eligible: its
soundness rests on an observation at two ranks. THE NEXT THEOREM, named: every
facet of the ordered wedge^N C^d moment polytope has cyclic-Schubert
coefficient one. That forward-validity statement is the whole remaining gap
between a measured 1.79x on the growing stage and a shippable gate.

🟦 CONFIRMED at merge, 8 of 8 maps, three particle numbers, reconciling this
with the ordered-stability campaign: a trailing-zero padded facet lands in the
higher facet list EXACTLY when the trailing-zero extension coincides with that
campaign's TIGHT extension `pi_pad` (counts 0, 4, 10, 37, 0, 8, 19, 0 on both
sides). The two results read as contradictory only until one notices the maps
differ. The trivial extension succeeds precisely where it was already tight, so
what is dead is the trivial extension, not extension of normals as such. An
earlier draft here overreached by concluding a stability theorem must act on
something other than the normal vector.

❌ DISPROVED: trailing-zero padding is not validity preserving on facets. Padded
lower-rank rows are VIOLATED at the next rank in 7 of the 8 published maps where
validity is decidable, including Borland-Dennis padded to rank 7. Expected once
stated properly, since the rank-`d` polytope is a FACE of rank `d+1` and a facet
of a face need not extend. Sharper structural fact, holding in every map across
all three particle numbers: a padded facet is never merely redundant, the
strictly-interior count is ZERO everywhere. So any stability theorem must act on
something other than the normal vector.

🟦 CONFIRMED across all twelve published systems: the normals compress. Facet
counts run 1 to 161 while distinct sorted-`tau` shapes run 1 to 21, a ratio at
or below 0.29 everywhere and at or below 0.22 wherever there are more than four
facets. The explosion is in placement, not in shape. The `N=3` entry bound of 7
at ranks 8 to 10 does NOT transfer across particle number (15 at `(4,10)`, 34 at
`(5,10)`), exactly as `tau_i = N a_i - b` implies, so no entry bound uniform in
`N` is claimed. Rank-9 and rank-10 rows are conditional on AK2008 completeness.

Cheapest next experiment: the same census at `(3,9)`, `(4,8)` and `(4,9)`, which
tests the forward direction at a third rank and a second particle number at
once, and costs one run of an already-written script.

## External constructor comparison: Bulois-Denis-Ressayre (2026-08-08)

The question was whether the Bulois-Denis-Ressayre fixed-representation
moment-cone algorithm (arXiv:2505.08812; code `ea-icj/moment_cone` pinned at
`7ab347d9a7837d68d92bde9fac74606913f395ca`) can remove exact work from the
constructor. It cannot be answered by running it here: the code host, the
landing page, arXiv and HAL are all refused by the session egress policy, the
package is not on PyPI, and it needs SageMath. The benchmark is labelled
INCOMPLETE and no external runtime or candidate count is recorded anywhere.

🟦 CONFIRMED, and this is what settles the practical question anyway. Under the
standing rule that every published row carries a locally replayed certificate,
the speedup available to ANY external candidate generator was bounded by direct
measurement at `(3,8)`: 1.82x for a free and perfect birationality prefilter,
6.27x if its candidate generation were also free, against a circular
oracle-handed-the-answer figure of 16.38x. The 6.27x is the non-circular bound,
because no generator can know which certified rows survive reduction without
running the reduction. RECOMMENDATION: do not integrate, do not spend another
sprint at ranks 7 and 8.
[docs/bdr_constructor_comparison.md; docs/prereg_bdr_constructor_comparison.md;
results/data/bdr_constructor_comparison.json;
tests/test_bdr_constructor_comparison.py]

🟦 CONFIRMED, 36 rows across 3 systems: every published nonstructural row of
`(3,6)`, `(3,7)` and `(3,8)` converts to primitive `tau`, screens `certified`,
replays and binds its certificate, and matches a local orbit. External output
is usable as candidates here, and no external verdict is taken on citation.

🔬 OPEN, and it is a correction to a shipped claim.
`docs/orbit_native_constructor.md` states that orbit enumeration "is now the
dominant stage". At `(3,8)` on this machine it is not: enumeration 6.25 s,
screening 7.64 s, composition 2.11 s, so enumeration is 39.1 percent across
three stable runs. That document measured 66.0 s for the same construction
against 16.00 s here, so this is a different machine, not a refutation, and
multi-machine results are never averaged. Re-measure on the original machine
before quoting the claim again.

🔬 OPEN trend, two points, explicitly NOT a statement about all `d`: the oracle
ceiling FALLS from 21.41x at rank 7 to 16.38x at rank 8 while the row-reduction
ratio rises from 137 to 213. The mechanism is that the oracle path is itself
dominated by composition, which scales with the retained facet count (1, 4, 31,
52, 93 at ranks 6 to 10). Smallest next experiment: the same three-path split at
`(3,9)`, one run of an already-written script.

Correspondence, from the in-repo call-graph audit in
`docs/fixed_n3_generator_methodology.md` as a SECONDARY source, not a fresh
source read: the two algorithms genuinely agree on the trace stage, where BDR's
`List_Inv_Ws_Mod` and the local direct survivor generator are both inversion-set
enumerations, and that is already the local pipeline's fastest stage. They
differ on candidate generation (1-PS versus closed flats), admissibility
(`is_sub_module` versus affine rank), and determinant handling (their
birationality filter versus exact fail-closed certification). Neither algorithm
claims anything uniform in `d`.

## Adjacent literature to respect

Liebert, Lemke, Altunbulak, Maciazek, Ochsenfeld, Schilling,
PRResearch 7, 023247 (2025), arXiv:2502.15464: spin-adapted GPCs for
larger systems. Different lattice of problems (spin sectors) from the
spinless census here; draw the boundary explicitly in any new paper.
Maciazek and Tsanov (J. Phys. A 50, 465304): doubly-excited inner bounds
coincide with the polytope for (3, d <= 7); the census confirms this
empirically (those systems are interference-free).
For the fiber program's outside-field siblings, see
docs/prior_art_roadmap.md.

# Part II. The inverse-fiber program (primary research product)

## Mission

Develop a mathematical theory of the **inverse geometry of the fermionic
one-body spectral (moment) map**, centered on fibers and sparse realization
rather than isolated representative states.

## Central objects

For a determinant support S:

- Incidence equations: A_S p = lambda
- Coherence equations: C_S(p, theta) = 0

Support-restricted fiber:

F_{S,lambda} = {(p,theta): A_S p = lambda, C_S(p,theta)=0}/gauge

**Gauge convention (mandatory precision).** "/gauge" means the quotient by
U(1)^d single-particle phases and the global phase ONLY. The further
quotient by the full degeneracy stabilizer U(m_1) x ... x U(m_k) of the
target spectrum defines the REDUCED fiber F_red, a different object. Every
fiber-dimension statement in this program must name which of the two it is
about, and must name FIXED-SPECTRUM vs FIXED-RHO. This is not pedantry: the
fixed-rho fiber vs fixed-spectrum fiber conflation was an actual bug
(levi.py audit, RESEARCH_ARCHIVE "BUGS FIXED" item 4), recurred in the P4
scoring (annotated in the archive), and the T1 upper bound fails precisely
because vertices are singular values of the moment map (THEOREMS.md T1).

Working philosophy:

- Image geometry (GPC polytope) describes which spectra are possible.
- Fiber geometry describes ambiguity, rigidity, coherence, deformation, and
  sparse realization complexity.

## Current state of the theory

### ✅ Proved

- The silence-rank lemma (first-order; docs/silence_rank_lemma.md, proof by
  degenerate perturbation theory). At a diagonal-1-RDM state over a vertex:
  (i) first-order fixed-SPECTRUM weight deformations are
  ker A_S ∩ ker J_within (within-block silent channels only; cross-block
  silence imposes no first-order fixed-spectrum condition);
  (ii) first-order fixed-RHO weight deformations are ker A_S ∩ ker J_all;
  (iii) the non-gauge fixed-spectrum tangent dimension is
  (dim K - rank J_within^Re) + (|S| - g - rank J_within^Im).
  Verified corollaries (scripts/reproduce_next_claims.py): v89 first-order
  fixed-spectrum rigid (formula gives 0, direct tangent kernel = gauge);
  v103 fixed-rho rigid (all-channel silence rank 5 = kernel dim 5) but
  fixed-spectrum DEFORMABLE (formula and direct tangent both give 6).

- Holonomy quadraticity at small channel classes (Theorem A,
  docs/holonomy_rationality.md). If a state has rational squared weights and
  rational channel targets, every one-hop class of size 2 gives a loop
  holonomy whose cosine lies in Q(sqrt(N)) for N a rational MONOMIAL in the
  weights. Proof by expanding the channel modulus; the v_B phase
  cos(gamma) = 3/(4*sqrt(14)) is exactly this formula.
  [docs/holonomy_rationality.md; results/data/holonomy_fields.json]

### 🟨 Theorem candidates

- Compact semialgebraic support layers Delta^(k). (Route: Delta^(k) is the
  image of a compact semialgebraic set under a polynomial map;
  Tarski-Seidenberg + compactness of continuous images. Near-free; write
  it.)
- Finite inverse stratification. (Route: already reduced in THEOREMS.md T7
  item 1 to Hardt trivialization + Prop 1 (semialgebraicity of the moment
  map, a polynomial map on the unit sphere). The written-proof gap is Prop 1
  plus a hypothesis check, not new mathematics. Cross-reference T7; do not
  re-derive.)
- Fixed-support realization is decidable by real quantifier elimination.
  (Route: the fixed-support fiber equations are polynomial; Tarski.
  Near-free.)
- Local constancy of unreduced fibers (via semialgebraic triviality).
  (Route: Hardt again, local triviality over each stratum. Same dependency
  as finite stratification; prove them together.)
- Unified first-order fiber-dimension formula, magnitude and cancellation
  regimes together: dim_1 = 2|S| - g - rank(D) with D the block-projected
  differential (docs/unified_fiber_dimension.md section 1). The
  cancellation-regime evaluation is the proved lemma above; the
  magnitude-regime evaluation (each independent loop contributes one weight
  and one phase direction; touched rigid 1-term classes subtract) is
  corpus-verified but needs a written proof. This is T1's first-order form.

🔬 Universal spectral-distance to fidelity lower bound. NOT yet a theorem
candidate: no precise statement exists in the repo (direction, distance
measures, and quantifiers unfixed). Candidate route: contractivity of the
partial trace + Weyl/Lidskii eigenvalue perturbation + Fuchs-van de Graaf.
Promote to 🟨 only once the inequality is written with quantifiers.

These should not be claimed as established until formal proofs are written
into the project.

### 🟦 Confirmed by project evidence

- Holonomy fields are 2-elementary, census-wide. All 82 gauge-invariant loop
  holonomies of all 63 loopy interference states have cosine of degree 1 (37),
  2 (32) or 4 (13) over Q, with Galois group C1, C2, V4 respectively; no C4,
  D4 or Q8 anywhere. Groups computed with sympy galois_group, not inferred
  from the polynomial having only even terms; loop-kernel bases certified
  saturated so nothing is missed; every PSLQ relation recognized at 120 digits
  and re-verified at 240 with the residual recorded (worst 1e-251).
  PRECISION PROTOCOL, now mandatory for this area: a recognized integer
  relation must be re-verified at at least double the working precision and
  the residual recorded; a relation that fails re-verification is DELETED, not
  downgraded. An exploratory scan at 60 digits produced a false counterexample
  without it. [scripts/holonomy_field_scan.py;
  results/data/holonomy_fields.json; tests/test_holonomy_fields.py]
- The theorem's hypotheses hold census-wide, checked rather than assumed: all
  84 one-hop channels of the loopy states have rational |rho_AB|^2, and for
  all 63 states the within-channel difference vectors generate the loop
  lattice, saturated over Z. The latter is a HYPOTHESIS, not a theorem: a
  support could carry a loop with no one-hop structure, and none does here.
  [results/data/holonomy_fields.json, channel_hypotheses]
- The TRISECTED phase tier is a printing artifact, not a phase-complexity
  tier. The tier is assigned syntactically ("atan" and "/3" both present in
  the printed closed form); the written phase (pi + 3*atan(t))/3 equals
  pi/3 + atan(t), and across all 13 states every representative amplitude has
  2-power degree over Q while every loop holonomy has degree 1, 2 or 4 with
  2-elementary group. No cubic subextension exists anywhere in the class.
  This is stronger than "the cube roots are a gauge artifact"; there was never
  any cubic algebra to be a gauge artifact of.
  [results/data/holonomy_fields.json, trisection_audit]
- Denominator minting is a two-stage arithmetic. The incidence solve
  A_S p = lambda, sum p = 1 and then the pinning solve by the channel
  conditions; the state denominator's new primes come from the second stage.
  At v89 the single silence condition restricted to the incidence line loses
  its quadratic terms and the surviving linear equation is literally
  130*p9 - 1 = 0, minting the prime 5. At v103 the same cancellation happens
  BETWEEN two channel equations, whose difference is 13*p12 = p9 + p13,
  minting 13; the incidence solve there is unimodular over the spectrum
  denominator 34 and mints nothing.
  [scripts/denominator_arithmetic.py; results/data/denominator_arithmetic.json]
- Exactly ONE of the 799 certified records has non-rational squared weights:
  the shipped v_B record at (4,9) index 65, whose weights are quadratic
  irrationals in Q(sqrt 35). It is a genuine non-rational-weight state, not a
  formatting inconsistency: it is the endpoint of a wall family cut out by the
  quadratic 85t^2 - 70t - 119, whose two roots are a Galois conjugate pair
  with both physical and one shipped, so the distinguished point has orbit
  size 2. Its own
  loop holonomy is a sign (minimal polynomial x + 1), consistent with the
  record being real up to gauge. It is excluded from denominator predictions,
  having no state denominator. [docs/holonomy_rationality.md;
  tests/test_holonomy_fields.py]

- Positive-dimensional REDUCED vertex fibers exist (v_B). [THEOREMS.md
  "v_B ANOMALY RESOLVED" (reduced spectrum-fiber is a surface);
  RESEARCH_ARCHIVE "THE FIBER IS THE 1-RDM -> 2-RDM INFORMATION GAP";
  scripts/fiber_sampler.py, scripts/vb_fiber_ideal.py.]
- Pair occupations separate every weight-sector direction on every shipped
  support. For a support S, let B_S be its pair-by-determinant incidence
  matrix. The identity sum_{j != i} (B_S w)_{ij} = (N-1)(A_S w)_i gives
  ker B_S contained in ker A_S, so full column rank of B_S makes pair
  occupations injective on all squared amplitudes and, in particular, on the
  incidence kernel. Exact nonzero maximal minors certify full column rank for
  799 of 799 supports, spanning 403 transport classes. The positive-kernel
  population is 63 records in 39 transport classes, with zero failures. This
  is a support-restricted fixed-frame statement, not a universal theorem:
  signed combinatorial 2-trades are abstract counterexamples, and pair
  occupations are not invariants of the full degeneracy-stabilizer quotient.
  [docs/fiber_symmetries.md; scripts/fiber_symmetries.py;
  results/data/fiber_symmetries.json;
  scripts/verify_fiber_symmetries_standalone.py]
- The first exact full-one-body-unitary invariant on the v_B family is now
  explicit. On its known one-parameter fixed-spectrum slice, with parameter t,
  tr(Gamma_2^2) = 1438/529 is constant but
  tr(Gamma_2^3) = 6(2899-8t)/12167. The cubic moment therefore recovers t and
  proves that distinct t-values survive as distinct points of the reduced
  fiber. It also separates the two Q(sqrt(35)) Galois wall states by
  1728*sqrt(35)/1034195, so Galois conjugacy is not one-body-unitary
  equivalence. Complex conjugation acts by theta -> -theta and fixes the two
  real walls modulo orbital-phase gauge. Whether the two interior conjugate
  sheets at fixed t are equivalent under the full U(1) x U(4) x U(4)
  stabilizer remains OPEN. At the exact interior point t=0, the full
  support-restricted fixed-spectrum constraint differential has rank 7 in 16
  real amplitude coordinates, its orbital-phase gauge has rank 7, and the
  support tangent modulo that phase gauge has dimension 2. Appending
  d tr(Gamma_2^3) and
  d tr(Gamma_2^4) raises the exact rank from 7 to 9; a nonzero 9 by 9 minor is
  299925504*sqrt(8211)/1035662780341225. Thus the cubic and quartic moments
  form a basis of the phase-gauge-quotient cotangent at t=0. Since both are
  full U(9) invariants, no additional one-body-orbit direction lies inside
  this support tangent. This is FIRST ORDER only:
  it does not prove integrability, germ smoothness, or global injectivity on
  the second surface direction.
  [docs/fiber_symmetries.md; results/data/fiber_symmetries.json;
  tests/test_fiber_symmetries.py]
- Cancellation rigidity exists (v89): the certified state is first-order
  fixed-SPECTRUM rigid with no touched 1-term class; its single silent
  channel is within-block and has full rank on the 1-dim incidence kernel;
  full tangent kernel = gauge (dim 9), verified independently 2026-07.
  v103 is NOT a second instance: it is first-order rigid only in the
  fixed-rho sense; three of its five silent channels are cross-block and
  the certified state is fixed-spectrum DEFORMABLE, with 6 non-gauge
  first-order directions (none of them in-span stabilizer-orbit directions)
  and three independent deformation curves continued to t = 0.10 at
  residual ~1e-17 under the exact (rho-aI)(rho-bI) = 0 certificate.
  [docs/silence_rank_lemma.md;
  scripts/reproduce_next_claims.py::v103_deformability; archive P4 sections
  annotated accordingly.]
- Corpus-wide first-order tangent census (all 156 interference states):
  tangent dimension = (surviving weight kernel) + (surviving free phases);
  the loopy/loop-free split matches the archive loop census exactly; v89 is
  the lone crossover of the two regimes (kdim 1, tangent 0); v_B's measured
  dim 2 independently reproduces the archive's "reduced spectrum-fiber is a
  surface". [docs/unified_fiber_dimension.md section 2.]
- v103 sparsification cascade: on-fiber continuation reaches exact-spectrum
  states of support 13 and then 12 (residuals 3.3e-16, 3.9e-15), so the
  certified support-14 state is NOT isolated in the fiber. The fiber germ at
  the certified point is a smooth 6-manifold whose closure is stratified by
  support. [docs/unified_fiber_dimension.md section 3.]
  SUPERSEDED 2026-07: the floor at 12 was a property of that continuation
  path. With a geometric schedule and fallback drop targets the cascade
  continues to support 10, and that endpoint is EXACT, not numerical: all
  ten squared weights are rationals (13/34, 5/34, 53/442, 195/1802, 5/68,
  5/68, 35/901, 1/34, 18/901, 84/11713), the single loop holonomy is pi so
  the state is real up to gauge, and det(rho - xI) = (x - 9/17)^4
  (x - 5/34)^6 holds in exact rational and radical arithmetic. So
  s_Q^free(v103) <= 10, certified. Its state denominator 46852 = 2^2*13*17*53
  carries a prime found in neither the spectrum denominator nor the
  library state, which is why no grid search over multiples of the natural
  denominator could reach it. Upper bound only, not a minimality claim.
  [docs/cascade_census.md; scripts/v103_endpoint_precision.py;
  scripts/v103_support10_certificate.py; tests/test_v103_support10.py;
  results/data/v103_support10_certificate.json]
- The natural-orbital GAUGE was never explored before 2026-07, and exploring
  it does not move the library column. Every vertex has a degenerate spectrum,
  so every certified state carries a gauge freedom U(m_1) x ... x U(m_k) whose
  every point has 1-RDM identically diag(lambda) and is an equally legitimate
  natural-orbital representative. Minimizing support over that family:
  0 of 799 states sparsify, and 670 of 799 are EXACTLY gauge-minimal, certified
  by a gauge-invariant lower bound. A block-diagonal unitary preserves each
  determinant's block profile. Within a profile class, exterior/Koszul
  contraction maps give support at least
  ceil(rank(C_q)/product_i binom(k_i,q_i)). The symbolic rank certificates
  strictly strengthen 82 records and newly close 41, and all 799 shipped bounds
  are exact (floating ranks only select the contraction map). The other 129 are
  open in both directions. The search's power is pinned by a recovery test (a provably
  rank-1 class must collapse to one determinant; it does, at residual 4e-16),
  and two earlier search methods that failed that test were removed rather than
  shipped.
  [docs/gauge_minimization.md; docs/prereg_gauge_minimization.md;
  src/gpc_census/gauge.py; results/data/gauge_min.jsonl]
- NONE of the 71 cascade endpoints is in a natural-orbital basis: 0 of 71 have
  a diagonal 1-RDM. The fixed-spectrum continuation preserves spec(rho), which
  is all it promises, and not diagonality. So every improved cascade bound is a
  bound on s_Q^free and NEVER on s_Q^NO, and none of them may be quoted beside
  s_I. At v103 the exact certificate is arithmetically sound but was labelled
  with the wrong quantity: its state carries rho_39 = -5/68 exactly, so it
  certifies s_Q^free(v103) <= 10 while s_Q^NO(v103) stays bounded by the
  library support 14. [results/data/orbit_check.json;
  scripts/v103_support10_certificate.py; tests/test_v103_support10.py]
- The cascade is mostly an ORBITAL ROTATION ENGINE, but not entirely. By the
  2-RDM spectrum (invariant under every one-body unitary), 63 of the 71
  endpoints ARE the certified library state in different orbitals (agreement
  ~1e-16), so their supports describe the STATE, not the vertex: the
  natural-orbital basis is not the sparsest basis for those states. The other
  8 differ by 3.4e-3 to 4.7e-2, which no one-body rotation can produce, so
  those are genuinely NEW extremal states over their vertices ((4,10) v93,
  (5,10) v65, v256, v144, v140, v230, v263, and (4,9) v65 = v_B). The split is
  bimodal with nothing between 1e-12 and 1e-3. At v103 the connecting unitary
  is exhibited explicitly (residual 7.5e-14) and its off-block entries reach
  0.196, so it is a U(10) rotation and not a gauge one.
  [results/data/orbit_check.json; scripts/orbit_check.py]
- Census-wide sparsification: the support column is an UPPER BOUND on
  s_Q^free, and it is STRICT at 71 of the 799 certified states,
  all INTERFERENCE. The cascade removes 89 amplitudes in total: 55 states
  drop one, 15 drop two, v103 drops four. All 71 endpoints are exactly
  recognizable, but on the SPECTRUM locus only: under the strengthened
  criterion (rho diagonal AND diag = lambda) 0 of 71 certify as attainers, 63
  are the library state in rotated orbitals and 8 are new states bounding
  s_Q^free. The final two exactifications, (5,10) v144 and v256, share
  Q(sqrt(1065)); their exact characteristic polynomials and nonzero
  off-diagonal 1-RDM entries are verified symbolically. [docs/cascade_census.md;
  src/gpc_census/cascade.py;
  scripts/census_cascade.py; scripts/exactify_cascade.py;
  results/data/cascade.jsonl; results/data/cascade_exact.json]
- Rigidity predicts sparsifiability across the census: every state with
  fixed-spectrum tangent 0 fails to sparsify (726 of 726) and 71 of the 73
  with positive tangent do, the two exceptions flooring at 2 to 3e-6, which
  is blocked-on-this-path rather than blocked. This is a POST-HOC pattern in
  the corpus the cascade itself explored; it needs a pre-registered
  out-of-sample test before it is believed. [docs/cascade_census.md]
- The fixed-rho / fixed-spectrum gap is NOT confined to the cancellation
  regime where it was discovered. It occurs in the magnitude regime at v_B
  itself (fixed-rho 1, fixed-spectrum 2) and at 73 of the 799 census states
  (41 transport classes), all INTERFERENCE. Any text presenting the two-fiber
  distinction as a cancellation-regime phenomenon is too narrow.
  PROVENANCE NOTE (2026-07): the (3,11) v43 instance once cited here as
  out-of-sample evidence is NOT independent; it is (3,10) v108 padded with an
  empty orbital, identical on every structural field. The claim stands on the
  census, not on the bracket.
  [results/data/cascade.jsonl; results/data/two_block_family.json;
  docs/prereg_bracket_3_11_3_12.md PROVENANCE CORRECTION]
- A CROSS-BLOCK one-hop channel is NECESSARY for the two fibers to disagree:
  gap > 0 implies X > 0, on all 799 states across 403 transport classes with
  ZERO exceptions. This promotes the silence-rank lemma's cross-block clause
  from a dimension formula verified at two states to a census-wide necessary
  condition. The CONVERSE IS FALSE: 82 states have a cross-block channel and
  no gap, every one of them with incidence kernel 0 and rigid in both fibers,
  so the condition exists but has no weight freedom to cut. The counting
  rules (gap == X, gap == 2X) are both FALSIFIED.
  [docs/two_block_family.md; docs/prereg_two_block_inversion.md;
  results/data/two_block_family.json; tests/test_two_block_family.py]
- The v103 inversion is not unique. 14 states are fixed-rho rigid with
  positive fixed-spectrum tangent, in 7 transport classes, patterns
  (0,1) x 11, (0,2) x 2, (0,6) x 1. v103 is the unique inversion of
  MAGNITUDE 6, not the unique inversion. Its six directions decompose as
  3 weight + 3 phase, reproducing the lemma's formula
  (5-2) + (14-9-2) term by term, with the sectors decoupling exactly as the
  proof's real-state hypothesis requires.
  [results/data/two_block_family.json, tangent_sector_split]
- Two blocks is NOT the operative invariant for fiber behaviour. The gap
  population spans block counts 2 through 6 and peaks at 4; two-block
  vertices are UNDER-represented (3 of 64 carry a cross-block channel against
  155 of 799 census-wide), and all 60 two-block DESIGN vertices have zero
  channels because design supports are one-hop free by definition. The four
  two-block interference vertices span the known behaviours because they span
  (X, dim K), not because two blocks causes anything.
  [docs/two_block_family.md]
- The cancellation regime (exactly diagonal 1-RDM with silent channels) has
  exactly TWO members across all 799 certified representatives, v89 and
  v103; the census-wide cascade adds none, so the silence-rank lemma's
  genericity question stays on a sample of size two. The regime is a
  property of the REPRESENTATIVE, not the vertex: v103 leaves it by cascade
  round 4. [results/data/cascade_summary.json]
- Real boundary states can exist inside fibers whose generic/dense
  representatives are complex. [Both rank-10 closures are all-real with
  rational squared weights, found by reweighted-L1 continuation after the
  contraction attack returned only dense points; RESEARCH_ARCHIVE "RANK-10
  CLOSED ... the near-theorem is FALSE".]
- Kernel-dimension = silent-channel-count coincidence at both closure
  vertices (1 = 1 at v89, 5 = 5 at v103); the lemma shows the v103 instance
  is a fixed-rho fact. [Same artifacts as above; see Problem 2b for the
  open genericity question.]
- Coherence obstruction at the incidence optimum exists (the witness-level
  content of "positive realization defect"), with in-repo certificates
  (scripts/reproduce_next_claims.py). DEFINITION (adopted): for a vertex
  spectrum lambda, s_I(lambda) = min |S| such that the diagonal incidence
  system A_S p = lambda, p >= 0, sum p = 1 is feasible; s_Q^NO(lambda) = min
  support of a pure state whose 1-RDM is EXACTLY diag(lambda); realization
  defect = s_Q^NO - s_I. (The defect is against the NATURAL-ORBITAL minimum.
  Taking it against s_Q^free instead, the minimum over all states with
  spec(rho) = lambda, is not meaningful: s_I does not bound s_Q^free. The
  original wording defined s_Q by the spectrum and then subtracted s_I from
  it, which is the conflation the support taxonomy above exists to stop.) CERTIFIED: s_I(v89) = 7 and s_I(v103) = 6, exact (MILP optimum
  + exact rational witness weights); AND both minimal incidence witnesses
  are PROVED quantum-infeasible in exact arithmetic: v89's by a
  projector-algebra contradiction ((P^2)_04 = P_08 P_84 with both factors
  forced nonzero by singleton channels while P_04 = 0), v103's by
  exhausting the eigenvalue-assignment branches (scalar-block kill;
  multiplicity kill; and the surviving rank-one branch's necessary
  magnitude system has zero nonnegative real solutions, exact solve). The
  defect hierarchy at the closures now reads
      v89 :  s_M = 4 <= s_Q^free <= 10;  s_I = 7 <= s_Q^NO <= 10
      v103:  s_M = 4 <= s_Q^free <= 10;  s_I = 6 <= s_Q^NO <= 14
  with s_M the (weak) majorization baseline. Only s_M <= s_Q^free and
  s_I <= s_Q^NO hold a priori. v89's upper bound serves both columns because
  its certified state is diagonal and did not sparsify; v103's support-10
  state serves only the free-basis column (its 1-RDM is not diagonal), and no
  gauge rotation improves its natural-orbital column.
  [docs/unified_fiber_dimension.md section 3; docs/gauge_minimization.md.]

🔬 STILL OPEN, do not overclaim: the GLOBAL statement s_Q^NO > s_I (defect > 0)
requires killing every support of size s_I..s_Q^NO-1, not just the minimal
witnesses. Ladder status for v89 at size 7: 41 supports enumerated
(incomplete), 34 killed exactly by the projector-pattern argument, 6
distinct survivors all floor numerically (1.5e-2..6e-2 over 24
multi-starts), exactifiable by the same branch analysis that killed the
v103 witness; sizes 8..9 untouched. Completion path: enumerate
incidence-feasible supports per size (MILP no-good cuts) and automate the
projector-pattern + branch-analysis kills (SOS pipeline,
docs/prior_art_roadmap.md section 7). CAVEAT on the baseline: Schur-Horn
permits diag(rho) majorized by lambda, so the airtight
necessary-condition baseline is s_M(lambda) = min |S| with A_S q majorized
by lambda, s_M <= s_I; the defect claim against s_M is strictly stronger
and untested.

### ❌ Disproved

- Holonomy radicands are square roots of rational MONOMIALS in the weights.
  This was the stronger sub-clause of the Holonomy Rationality sketch; it
  fails on 15 of the 82 loops, whose radicands (42, 235, 517, 3266, 4899,
  5538) lie outside the group generated by the weights and the denominator in
  Q*/(Q*)^2. Mechanism: at a channel class of size 3 or more the elimination
  adjoins sqrt(alpha^2 + beta^2 - delta^2), a POLYNOMIAL and not a monomial in
  the weight data, and later loops see the sines of earlier ones. The
  mechanism predicted exactly where the failures must live (only at states
  with a class of size >= 3) and that prediction is confirmed 63 of 63. The
  conjecture's actual conclusion, that the field is multiquadratic, is
  untouched. [docs/holonomy_rationality.md;
  results/data/holonomy_fields.json, radicand_weight_group]
- Ratio-1 states have unimodular incidence solves (P-DEN-2, pre-registered).
  FAIL, 18 exceptions. Unimodularity of A_S is SUFFICIENT for the incidence
  solve to stay over the spectrum denominator but not NECESSARY: an invariant
  factor f > 1 obstructs only right-hand sides outside the image sublattice,
  and (lambda, 1) can lie inside it anyway. P-DEN-3 (the converse) also
  failed, as pre-registered, with its predicted shape: the single exception
  v103 is pinned. P-DEN-4, the version with the escape route named in advance,
  PASSED on 736 states with 0 exceptions.
  [docs/prereg_denominator_arithmetic.md;
  results/data/denominator_arithmetic.json]
- Fiber dimension = incidence-kernel dimension, AS A UNIVERSAL LAW. [P4,
  both closures.] SCOPE NOTE: the law was fit on magnitude-target states
  and fails in the cancellation regime; the scoped (magnitude-sector)
  version is the surviving candidate and is exactly the unobstructed sector
  of the T1 Jacobian program. Do not discard the scoped version with the
  universal one.
- Singleton channels are the only rigidity mechanism. [P4: both closures
  rigid via multi-term cancellation constraints, no 1-term class. v89's
  rigidity is fixed-spectrum; v103's is fixed-rho.]
- Universal <=2-channel interference law. [P2: v103 carries 5 exchange
  channels; independently re-verified.]
- State denominator = spectrum denominator (and its ratio-in-{1,2}
  refinement). [P3: ratios 5 (v89) and 5746 (v103), rational-weight states,
  no escape clause. LEDGER NOTE 2026-07: the regenerated ledger ships all
  12 DESIGN-REAL states at ratio 2 (the archive's mined "9 of 12" is stale;
  representative-dependent, as the P2/P3/P4 caveat predicts).]
- Complex representative implies intrinsically complex fiber. [Both
  closures real; archive "the near-theorem is FALSE".]
- Positive-dimensional fibers are literal torus orbits. [Archive
  retraction: rational curve fibered over the gauge torus, not a torus
  orbit.]
- Generic one-loop fibers are elliptic. [Archive retraction; census of
  kdim-1 families: 18 genus-1 vs 26 genus-0.]
- Dense numerical representatives imply no sparse exact representative.
  [The contraction attack lands on dense generic fiber points; reweighted-L1
  over the full fiber found support 10 and 14; archive "WHAT OVERTURNED THE
  NEAR-THEOREM". The v103 cascade sharpens this: even a certified sparse
  point can fail to be the sparsest in its fiber.]
- Two spectral classes require higher-order interference, under the ADOPTED
  reading "spectra with exactly two eigenvalue blocks force complex (or
  higher-degree algebraic) interference phases". DISPROVED decisively by
  ledger census (scripts/reproduce_next_claims.py): ALL 64 two-block
  vertices ship all-real certified states, including all four two-block
  interference vertices (v89, v103, v108, and one (5,10) vertex); trisected
  (cube-root) phases occur only at 3-5 blocks. Existence of a real state
  refutes necessity, so this direction is representative-safe.
- Spectral border support differs from exact support, under the ADOPTED
  reading "the minimal support attainable as a LIMIT of exact-spectrum
  states can be smaller than the minimal exactly-attained support".
  DISPROVED at the tested instance (v89): the certified support-10 state is
  the c_(036) -> 0 endpoint of a support-11 one-parameter exact-spectrum
  family (continuation verified to residual ~1e-17 over c_x in [0, 0.2];
  the other five first-order openers are second-order obstructed), and the
  endpoint IS attained; border = exact here. NOTE: both adopted readings
  are reconstructions; if the parallel instance's original definitions
  differ, re-tag against those. Formal definitions of support-rank and
  border-support-rank, with the tensor rank vs border rank analogy,
  are in docs/prior_art_roadmap.md section 5.

## Highest-value open problems

(Each problem cites its existing groundwork; new agents start there, not
from zero. docs/prior_art_roadmap.md maps each to its outside-field sibling
and names the tools to import.)

0. 💡 TWO-TORSION IN THE FIBER. Relocated here from the manuscript outlook
   before the freeze: it is an interpretation resting on a single vertex, so
   it belongs in the open-questions register and not in the paper. The two
   counts it rests on ARE certified and stay in the manuscript's results
   section as data (results/data/vb_galois_orbit.json). Verbatim as it stood:

   > One observation ties the arithmetic to the geometry, and we flag it as
   > speculative because it rests on a single vertex. At v_B the amplitude
   > signs that attain the vertex are not arbitrary: for each of the two
   > Galois embeddings, exactly 64 = 2^6 of the 128 sign patterns modulo
   > global sign attain, so the attaining set is a 2-torsion coset over a
   > single weight vector rather than a scattered subset [E]. That is the
   > Kummer structure of the holonomy fields, every field in the census being
   > multiquadratic, hence generated by square roots, appearing geometrically,
   > as 2-torsion in the fiber over a rational point. Whether the fixed-rho
   > fiber of a rigid vertex is in general a finite 2-torsion set over a
   > unique weight vector is open, and we do not assert it: the one other
   > vertex where multiple isolated fiber points have been exhibited, the
   > denominator-34 vertex, has points that share a single weight vector but
   > are separated by holonomies taking values across (-1,1) rather than by
   > signs, so it does not fit the picture as stated.

   Status: the v_B leg is certified. The v103 leg no longer refutes anything,
   because v103's fiber has exactly ONE point modulo gauge, holonomy
   (0, pi, pi, 0, 0), certified by exhausting the 5-torus of gauge-invariant
   holonomies (16807 grid starts, all converging, one point). The earlier
   report of dozens of points at v103 was an artifact of a non-integer loop
   basis and is withdrawn; see correction 11 in docs/closing_ledger.md. So the
   open question is not whether v103 fits the 2-torsion picture, it is whether
   the picture says anything at all beyond v_B, where it is a single certified
   data point. Attack it as a mechanism, not by mining two vertices.
   [results/data/vb_galois_orbit.json; results/data/v103_fiber_count.json] So whatever unifies the two vertices is not
   2-torsion in the holonomy. This is arc two's first question, and it should
   be attacked as a mechanism, not mined from the two data points.
   [results/data/vb_galois_orbit.json;
   results/data/uniqueness_census.json, v103_discreteness_probe]

1. Jacobian formula for local fiber dimension. [= THEOREMS.md T1. Proof
   route on file: MGS local normal form in the Sjamaar-Lerman stratified
   setting; the obstruction quadric is computable per state. First-order
   form now on file with corpus verification
   (docs/unified_fiber_dimension.md); the cancellation regime is proved
   (silence-rank lemma); remaining: the magnitude-regime proof and the
   second-order/integrability structure.]
2. Explicit reduced-fiber invariant for v_B, and exactification of the v103
   fiber. [PROGRESS: on the known fixed-spectrum v_B slice,
   tr(Gamma_2^3) = 6(2899-8t)/12167 exactly recovers the family parameter and
   separates the two Galois walls under the full one-body group. This gives
   one exact reduced coordinate. At t=0 its differential together with
   d tr(Gamma_2^4) exactly separates the complete two-dimensional support
   tangent modulo orbital-phase gauge. Global integration and separation of
   the second surface direction,
   and the full-stabilizer status of the conjugate phase sheets, remain open.
   Groundwork: docs/fiber_symmetries.md, the v_B anomaly resolution,
   holonomy/wall data, vb_holonomy.py, vb_fiber_ideal.py. New v103 lead: the
   (0,2,7) weight
   equals 5/34 (the small spectrum value itself) at the certified point and
   at every cascade point to err < 3e-17, a candidate exact invariant and a
   coordinate to quotient out before computing the fiber ideal; PSLQ over
   quadratic surd lattices, seeded from the support-12 endpoint, is the
   next tool (docs/unified_fiber_dimension.md section 4). Witness sets +
   numerical irreducible decomposition + monodromy Galois groups are the
   industrial method (prior_art_roadmap section 4).]
   2b. Silence-rank genericity. The lemma is PROVED
   (docs/silence_rank_lemma.md); what remains open: the second-order
   structure at v103 (which of the 6 tangent directions integrate; at least
   3 do), an exact constant-rank argument upgrading the numerical
   germ-smoothness, and whether within-block silence rank = kernel dim is
   forced for census cancellation states or accidental (sample of size 2).
   The rigidity-matroid formulation (coherence rigidity matroid,
   generic-rank theorem or counterexample search) is the imported frame
   (prior_art_roadmap section 2).
3. Interior inverse walls. [Groundwork: wall_solver.py, wall_test.py,
   per-family endpoint theorem (T4).]
4. Graph-theoretic coupled-phase realizability criterion. [Groundwork:
   one-hop class combinatorics in classify.py; loop census in the archive.]
5. Infinite positive-defect family. [The defect definition is now fixed and
   witnessed (s_I / s_Q^NO / s_M above); still BLOCKED on the global defect
   inequality at one instance (the 🔬 ladder above) before hunting a
   family. Crystallographic sparse phase retrieval is the template
   (prior_art_roadmap section 1).]
6. Infinite deformable and cancellation-rigid families. [Candidate seeds:
   the 44 kdim-1 deforming families; the two closure states.]
7. Levi-theoretic lower bounds. [Groundwork: levi.py post-audit, the
   rank_deficiency == stab - 1 identity (799/799), the hard/soft plethysm
   criterion.]
8. Holonomy Rationality at channel classes of size 3 or more. [Theorem A
   settles size <= 2; Proposition B gives only a 2-tower, and a 2-tower can be
   dihedral (sqrt(1 + sqrt 2) has D4 closure). What is needed: an argument
   that the eliminated discriminants alpha^2 + beta^2 - delta^2 are always
   squares times elements already in the field, or a structural reason the
   towers cannot go dihedral. A single C4 or D4 holonomy anywhere settles it
   the other way; none occurs in 82 loops. Second open piece: is hypothesis
   (H3) forced, i.e. can a vertex support carry a loop with no one-hop
   structure? docs/holonomy_rationality.md.]
9. Complexity classification. [BLOCKED on choosing the complexity measure:
   support size, state-denominator, phase field degree, and clique budget
   are all in play and provably non-equivalent (den-ratio falsification).
   Define the measure first.]

## Dependency graph

Incidence geometry
    -> Weight polytope
    -> Coherence equations
    -> Jacobian
    -> Fiber dimension / rigidity
    -> Inverse walls
    -> Support complexity atlas

## Dead ends

Future agents should NOT assume:

- kernel dimension predicts fiber dimension (outside the magnitude-target
  scope);
- representative properties characterize an entire fiber (P2/P3/P4 and the
  projective-stabilizer demotion are all representative-dependent);
- fixed-rho rigidity implies fixed-spectrum rigidity (v103; name the fiber
  in every rigidity claim);
- the certified ledger support is the minimal support (it bounds s_Q^NO from
  above; the natural-orbital gauge does not improve it anywhere in the census,
  and it is NOT known strict at v103, whose support-10 witness is a
  free-basis one);
- that a support bound may be compared with s_I without naming which minimum
  it is (s_I bounds s_Q^NO only, never s_Q^free);
- that a sparser state reached by continuation is a NEW state (63 of 71
  cascade endpoints are the library state rotated; check the 2-RDM spectrum);
- toric-orbit descriptions;
- denominator preservation;
- two-channel universality;
- that a dense numerical fiber point says anything about sparse exact
  points (the contraction-attack artifact that nearly closed the census
  wrongly).

## Immediate tasks

Done and landed (2026-07): scripts/reproduce_next_claims.py (exact s_I
certificates, witness infeasibility proofs, two-block phase census, v89
support-11 family, v103 deformability); the silence-rank lemma proof
(docs/silence_rank_lemma.md); the corpus tangent census and v103 cascade
(docs/unified_fiber_dimension.md).

Also done and landed (2026-07, fiber session): the two fiber notions made
structural in code (src/gpc_census/fiber.py, tests/test_fiber_notions.py);
the census-wide sparsification cascade (src/gpc_census/cascade.py,
scripts/census_cascade.py, docs/cascade_census.md); exact certification of
the improved support bounds (scripts/exactify_cascade.py); the scored
out-of-sample bracket pre-registration (docs/prereg_bracket_3_11_3_12.md);
and the reconciliation artifact (scripts/reconcile_claims.py).

Also done and landed (2026-07, gauge session): the natural-orbital gauge
minimization and its exact lower bound (src/gpc_census/gauge.py,
scripts/gauge_census.py, docs/gauge_minimization.md); the 2-RDM acceptance gate
and the census-wide orbit check of the cascade endpoints
(scripts/orbit_check.py); the support taxonomy above; and the relabelling of
every cascade bound from s_Q to s_Q^free.

Also done and landed (2026-07, arithmetic session): the census-wide holonomy
field scan (scripts/holonomy_field_scan.py,
results/data/holonomy_fields.json, tests/test_holonomy_fields.py); the
Holonomy Rationality theorem note with Theorem A proved and the conjecture
scoped (docs/holonomy_rationality.md); the two-stage denominator derivation
with its pre-registered census test (scripts/denominator_arithmetic.py,
docs/prereg_denominator_arithmetic.md,
results/data/denominator_arithmetic.json); the trisection dissolution; and
the manuscript's Galois section updated from the stale 142-state tally to the
full 82-loop result.

Priority:

0c. Extend the exact global no-design certificate atlas from the 12 numbered
   theorem/proposition cases to the other 144 INTERFERENCE verdicts. The
   disjunctive Farkas method is complete and fast on all 12 published targets;
   a highly symmetric rank-10 two-block case shows that orbit clauses must be
   stored and verified symbolically rather than expanded naively.
0. Continue the gauge lower-bound campaign on the 129 open rows. The
   exterior/Koszul contraction family removed the former Lambda^a fallback,
   improved 82 records, and certified 41 more; the remaining slack needs a
   different invariant.
0b. Exactly recognize the remaining 6 genuinely new extremal states found by
   the orbit check. The v144/v256 Hodge pair is now exact over Q(sqrt(1065)).
1. Formalize the theorem-candidate proofs (the two Hardt-dependent ones
   together).
2. Alpha-certify what exact recognition could not reach. OVERTAKEN for the
   cascade: all 71 endpoints are exactly recognized on the spectrum locus by
   high-precision refinement, degree-at-most-two integer relation detection,
   and an exact characteristic-polynomial check. The v89 (0,3,6) family points
   remain candidates for alpha theory.
   [results/data/cascade_exact.json]
3. Witness set + numerical irreducible decomposition of the v103 fiber and
   the v_B surface; monodromy Galois groups feed the holonomy program
   (prior_art_roadmap action 2).
4. Prove the Jacobian dimension theorem (T1), unifying the two regimes
   (Problem 1); coherence rigidity matroid for Problem 2b's remainder.
5. Complete the v89 defect ladder (sizes 7..9) with the automated
   projector-pattern + branch-analysis + SOS pipeline; then the global
   defect theorem.
6. Integrate the exact cubic-quartic tangent coordinates over the second
   direction of the v_B reduced surface, test global injectivity, decide the
   full-stabilizer status of the conjugate sheets, and exactify the v103 fiber
   (Problem 2; the 5/34 invariant coordinate is the v103 seed).
7. Search for interior inverse walls.
8. Develop practical realization algorithms from the corrected theory
   (reweighted-L1 full-fiber continuation -> PSLQ exactification is the
   validated pipeline; isospectral-flow explorers are the principled
   replacement for ad-hoc continuation; productionize
   scripts/repro_and_sparsify.py-class tooling out of scratch).

## Rigidity and rationality (2026-07)

Full write-up: docs/rigidity_rationality.md. This is where the weight-rationality
hypothesis that docs/holonomy_rationality.md and the denominator arithmetic both
assume gets its grounding.

- ✅ DESIGN CASE PROVED. On a one-hop-free support every off-diagonal 1-RDM
  entry vanishes termwise, so rho = diag(lambda) is exactly the rational
  linear system A_S p = lambda with sum p = 1; when its solution is unique the
  weights are rational by linear algebra alone. All 643 design states have
  incidence kernel 0, so the hypothesis holds throughout and the theorem
  covers every one of them. No rigidity input and no Galois argument needed.
  The 156 interference states are not reachable this way, because phases make
  the system nonlinear in p. [docs/rigidity_rationality.md]
- 🟦 FIXED-RHO RIGID IMPLIES RATIONAL WEIGHTS, 740 of 740, zero
  counterexamples. The single irrational record (v_B, weights in Q(sqrt 35))
  is deformable. Deformability does not force irrationality: 58 of the 59
  deformable records are rational anyway.
  POWER WARNING, not optional when quoting this: all 645 natural-orbital
  representatives are rigid and all 59 deformable records are among the 154
  that fail the diagonality gate, so the discriminating power lives entirely
  in the non-natural-orbital population and there is no rigid-versus-
  deformable contrast among gate-passing records. The count is right and the
  independent content is smaller, exactly as with the bracket holdout.
  [scripts/rigidity_rationality_census.py;
  results/data/rigidity_rationality.json; tests/test_rigidity_rationality.py]
- 🎯 TARGET, not proved: if the support-restricted fixed-rho fiber modulo
  gauge is a single point (true isolation, not first-order rigidity), that
  point is defined over Q. Two open gaps, both stated in the write-up:
  (a) first-order rigidity is not isolation, needing stress-matrix or exact
  isolation certificates; (b) isolation gives finiteness, not uniqueness, so
  the honest form is "algebraic of degree at most the number of isolated
  fiber points modulo gauge, rational when unique". v_B instantiates (b):
  degree 2 over Q, orbit size 2, from the quadratic wall 85t^2 - 70t - 119.
- ❌ CORRECTION: the diagonality gate fails on 154 records, not one. The
  handoff that opened this line expected exactly one violation (v_B) and
  instructed that the expectation be verified; it is wrong. Every
  INTERFERENCE record except v89 and v103 has a non-diagonal 1-RDM. All 154
  still attain the spectrum exactly at unit norm, the library certificate
  never asserted diagonality for them (it is the characteristic-polynomial
  identity), and RESEARCH.md already recorded the same fact from the other
  side in the cancellation-regime count. So this is how the interference
  library was built, not a regression. What genuinely needed repair was the
  support-taxonomy claim, corrected above.
  [scripts/attainer_gate_audit.py; results/data/attainer_gate_audit.json]
- 🟦 REPAIRED (2026-07), additively. All 156 interference records now have a
  natural-orbital representative in results/data/states_natural_orbital.jsonl,
  each with 1-RDM exactly diag(lambda) (worst off-diagonal 2.8e-16) and each
  passing the gate, re-verified from the shipped amplitudes. The 2-RDM
  spectrum is unchanged to 1.6e-15, so the rotation moved the basis and not
  the state. Convention note: the compound is of U^T, NOT U^H; they agree for
  real rotations, so a wrong convention fails only on the complex records.
  THE COST: support grows on 142 of 156 (mean 6.8 -> 12.1, max 14 -> 79) and
  the exact closed forms do NOT survive the rotation, so the new ledger is
  NUMERICAL. states.jsonl is therefore NOT replaced. The two ledgers bound
  different quantities: the census column bounds s_Q^free, the new ledger
  bounds s_Q^NO, and any s_Q^NO claim must quote the latter.
  [scripts/natural_orbital_ledger.py;
  results/data/states_natural_orbital.jsonl;
  tests/test_natural_orbital_ledger.py]
- 🟦 CLOSED, PREDICTION CONFIRMED: the v_B Galois embedding test. Both roots
  of the wall polynomial 85t^2 - 70t - 119 give positive weight vectors that
  attain the vertex EXACTLY, each by 64 of 128 sign patterns modulo global
  sign. Orbit size 2, weight field degree 2: degree equals orbit size, which
  instantiates gap (b) of the Rigidity-Rationality target. The test resolves
  the Galois action on the wall parameter t rather than on nested surds,
  which is precisely what the naive attempt got wrong.
  [scripts/vb_galois_orbit.py; results/data/vb_galois_orbit.json]
- (superseded, kept for the record) PREREG-PENDING: the v_B embedding test. The naive substitution
  sqrt(35) -> -sqrt(35) fails attainment (spectrum error 0.068), and that
  failure is DIAGNOSTIC, not refuting: the amplitudes carry nested surds like
  sqrt(120 - 18 sqrt(35)) and independent principal-branch substitutions are
  not a field embedding. Since the fiber equations are polynomial over Q, a
  genuine embedding carries exact solutions to exact solutions, so a
  substitution that breaks attainment proves only that it was not an
  embedding. The proper test (build the field, take its real embeddings,
  apply each consistently, test attainment, orbit size pre-registered) is
  specified and NOT YET RUN.

## Reconciliation ledger (2026-07 fiber session)

Each item is re-derived by scripts/reconcile_claims.py into
results/data/reconciliation.json, so a reader can check it rather than trust
this table.

| superseded reading | corrected reading |
|---|---|
| v103 is first-order rigid | rigid for the FIXED-RHO fiber, deformable (tangent 6) for the FIXED-SPECTRUM fiber; v89 rigid in both |
| v103 is an isolated support-14 point | true of the route-B search landscape and gauge slice, false of the fiber: support 10 is reached and certified exactly |
| the census support column is the minimal support | it is the certified support, an upper bound on s_Q^NO; the cascade's 71 improvements are improvements of s_Q^free, a different quantity |
| the certified support column was never minimized over the natural-orbital gauge, so it is loose | it was never minimized, and minimizing it changes nothing: 0 of 799 sparsify, 670 of 799 are certified gauge-minimal |
| s_Q^NO(v103) <= 10 (the support-10 endpoint read as a natural-orbital state) | FALSE. That endpoint's 1-RDM carries rho_39 = -5/68 exactly, so it certifies s_Q^free(v103) <= 10 only; s_Q^NO(v103) stays bounded by the library support 14 |
| the cascade endpoints are all the library state in another basis | 63 of 71 are; the other 8 have 2-RDM spectra no one-body rotation can reach, so they are genuinely new extremal states |
| 9 of 12 DESIGN-REAL states at denominator ratio 2 | the shipped ledger has 12 of 12; the mined figure is representative dependent |
| the trisected states are a cube-root phase-complexity tier | a printing artifact: (pi + 3*atan(t))/3 is pi/3 + atan(t), and no cubic subextension exists in any of the 13, invariant or representative |
| holonomy radicands are square roots of weight monomials | true only for channel classes of size at most 2 (proved); false at 15 of 82 loops, all at classes of size 3 or more, exactly as the mechanism predicts |
| the manuscript's Galois tally covers 142 states / 62 holonomies | superseded by the full corpus: 63 states, 82 loops, all 2-elementary |

REPRESENTATIVE-SWAP CHECK (v0.12, the real v_B record). Shipping the real
theta=0 endpoint in place of the phase-carrying fiber point over
(4,9) index 65 leaves every fiber-program result unchanged. Re-run on the new
ledger: identical census totals (71 sparsifying states, 89 amplitudes, 33
multiplicity-drift stops, 2 cancellation-regime members), identical per-vertex
verdicts and tangent dimensions for all 799 records, and identical exact
certification (all 71 after quadratic recognition, same numerical endpoints
throughout). v_B's own numbers are
stable too: incidence kernel 1, fixed-spectrum tangent 2, fixed-rho tangent 1
for BOTH representatives, and both cascade to the same support-7 point,
dropping determinant (0,1,2,4) and certifying at the same squared weights
(9/23, 4/23, 3/23, 12/115, 8/115, 1/23, 2/23), state denominator 115. The two
endpoints differ only by a gauge phase, so the real representative simply
lands on the real form directly. The one visible change is elsewhere: v_B's
phase tier in the block-count census moves from COMPLEX to REAL, so the
3-block row reads 152 REAL / 13 COMPLEX rather than 151 / 14
(scripts/reproduce_next_claims.py claim 2).

## Guiding principle

The census should now be viewed as **experimental evidence**.

The primary research product is the mathematical theory of inverse fibers,
lifting complexity, and coherence obstructions. The inverse-fiber program
imports its questions from phase retrieval (ambiguity geometry), its
rigidity from framework theory, its algorithms from structured inverse
eigenvalue problems and numerical algebraic geometry, and its pathology
expectations from tensor rank, applied to a map on which none of those
fields has yet worked (docs/prior_art_roadmap.md).

## Exact RDM observability atlas (2026-07)

The census now has a second role: a certified sparse correlator-reconstruction
atlas with a phase-synchronization interpretation.
For 771 of 799 shipped determinant supports, including all 156 INTERFERENCE
supports, the 2-RDM rationally reconstructs every pure full-support state on
the fixed carrier and therefore every higher RDM. The proof uses two exact
support gates: full column rank of pair incidence and a connected spanning
tree of collision-free 2-RDM coherence channels. Pair incidence has full
column rank on all 799 supports.

This is a fixed-frame, carrier-restricted, pure-state statement. It does not
extend to arbitrary mixed states, infer a carrier, prove non-injectivity for
the 28 unresolved supports, or supply fixed-rho dynamics. Its immediate use
is exact higher-RDM reconstruction and carrier-preserving reduced-dynamics
benchmarks. The high-value dynamic gate, a correlation-changing two-body
trajectory that remains in a certified carrier and is not a one-body orbit,
is met below for \(v_B\).
The remaining gates are an infinite-family and stability theorem.

[docs/rdm_observability_census.md;
scripts/rdm_observability_census.py;
scripts/verify_rdm_observability_standalone.py;
results/data/rdm_observability_census.json;
tests/test_rdm_observability_census.py]

## Exact scheduled v_B carrier dynamics (2026-07)

The main dynamic gate above is now closed for one exact finite model. A
Hermitian pair kernel with a seven-entry symbolic support drives the known
\(v_B\) conic, preserves its complete eight-determinant carrier in the full
126-dimensional four-particle sector, and satisfies the Schrödinger equation as a vector
identity. A diagonal orbital lift gives a second at-most-two-body schedule
whose reference trajectory has one constant complete 1-RDM, not merely a
fixed diagonal or spectrum.

The basis-free invariant
\(\operatorname{Tr}(\Gamma_2^3)=6(2899-8t)/12167\) varies strictly with
\(t\). Therefore the loop is not a one-body orbital orbit, including an orbit
inside the stabilizer of the fixed 1-RDM. Since the disconnected two-body
contribution is fixed by that complete 1-RDM, the changing part is the
two-body cumulant.

The carrier has a unimodular pair-incidence minor and a connected seven-edge
2-RDM coherence tree. Therefore the scheduled flow remains on a chart with an
exact rational \(\Gamma_3(\Gamma_2)\) map and gives exact closed nonautonomous
TD-2RDM evolution. This is a reference-schedule result, not an autonomous
Hamiltonian, a mixed-state closure, or a scalable speedup. The package proves
closure but does not yet implement an amplitude-free reduced propagator. The
next gates are that executable reconstruction, conditioning, additional
nontransport carriers, an infinite family, and a beyond-full-wavefunction
benchmark.

[docs/vb_carrier_hamiltonian.md;
scripts/vb_carrier_hamiltonian.py;
scripts/verify_vb_carrier_hamiltonian_standalone.py;
results/data/vb_carrier_hamiltonian.json;
tests/test_vb_carrier_hamiltonian.py]

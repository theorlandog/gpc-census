# AGENT_PLAN.md: mission queue and sync invariant

Read `AGENTS.md` first; it owns house rules, layout, commands, and versioning.
This file owns WHAT to work on and the drift policy. If they conflict,
`AGENTS.md` wins on engineering conventions, this file wins on task selection.

State moves daily. At session start: `git pull`, re-read `docs/RESEARCH.md`
(current campaigns), the newest `docs/HANDOFF*.md`, and this file. Re-derive
every count from artifacts before quoting it. Never trust a snapshot over the
repo.

## Session protocol

1. Pick exactly ONE task from the queue below, highest eligible priority first.
2. Any run that measures something gets a pre-registration from
   `docs/prereg_template.md` BEFORE the run, and the prereg is scored before
   any narrative is written. Failed predictions are findings; record them.
3. A claim is not a certificate. Every new claim ships with an exact artifact
   in `results/data/`, a generating script, a guard test, and where feasible a
   standalone verifier that imports nothing from `gpc_census`.
4. Later HANDOFFs supersede earlier ones. Never reinstate a retracted claim
   (example: `s_Q^NO(v103) <= 10` and any diagonality assumption on the v103
   support-10 endpoint are dead; see HANDOFF 4 correction A).
5. Session end: update `docs/RESEARCH.md` with tagged findings and artifact
   pointers, run `make emit`, `make test` and `make lint` green, PR per
   template.
6. Writing style: no em dashes, no filler, evidence tags on every claim using
   the existing glyphs and tiers (PROVED / CONFIRMED / OPEN, [E] / [CA],
   and the RESEARCH.md status glyphs).

## SYNC INVARIANT (standing; blocks merge)

Every number that appears in more than one place has exactly one generator.

1. Emitted blocks. Any table or headline count in a doc that restates artifact
   data is machine-emitted between markers
   `<!-- sync:NAME:start -->` ... `<!-- sync:NAME:end -->`.
   `scripts/emit_doc_tables.py` owns them: `--write` regenerates in place,
   `--check` exits nonzero on drift, and `tests/test_doc_sync.py` runs the
   check. Never hand-edit between markers; run `make emit`.
   Current coverage: the README coverage table, the (3,11) bracket tally, and
   the four benchmark tables.
2. Benchmark runs are artifacts, not prose. Multi-machine results are reported
   per machine and never averaged. The emitter states the machine, or states
   that the artifact records none, and carries the unmet
   equally-optimized-control caveat in every benchmark footer.
3. Completeness audit as a test. Every file under `results/data/` must have a
   `SHA256SUMS` entry, a `PROVENANCE.md` row, a named generating script, and
   at least one guard test. `scripts/audit_data_completeness.py` enforces it
   and `tests/test_data_completeness.py` fails the suite on an orphan.
   Relaxations live in one declared `FAMILY_RULES` table with the exact marker
   each waiver depends on, and a rule matching nothing is itself an error.
   Frozen zones (`results/report/anonymized/`) are checked for drift but never
   rewritten by the emitter.
4. Narrative numbers in `docs/RESEARCH.md` are exempt from emission but must
   cite the artifact path inline, which is already the convention.

`make verify-data` runs all three gates.

## One-time cleanup: DONE (2026-08-05)

Superseded, kept as the record of what it found. Do not rerun as a task.

- README coverage table: it already reconciled at 799, but the ledger it
  restated did not. It is now emitted from `census_master.csv`.
- `census_master.csv` had NO generator, was malformed CSV (the system field
  `(3,6)` shipped with its comma unquoted, so every row parsed as seven fields
  against a six-field header), and was missing the entire `(3,8)` block: 761
  rows against 799 vertices, with the pre-v0.11 DESIGN-REAL count of 11.
  `scripts/build_census_master.py` now derives it from `states.jsonl`.
- `docs/CLOSURE_PLAN.md`: the 797/799 header is marked HISTORICAL, and the
  parallel SOS cand-44 item is marked SUPERSEDED, do not run.
- The four benchmark docs are emitted from their artifacts; hand-typed metrics
  are gone.
- Completeness audit: 52 orphans on first run, 0 after. Two artifacts,
  `cyclic_window_benchmark.json` and `entangling_window_hamiltonian.json`, have
  no surviving generator, no doc and no consumer. They are declared in
  `PROVENANCE.md`, waived by name in `FAMILY_RULES`, and pinned by digest in
  `tests/test_orphan_artifacts.py`. OPEN DECISION: regenerate from a recovered
  method, or delete. Do not delete incidentally.

## Priority queue

### T0. Ship-prep for paper 1 (highest; human gates flagged)

- Verify CLOSURE_PLAN steps 2 to 4 are actually complete, with evidence:
  P1-P6 scored in `docs/RESEARCH_ARCHIVE.md`, blog line softened, final-count
  markers gone from `results/report/main.md`.
- Run `make verify-paper`, `make verify-data`, `make report`,
  `make anonymize` + drift check, `make paper1-supplement`. All green.
- Do NOT expand paper scope. Do NOT fold the level-5 proofs into paper 1;
  its frozen open-problem statement was true at census freeze and stays.
- Output: `SHIP_CHECKLIST.md` listing the HUMAN-ONLY steps in order:
  arXiv upload via the standing endorsement, Zenodo release tag, the one
  email (paper, PhD question, spin-fiber direction), EV application.
  The agent prepares and validates; it never submits.

### T1. Companion note: the level-5 certificates. DRAFTED 2026-08-05.

- `results/report/level5/main.md`, built by `make note` on the same pinned
  pandoc image and citation discipline as the paper.
  `tests/test_level5_note.py` re-derives every load-bearing number from the
  artifact and asserts it appears in the text, checks that every bib entry is
  cited and every citation resolves, and asserts the note states what it does
  NOT prove. OUTSTANDING: the PDF has never been rendered (no container
  runtime in the drafting environment), and the bibliography entries for
  Ressayre 2010, Belkale-Kumar 2006 and Berenstein-Sjamaar 2000 were written
  from memory and need a citation check before submission. Human gate: arXiv.
- Original scope, for reference. Short note (4 to 8 pages) from
  `docs/level5_ressayre_3_11.md`,
  `docs/rank11_open_candidate_designs.md`, and
  `scripts/verify_level5_ressayre_3_11_standalone.py`. First proofs of the
  claimed Klyachko level-5 rows; bracket closes 19 TRUE, 31 REFUTED, 0 OPEN.
- Every claim maps 1:1 to a repo artifact. State explicitly that validity is
  proved and facetness, irredundancy, and completeness of (3,11) are not.
- Reuse the report toolchain (new directory beside `results/report/main.md`,
  same pandoc container, same citation discipline). Human gate: arXiv.

### T2. Generator ladder toward the first complete new GPC system since 2008

Strict order, one gate per session minimum, prereg per gate:

1. (3,7). CLOSED TWICE, and no longer a gate. Corrected 2026-08: an earlier
   version of this step read "candidate exhaustiveness is bypassed, not
   proved", which is true of the inner/outer closure alone and made rank 7 look
   open when it is not.
   (a) INNER/OUTER CLOSURE, DONE 2026-08-05, all seven pre-registered
   predictions PASS, prereg frozen at 9c55f1b before the run: 10 vertices
   matching the reference, every row certified on independent replay, every
   vertex carrying an exact attaining state, and controls confirming the four
   rows cut the 13-vertex ordered slice to 10 with each row individually
   changing the polytope. This route bypasses exhaustiveness.
   (b) CANDIDATE EXHAUSTIVENESS, proved separately by the closed-flat
   enumerator (`docs/fixed_n3_generator_methodology.md`): 5,341 hyperplanes in
   place of 1,623,160 candidate bases, 10,682 oriented rows, screening
   resolving every one of them with none unresolved (10,133 trace failures, 499
   exact determinant zeros, 1 structural Pauli row, 49 nonstructural Ressayre
   rows), and exact reduction returning the four published facets. So rank 7
   has system completeness, not just a relative system.
   [docs/generator_gate1_3_7.md; results/data/generator_gate1_3_7.json;
   results/data/stage1_reference_series.json]
2. WEYL-ORBIT CANONICAL ENUMERATOR, before grinding rank-8 flats. Orbit
   reduction is the one pruning the standing rule permits, since the moment
   polytope is Weyl invariant and the ordered slice is a fundamental domain, so
   nothing is discarded. Measured payoff: 362 hyperplanes in 8 orbits at rank 6
   and 5,341 in 19 at rank 7, reductions of 45x and 281x, so the raw count
   explodes while the orbit count barely moves. The right algorithm never
   materializes the flat lattice; it enumerates orbit representatives by
   canonical augmentation.
   [docs/prereg_orbit_canonical_flats.md; scripts/orbit_canonical_flats.py;
   results/data/orbit_canonical_flats.json]
   DONE, and then made cheap enough to matter. Canonicalization was measured,
   not assumed, and its worst per-flat cost was exactly `d!` at every rank,
   because plain colour refinement splits nothing on the most symmetric flats.
   Individualization removes that tail: (3,8) peak worst case 40,320 to 1,440
   and the (3,8) enumeration 401 s to 16.7 s, with leaves now equal to the
   automorphism group order on 300 of 313 flats, which is the floor. Children
   are also STREAMED rather than materialised, so peak memory is the orbit count
   of the next level, which is the shape problem that killed the (3,11) run at
   6.2 GB.
   [docs/orbit_canonical_flats.md; tests/test_orbit_canonical.py]
2b. DIRECT SURVIVOR GENERATION, which removes the row scan. The Ressayre trace
   test is `inv(sigma h) == B` with `B` an orbit invariant, so the survivors of
   an orbit are the arrangements of a multiset with a prescribed inversion
   number, and those are generated directly in output-linear time instead of
   filtered out of the orbit. Verified SET against set on all 35 rank-7 orbits,
   549 of 10,682 rows, 133x. This is what
   `docs/screening_orbit_structure.md` named as the step that sets the reachable
   rank.
   [docs/prereg_direct_survivor_generation.md;
   docs/direct_survivor_generation.md; tests/test_orbit_canonical.py]
3. (3,8) gate: first rank with interference vertices; the honest stress test.
   Flat lattice already enumerated,
   [1, 56, 1540, 21420, 147630, 467082, 565208, 166420], 166,420 admissible
   hyperplanes against C(56,7) = 2.3e8 for brute force.
4. (3,9)/(3,10) blind reproduction against the census's own tables.
5. Full (3,11): candidates, Ressayre screen, relative H-system, then the
   exhaustiveness proof per the claim ladder in
   `docs/fixed_n3_generator_methodology.md`. No level of the ladder implies
   the next without its named certificate.
- New-rank vertices feed the fiber out-of-sample tests using the two
  targeting rules in `docs/GENERATOR_PLAN.md` (two-block for the cancellation
  regime, multi-block unequal for the gap sample; never collapse them).

### T3. Physics: complex engine benchmark. DONE 2026-08-05.

The standing flag is cleared. The asymmetric spinless t-V ring with exact
boundary phase `(4+3i)/5` gives 10 of 12 usable scan points, five with genuine
carrier flux. Independently reverified here from scratch in occupation space:
ground energy -0.1363514585, gap 1.434787, min occupation gap
0.0016711196969708686, all matching the artifact. The run also corrected a real
bug invisible to the zero-flux benchmark: under `gamma[p,q] = <a_p^dag a_q>`
the passive CI transform is the third exterior power of `U.T`, giving 2.094e-16
off-diagonal, where `U.conj().T` gives 0.1694.
[scripts/verify_complex_flux_physical_benchmark_standalone.py;
results/data/conditional_2rdm_program.json]

Still open, and not a repeat of the above: extension to molecular magnetic
fields, complex Bloch orbitals, or spin-orbit terms.

### T4. Backlog (eligible when T0 is done and nothing above is in flight)

- The 129 open gauge lower-bound rows: a NEW invariant family, not more of
  the exterior/Koszul contraction family.
- Exact recognition of the remaining new attainers from the endpoint set of 8.
- m >= 4 conditional-body posture: chase exactness with different structure or
  freeze outer-bounds-with-stated-status; record the decision as a doc either
  way (the rank-two extreme point closed the elliptope route).
- Decide the two unregenerable artifacts named in the cleanup section.

## Do-not-do list

- No SOS cand-44 session (superseded by the level-5 certificates).
- No census-wide additions to paper 1 pre-submission.
- No hand edits to `results/report/anonymized/` or to any emitted block.
- No count quoted anywhere before ledger regeneration and green CI.
- Nothing that a later HANDOFF retracted comes back without new evidence.

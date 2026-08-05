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

1. (3,7) exhaustiveness via the (3,6) sandwich template (exact plethysm inner
   points plus H-to-V closure). Until then (3,7) stays labeled relative.
2. (3,8) gate: first rank with interference vertices; the honest stress test.
3. (3,9)/(3,10) blind reproduction against the census's own tables.
4. Full (3,11): candidates, Ressayre screen, relative H-system, then the
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

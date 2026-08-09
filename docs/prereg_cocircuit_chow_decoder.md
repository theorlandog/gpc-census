# PRE-REGISTRATION: direct cocircuit discovery and constructive signed-Chow decoding (2026-08-09)

Written BEFORE the in-repo measuring script was run. Standing rules:
`docs/prereg_template.md`. Certifying artifacts:
`scripts/cocircuit_chow_decoder.py`,
`scripts/verify_cocircuit_chow_decoder_standalone.py`,
`results/data/cocircuit_chow_decoder.json`.

## Disclosure, because this is not a blind pre-registration

An external bundle (`gpc_cocircuit_chow_experiment`, a dependency-free C++
cocircuit prototype plus a Python profile decoder) was replayed in full before
this file was written: its binary was rebuilt from source here, its four
cocircuit runs and both decoder verifications were rerun, and its output CSVs
were byte-compared against the copies it shipped. So the predictions below are
INFORMED, not blind, and P1 to P4 are best read as regression targets for an
independent reimplementation rather than as discoveries. R4 is satisfied in the
strongest possible way: the numbers were checked against the repository corpus
first, and the two that reproduce sealed repository populations are named as
such.

What is genuinely untested at the time of writing is everything about the
IN-REPO code: no line of `gpc_census.generation.cocircuit` or
`gpc_census.generation.chow_decoder` existed when the external replay ran, they
share no code with the bundle, and the cost comparison against the repository's
own two enumerators has never been measured on one machine.

**R1 is only partly met, and the shortfall is stated rather than papered
over.** The predictions were written before the in-repo measuring script was
run, but this file and the script were committed together AFTER that run, in
one commit, so there is no earlier hash proving the predictions unedited. The
sequence was: replay the bundle, write the module, write this file, run,
score, commit. Read P1 to P4 accordingly: they are regression targets whose
values were already known, and the disclosure above is what carries their
weight, not the commit order.

## What is being tested

Three claims, in decreasing order of how much they would change the generator.

**C1 (cocircuit equivalence).** Admissible-hyperplane enumeration for the
exterior-weight configuration is the same problem as enumerating cocircuit zero
sets, that is maximal coplanar subsets of linear rank `d-1`. The operational
content is that a search may go straight to codimension one without
constructing any lower-rank flat.

**T4 (block symmetry).** If `c_+(tau)_i == c_+(tau)_j` then the transposition
of orbitals `i` and `j` maps `B_+(tau)` to itself, so `B_+` is a union of
occupancy profiles over the equal-degree blocks of `c_+`.

**T5 (block-constant separator).** Averaging a separator over the product of
symmetric groups on those blocks preserves the strict-positive family and
yields a separator constant on each block, whose levels are ordered by `c_+`.
Consequently the selected occupancy profiles form an upper ideal in prefix-sum
dominance.

Together T4 and T5 are meant to close the caveat `docs/signed_chow_projection.md`
records as open: "Recovering `B_+` from a bare `c_+` is the Chow-parameters
inversion problem and is not solved here."

## Holdout and its provenance

The populations are the complete admissible-hyperplane sets at `(3,5)` through
`(3,8)` and the complete nonstructural trace-survivor populations at `(3,6)`,
`(3,7)` and `(3,8)`. These are exhaustive enumerations, not samples, so R2's
transport question does not arise in its usual form. The relevant
non-independence is different and is stated here instead:

- The `(3,6)`, `(3,7)` and `(3,8)` hyperplane counts (362, 5,341, 166,420) and
  the orbit counts (8, 19, 56) are SEALED REPOSITORY POPULATIONS, already
  derived three ways. A new backend reproducing them is a regression test with
  an independent sample size of ZERO. It is still worth running, because the
  algorithm is new, but it must not be reported as new information about the
  polytopes.
- `(3,5)` (30 hyperplanes) is not separately sealed in the repository and is
  scored as a consistency point only.
- The decoder is tested on every nonstructural trace survivor, so the decoder
  rows ARE new information: nothing in the repository has previously inverted a
  shadow. Independent scope equals scope for the decoder predictions.

## Power caveat, stated before scoring

A PASS on P1 and P2 establishes that a second, structurally unrelated algorithm
reaches the same populations. It establishes NOTHING about facetness,
redundancy, birationality, or completeness of any system, and it does not make
the cocircuit route the production enumerator: the search here is unsymmetric
and therefore emits every hyperplane rather than every orbit, which is the
wrong scaling for rank 9 and above where `orbit_canonical` already works.

A PASS on P5 and P6 closes the decoding caveat as an EXISTENCE-AND-PRACTICE
statement on the populations measured. It is not a complexity theorem, and P7
is written specifically so that the absence of a complexity bound is scored
rather than glossed.

## Predictions

P1. The in-repo cocircuit enumerator returns 30, 362, 5,341 and 166,420
admissible hyperplanes at `d = 5, 6, 7, 8`. Expectation: PASS, because the
external prototype did and C1 says the objects are the same.

P2. At `d = 5, 6, 7` the cocircuit output agrees with
`enumerate_admissible_hyperplanes_by_flats` PAIR FOR PAIR in the `(h, z)`
convention, not merely in count. Expectation: PASS. This is the prediction that
would catch a convention error, which `docs/orbit_canonical_flats.md` records
as the failure mode that once produced 15 orbits instead of 8.

P3. Reverse-search node counts are 154, 2,775, 63,136 and 2,042,570.
Expectation: PASS, and this is the sharpest test of the reimplementation: node
counts depend on both prunings firing at exactly the same nodes, so an
agreement here is much stronger evidence of algorithmic identity than an
agreement on outputs.

P4. Trace survivors are 10, 64, 549 and 6,607, with 8, 62, 547 and 6,605
nonstructural. `S_d` unoriented orbit counts computed from primitive-normal
coordinate multisets are 8, 19 and 56 at `d = 6, 7, 8`. Expectation: PASS.

P5. The decoder reconstructs `B_+`, `B_-` and the primitive `tau` on every
nonstructural trace survivor at `(3,6)`, `(3,7)` and `(3,8)` with zero
failures. Expectation: PASS.

P6. Decoding is by a small finite search, not brute force: the maximum number
of branch states per signed side is below 100 at `(3,8)`, against a slice of 56
determinants and hence `2^56` naive subsets. Expectation: PASS.

P7. Search size grows with rank on random threshold shadows at `d = 8, 10, 12`
and `20`, and the growth is fast enough that no polynomial bound may be
claimed from the data. Expectation: PASS as a NEGATIVE result. If the maximum
state count instead stayed flat, that would be evidence toward a bound and
would be recorded as such.

P8. The direct cocircuit enumerator is SLOWER at `(3,8)` than
`orbit_canonical`'s augmentation enumerator on the same machine. Expectation:
PASS, meaning the new route does not replace the existing one. A FAIL here
would be a genuine finding and would change the recommendation in
`AGENT_PLAN.md` T2.

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED per prediction, decided by
`results/data/cocircuit_chow_decoder.json` alone. Each verdict reports scope
and independent scope. No prediction is edited after the run. Timings are
reported per machine and never averaged, per the sync invariant.

## SCORED (2026-08-09, after running; predictions above are unedited)

| id | verdict | scope | independent scope | measured |
| --- | --- | ---: | ---: | --- |
| P1 | **PASS** | 4 systems | 0 | 30, 362, 5,341, 166,420 admissible hyperplanes |
| P2 | **PASS** | 4 systems | 0 | pair-for-pair `(h,z)` agreement at `d = 5, 6, 7` and also at `d = 8` |
| P3 | **PASS** | 4 systems | 0 | 154, 2,775, 63,136, 2,042,570 reverse-search nodes |
| P4 | **PASS** | 4 systems | 0 | survivors 10, 64, 549, 6,607; nonstructural 8, 62, 547, 6,605; orbits 8, 19, 56 |
| P5 | **PASS** | 7,214 rows | 7,214 | 0 decode failures, 0 primitive-`tau` failures |
| P6 | **PASS** | 6,605 rows | 6,605 | 37 branch states at worst per signed side at `(3,8)`, mean 6.54 |
| P7 | **PASS**, negative | 4 systems | 4 | max states 53, 137, 253, 1,717 against 56, 120, 175, 442 profiles |
| P8 | **PASS** | 1 system | 0 | cocircuit 278.7 s, flat lattice 140.3 s, orbit augmentation 10.1 s |

Eight for eight, and the two that matter most for what happens next are the
ones scored as negatives.

**P2 exceeded its own statement.** The pair-for-pair comparison was predicted
for `d = 5, 6, 7`, because a rank-8 comparison needs the materialising
enumerator to finish. It did finish, in 140.3 s, so the comparison also ran at
`d = 8` on all 166,420 pairs. The trace test was additionally cross-computed in
the repository's centered convention at every rank, agreeing everywhere.

**P7 is the one to read carefully.** State counts grow 53, 137, 253, 1,717
while profile counts grow 56, 120, 175, 442, so the ratio moves from 0.95 to
3.89 across `d = 8` to `d = 20`. Four points do not fix an asymptotic, which is
exactly why this is scored as a negative result: nothing here licenses a
polynomial bound.

**P8 carries the architectural consequence.** The handoff that prompted this
work recommended making the direct route the preferred pipeline. It is 27x
slower than the orbit-canonical enumerator at `(3,8)` and holds 166,420 objects
where that enumerator holds 314, and the object gap widens by three more orders
of magnitude by rank 10. The recommendation is declined in
`docs/cocircuit_chow_decoder.md` and in `AGENT_PLAN.md` T2d.

**One measurement error was made and caught before scoring**, recorded per R5.
The first run of the P5 gate compared the decoded normal against a
sign-normalized representative of the true normal, which discards the
orientation the decoder is supposed to recover, and reported 31, 365 and 4,956
"failures" at ranks 6, 7 and 8. Those were an artifact of the comparison, not
of the decoder: with oriented equality the count is 0 at every rank. The gate
now compares oriented tuples, which is the stronger check, and the wrong first
reading is recorded here rather than deleted.

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
(A second one, in round 2, is recorded there.)
The first run of the P5 gate compared the decoded normal against a
sign-normalized representative of the true normal, which discards the
orientation the decoder is supposed to recover, and reported 31, 365 and 4,956
"failures" at ranks 6, 7 and 8. Those were an artifact of the comparison, not
of the decoder: with oriented equality the count is 0 at every rank. The gate
now compares oriented tuples, which is the stronger check, and the wrong first
reading is recorded here rather than deleted.

---

# ROUND 2: arbitrary particle number (2026-08-09)

Round 1 above is closed and unedited. This is a separate round with its own
predictions, written after round 1 was committed and before
`scripts/cocircuit_chow_higher_n.py` was run. Certifying artifacts:
`scripts/cocircuit_chow_higher_n.py`,
`results/data/cocircuit_chow_higher_n.json`,
`tests/test_cocircuit_chow_higher_n.py`.

## Disclosure

A second external bundle arrived claiming the whole route generalizes beyond
`N = 3`, with its own C++ backend and decoder. As in round 1 the predictions
are INFORMED rather than blind: the bundle's headline numbers were read before
this was written. Unlike round 1, no part of the bundle was executed here. The
in-repo modules were already written for general `n`, so what is under test is
whether that generality is real or merely typed, and one prediction below is
written specifically to contradict the bundle.

R1 shortfall is the same as round 1 and for the same reason: this section and
the script are committed together after the run.

## What is being tested

**G1 (generality).** C1, T1 to T5 hold at arbitrary `N`, since C1 needs only
`1 . chi_I = N != 0` and T1 to T5 are fixed-weight statements. Operationally:
the existing modules, written with an `n` parameter, are correct at `n != 3`.

**G2 (particle-hole shadow).** An `n`-uniform family `F` of size `m` with
shadow `c` has complement family `F^c` in the `(d-n)` layer with shadow
`m - c`, so decoding either layer is exact.

## Holdout and its provenance

`(4,7)` and `(5,8)` are particle-hole duals of `(3,7)` and `(3,8)`, so their
populations are DETERMINED by round 1 and carry independent sample size zero.
The published higher-N representatives at `(4,8)`, `(4,9)`, `(4,10)` and
`(5,10)` are census rows this work has never touched, so they are independent
of everything in round 1, and they are the only rows here that are. Three of
those four systems store one representative per Levi signature class rather
than every facet row; that is stated here, before scoring, so a PASS is not
read as a complete-population result.

## Power caveat, stated before scoring

A PASS on G1 says the code is `N`-general on the systems tested. It does not
produce a single new higher-N population, and it must not be read as progress
toward a new complete system at any `N`.

## Predictions

P9. `(4,7)` and `(5,8)` reproduce 5,341 and 166,420 admissible hyperplanes, 19
and 56 orbits, 549 and 6,607 trace survivors, agreeing with the flat lattice
pair for pair. Expectation: PASS, forced by particle-hole duality if the code
is correct.

P10. Every published higher-N representative passes the whole chain, including
ORIENTED reconstruction of the primitive normal. Expectation: PASS. This is the
prediction with real content, because these rows come from a different campaign
and were never used to develop the decoder.

P11. The particle-hole reduction leaves the profile VARIABLE COUNT unchanged on
every row. Expectation: PASS, and it CONTRADICTS the handoff, which describes
decoding in the smaller layer as something that can dramatically reduce the
profile search. `k -> m - k` is a bijection between the two layers' profiles
over the same blocks, so the count cannot change. If instead the counts differ
anywhere, my reading of the bijection is wrong and the handoff is right.

P12. The particle-hole route nonetheless reduces total branch states, because
the coefficients and target differ. Expectation: PASS as a constant-factor
effect only, with some rows getting dearer rather than cheaper.

P13. Exact trace survivors at `(5,10)` and `(6,12)` decode with zero failures,
and the largest profile problem stays far below the slice size. Expectation:
PASS.

P14. `(4,8)` does NOT complete inside a 2,000,000-node reverse search.
Expectation: PASS, meaning the wall arrives earlier at half filling than the
`N=3` ladder suggests. A FAIL here would be good news and would reopen the
direct route.

## Scoring rules

As round 1. Verdicts are decided by
`results/data/cocircuit_chow_higher_n.json` alone.

## SCORED (2026-08-09, after running; predictions above are unedited)

| id | verdict | scope | independent scope | measured |
| --- | --- | ---: | ---: | --- |
| P9 | **PASS** | 2 systems | 0 | `(4,7)` 5,341/19/549 and `(5,8)` 166,420/56/6,607, both pair-for-pair |
| P10 | **PASS** | 61 rows | 61 | 61/61 on span, trace, both sign families, oriented `tau` |
| P11 | **PASS**, contradicting the handoff | 104 shadows | 104 | variable counts differ on 0 rows at `(5,8)`, `(6,9)`, `(7,10)` |
| P12 | **PASS**, constant factor | 104 shadows | 104 | states 317 to 183, 499 to 269, 1,578 to 354; 51 rows cheaper, 4 dearer |
| P13 | **PASS** | 30 rows | 30 | 0 failures; 38 profile variables against a 252 slice at `(5,10)`, 40 against 924 at `(6,12)` |
| P14 | **PASS** | 1 system | 1 | `(4,8)` did not complete in 2,000,000 nodes (304.5 s), no partial count recorded |

Six for six, and again the informative ones are the ones that push back. P11
was written to be falsifiable against the handoff's own description and the
handoff loses: the reduction is exact and is not a size reduction. P14 says the
direct route's wall arrives two ranks earlier at half filling than the `N=3`
ladder suggests.

**An unpredicted finding, which is the most useful thing in round 2.** P10
passed, but the profile counts inside it are the story: the worst `(4,10)`
representative needed 210 profile variables against a slice of `C(10,4) = 210`,
and the worst `(5,10)` one needed 252 against `C(10,5) = 252`. T4 compresses a
shadow exactly as much as that shadow is degenerate, and on those two published
rows it compresses nothing at all. No prediction covered this because the
`N=3` populations are degenerate enough that the question never arose. It is
recorded in `docs/cocircuit_chow_decoder.md` and in the "what is not claimed"
list rather than left in the artifact for someone to notice.

**A second measurement error, caught before scoring**, recorded per R5. The
P13 harness first reported 8 failures at `(5,10)` and 4 at `(6,12)`. Every one
was a sampled normal whose content exceeded one being compared against the
decoder's primitive output; the decoder had recovered the correct ray with the
correct orientation on all of them. The harness now divides the content out
before comparing, which is what the comparison always meant, and the wrong
first reading is recorded here rather than deleted. This is the same class of
mistake as round 1's orientation slip: both times the decoder was right and the
comparison was wrong.

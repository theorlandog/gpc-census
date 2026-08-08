# The Bulois-Denis-Ressayre constructor comparison

Pre-registration: `docs/prereg_bdr_constructor_comparison.md`, committed at
`2e0216e` before any measurement. Generating script:
`scripts/bdr_constructor_comparison.py`. Artifact:
`results/data/bdr_constructor_comparison.json`. Guard test:
`tests/test_bdr_constructor_comparison.py`.

> ## ❗ CORRECTION (2026-08-08), read before anything below
>
> This document claimed the pinned external implementation "could not be
> fetched or run in this environment", and listed four blockers. **That
> conclusion is false, and it was never tested.** The code is mirrored at
> `github.com/ea-icj/moment_cone`, the session's repository-attachment path
> serves it, the clone lands on the exact pinned revision
> `7ab347d9a7837d68d92bde9fac74606913f395ca`, its licence is **MIT** (not
> `not_verified`), and the upstream `passagemath` extra installs from PyPI into
> a throwaway venv. The backend has been run. See
> `docs/bdr_linear_triangular_gate.md` for the recipe and the results.
>
> The **INCOMPLETE** label below is still correct, for a different reason: this
> benchmark measures wall clock and stage counts, and no BDR timing has been
> taken even now. Everything measured here remains valid. What was wrong is the
> stated reason, and the four blockers are restated accurately in the section
> below.

## Headline

🔬 The benchmark is **INCOMPLETE**: no BDR runtime, candidate count, or stage
figure appears anywhere in this document or its artifact, because this
benchmark never ran the external path.

🟦 What was measured instead settles the practical question without it. Under
the standing rule that every row entering a published system carries a locally
replayed exact certificate, the maximum speedup available to ANY external
candidate generator was measured at ranks 7 and 8, and it is much smaller than
the raw candidate-count reduction suggests:

<!-- not emitted: three derived ratios, all from one artifact, read below -->

| bound at `(3,8)` | speedup | what it assumes |
|---|---|---|
| perfect birationality prefilter alone | **1.82x** | local enumeration kept, every determinant-zero row removed for free |
| plus free external candidate generation | **6.27x** | non-circular: the generator must still emit at least the certified rows |
| oracle handed the final 31 facets | 16.38x | circular, unreachable at a new rank |

The middle row is the honest ceiling. It assumes the external generator is both
free and perfect, and it is still only 6.27x, because the irreducible work is
screening and reducing the 137 Ressayre-certified rows, which no external tool
can remove without also deciding redundancy.

❌ The one BDR step with no local counterpart, the birationality filter, is
worth **1.82x at rank 8 even when it is free and perfect**. That is the number
that decides the sprint.

## What was not done, precisely

The external object is pinned:

```text
paper       An Algorithm to compute the Kronecker cone and other moment cones
authors     Michael Bulois (ICJ, AGL), Roland Denis (INSMI-CNRS, ICJ),
            Nicolas Ressayre (ICJ, AGL)
preprint    arXiv:2505.08812 (math.AG); HAL hal-05044759 version v3
code        ea-icj/moment_cone at plmlab.math.cnrs.fr; https://ea-icj.github.io/
mirror      github.com/ea-icj/moment_cone  (this is the one that works here)
revision    7ab347d9a7837d68d92bde9fac74606913f395ca
runtime     Python-Sage (SageMath), or the passagemath fork
license     MIT, read from LICENSE at the pinned revision
case        fermionic cone, GL_d(C) on wedge^r C^d, here r = 3
```

Four blockers were listed. Restated accurately, after testing each:

1. `plmlab.math.cnrs.fr` is refused by this session's egress proxy with a 403
   on CONNECT. **True, and not load bearing:** the code is mirrored on GitHub.
2. `ea-icj.github.io`, `arxiv.org` and `cnrs.hal.science` are refused by the
   same policy, so the paper text could not be read here. **True, and not load
   bearing:** none of them is needed to run the code.
3. `moment-cone` is not published on PyPI. **True (404), and not load
   bearing:** `pip install .` from the clone works.
4. SageMath must never become a package runtime dependency. **True, and
   unchanged:** the upstream ships a `passagemath` extra that installs into a
   separate throwaway venv, touching neither `pyproject.toml` nor `uv.lock`.

The inference drawn from those four, that the backend could not be fetched or
run, does not follow and is false here. Consequence for the correspondence
table below: at the time it was written the BDR column was a
**secondary source**, drawn from the call-graph audit already in
`docs/fixed_n3_generator_methodology.md` (written in an earlier session that
did have the code installed) plus the paper's published description. It is
labelled as such per cell and must not be read as a fresh source read. Those
cells are now checkable against the source, and one of them has been checked:
the `linear_triangular` filter, in `docs/bdr_linear_triangular_gate.md`. The
rest stay marked secondary until someone reads them.

## Task 1: algorithm correspondence

| Concept | BDR implementation | gpc-census implementation | Same theorem? |
|---|---|---|---|
| candidate generation | `TauCandidatesStep` calls `tau.find_1PS`, a recursive enumeration of one-parameter subgroups [secondary] | orbit-canonical closed-flat enumeration, one exact `(h,z)` per top-level `S_d` orbit (`generation.orbit_native`) | **No.** Different constructions of the same target set. Ours enumerates flats of the augmented exterior weights; theirs enumerates 1-PS. Neither is derived from the other here. |
| admissibility | `SubModuleConditionStep` filters by `Tau.is_sub_module` [secondary] | affine rank `d-1` of the on-hyperplane exterior weights, by a modular nonzero minor with exact rational fallback | **No.** Both exact, different criteria. The local one is Ressayre admissibility directly. |
| trace condition | `InequalityCandidatesStep` calls `List_Inv_Ws_Mod`, a deterministic partition recursion over compatible inversion sets [secondary] | `inv(sigma h) == B` with `B` an orbit invariant; survivors generated directly as arrangements of the multiset `h` with `B` inversions | **Yes, and this is the striking one.** Both are inversion-set enumerations over the same combinatorial object. |
| determinant / birationality | `PiDominancy` and birationality filters, characterised by the in-repo audit as explicitly probabilistic in the pinned revision [secondary] | exact Bareiss determinant at a deterministic integer point, Hall matching pre-test, exact symbolic fallback; failure is `unresolved`, never a rejection | **No, and not adoptable as-is.** Theirs decides; ours must be fail-closed and exact. |
| Weyl chamber orientation | dominant chamber, `H.lambda >= z` form | ordered slice `lambda_1 >= ... >= lambda_d`, `lambda_1 <= 1`, `lambda_d >= 0`, as a fundamental domain; negative type-A roots `e_j - e_i`, `i < j` | **Equivalent** under the convention map below. |
| redundancy / minimality | produces the minimal list of inequalities for the representation | exact Farkas multipliers plus rational facet witnesses on the ordered slice, replayed in removal order | **Same claim, different evidence.** Only the local one ships a replayable certificate per removal and per retained row. |
| completeness statement | minimal list for one fixed representation | fixed `(N,d)` completeness from `S_d` orbit coverage plus closed-form survivor count conservation | **Both fixed-representation.** Neither is uniform in `d`. See Task 6. |

### Convention map, exercised not asserted

The bridge is `src/gpc_census/generation/bdr.py`, and
`tests/test_bdr_constructor_comparison.py` runs the round trip on every
published row of `(3,6)`, `(3,7)` and `(3,8)`.

```text
BDR homogeneous     tau . lambda <= 0, tau in Z^d
local affine        a . lambda <= b on the slice sum(lambda) = N
forward             tau_i = N * a_i - b, divided by its content
inverse             constraint_from_tau, deterministic nonnegative
                    representative with least coefficient zero
orientation         positive rescaling only; a sign change is a DIFFERENT
                    row and is never used to make a match succeed
Ressayre H form     S = sum(tau), H_i = S - d * tau_i, z = 3 * S
root orientation    negative type-A roots e_j - e_i with i < j, H_j - H_i < 0
particle number     N = 3, trace slice sum(lambda) = 3
structural row      (0,...,0,-1) is lambda_d >= 0, separated from GPC rows
primitive gcd/sign  divide by the positive gcd; zero content is an error
```

One convention trap was found and fixed while building the ledger, and it is
worth recording because it is the kind of thing that silently produces a
"no match" result. The orbit enumerator stores **one** orientation per
unoriented orbit, and which one is not a documented convention. Choosing a
canonical key by an arbitrary rule (the minimum of the two sorted `(h,z)`
keys) failed to match the Borland-Dennis row at rank 6, whose stored orbit
carries the other orientation. The ledger now tests **both** keys for
membership. With that fix every published row matches a local orbit at all
three ranks.

## Task 2: exact set comparison on the overlapping systems

`(3,6)` is present but **screening only**: `compose_candidate_system` guards
ranks 7 through 13, and rank 6 is covered by the separate exhaustive reference
gate in `docs/stage1_ressayre_3_6.md`. It is excluded from the ceiling
explicitly, with the reason stored in the artifact, rather than dropped.

`(3,9)` and `(4,8)` were excluded in the pre-registration, before the run.

Full local construction, reproducing the counts pinned in
`docs/orbit_native_constructor.md` exactly:

| stage | `(3,7)` | `(3,8)` |
|---|---|---|
| unoriented orbits | 19 | 56 |
| admissible hyperplanes | 5,341 | 166,420 |
| oriented rows total | 10,682 | 332,840 |
| never materialised (trace) | 10,133 | 326,233 |
| survivors generated | 549 | 6,607 |
| survivors predicted, closed form | 549 | 6,607 |
| structural rows | 1 | 1 |
| screened candidates | 548 | 6,606 |
| certified | 49 | 137 |
| exact determinant zeros | 499 | 6,469 |
| unresolved | 0 | 0 |
| redundant | 45 | 106 |
| retained | 4 | 31 |

The oracle path, supplied only the published minimal rows:

| stage | `(3,6)` | `(3,7)` | `(3,8)` |
|---|---|---|---|
| supplied rows | 1 | 4 | 31 |
| certified | 1 | 4 | 31 |
| rejected at any stage | 0 | 0 | 0 |
| removed by reduction | n/a | 0 | 0 |
| retained | n/a | 4 | 31 |

The stages do not align one-to-one with BDR's, and no attempt is made here to
force them into a shared count. Three of the seven rows in the correspondence
table are marked "No", so a stage-by-stage candidate-count table against BDR
would be a fabricated comparison even if the code had run.

## Task 3: independent recertification, no verdict by citation

Every published row of all three systems was converted to primitive `tau` and
put through the local exact verifier. The artifact stores, per row: the
external row, the converted `tau`, the local status, the local resolution
method, the local certificate SHA256, and the matched local orbit.

🟦 **CONFIRMED, 36 rows across 3 systems.** Every one screens `certified` by
deterministic evaluation, every certificate replays and binds back to its
oriented `tau`, and every row matches a local orbit. No row is
`trace_rejected`, `inadmissible`, `symbolic_zero_rejected`, or
`evaluation_unresolved`.

This is the prediction P1 evidence, and it is also the answer to research
question 5: yes, external output converts into the local `tau` convention and
is independently certifiable here, with a per-row certificate hash.

## Task 4: work avoided, by stage

Wall clock, single worker, one machine (`Linux-6.18.5-fc-v20-x86_64`, CPython
3.14.0rc2). Repeat runs agreed to within 2 percent on the stage shares. Not
averaged with any other machine, per the standing benchmark rule.

`(3,8)`, the decisive system:

| path | enumeration | screening | composition | total | speedup |
|---|---|---|---|---|---|
| full local | 6.25 s | 7.64 s | 2.11 s | 16.00 s | 1.00x |
| perfect birationality prefilter | 6.25 s | 0.68 s | 1.87 s | 8.80 s | 1.82x |
| prefilter and free enumeration | 0 | 0.68 s | 1.87 s | 2.55 s | 6.27x |
| oracle handed the answer | 0 | 0.15 s | 0.83 s | 0.98 s | 16.38x |

Does the external method avoid the named work?

- **Lower-flat enumeration.** Unknown, and the honest answer is 🔬 OPEN: BDR
  enumerates 1-PS, not flats, so its cost is not comparable without running it.
  The bound above prices it at zero and still lands at 6.27x.
- **Top-level orbit canonicalization.** No evidence either way [secondary
  source silent]. 🔬 OPEN.
- **Trace-rejected orientations.** ❌ Already removed locally. The local
  constructor never materialises 326,233 of 332,840 rows at rank 8. There is no
  work here for BDR to take.
- **Hall-structural zeros.** Mostly already removed. Of the 6,469 rows reaching
  the symbolic fallback at rank 8, 6,414 are decided by Hall matching and only
  55 reach the exponential dynamic program
  [`docs/orbit_native_constructor.md`]. A birationality prefilter would be
  displacing a matching test that already runs fast.
- **Redundant valid rows.** ❌ Cannot be avoided by any external generator that
  is not circular. Deciding that 106 of the 137 certified rows are redundant
  IS the reduction step.
- **Determinant algebra.** Partly, and this is BDR's real offer: 1.82x at rank
  8 when free and perfect.

The single most useful number: the local pipeline already avoids **98.0
percent** of the oriented rows at rank 8 before any screening happens. The
candidate-count reduction an external generator can still offer is bounded by
the 6,606 rows that remain, and 6,469 of those are determinant zeros.

## Task 5: the hybrid, and why it is not integrable

The oracle path IS the hybrid the brief names, "BDR minimal rows into local
certification", with a perfect oracle standing in for BDR. It passes the gate
the brief demands:

🟦 **CONFIRMED.** At ranks 7 and 8 it reproduces the local retained system
exactly: oriented `tau` sets equal in both directions, `generated_only` and
`published_only` both empty, `known_system_regression` passed, and the retained
`tau` set identical to the full path's. Zero rows removed by reduction, zero
unresolved.

❌ **It must not be integrated, and not because it fails a gate.** It is
circular. It is handed the final system, so it can only run where the answer is
already known. Its value is exclusively as an upper bound. Reporting its 16.38x
as an achievable speedup would be exactly the kind of claim this repository
exists to prevent.

The non-circular hybrid is the prefilter path, which does reproduce the retained
system exactly at both ranks, and is worth 1.82x at rank 8.

## Task 6: completeness scope audit

What a fixed-representation algorithm proves, separated as the brief requires:

| statement | BDR, for a fixed representation | gpc-census, at ranks 7 and 8 |
|---|---|---|
| candidate exhaustiveness | claimed by the algorithm's design; **not proved in this repository**, and the in-repo audit names three unclosed coverage gaps in `find_1PS`, `is_sub_module` and `List_Inv_Ws_Mod` | 🟦 proved, by `S_d` orbit coverage plus closed-form survivor count conservation |
| validity of retained rows | asserted by the method | 🟦 proved, exact Ressayre determinant certificate per row, replayed |
| minimality of a supplied system | the stated output | 🟦 proved, exact Farkas removals replayed in order |
| facetness | not separated in the secondary source | 🟦 proved, rational facet witness per retained row on the ordered slice |
| full moment-cone completeness | for the fixed representation | 🟦 proved at `(3,7)` and `(3,8)` |
| uniformity in `d` | **not claimed, and does not follow** | **not claimed, and does not follow** |

The last row is the load-bearing one for the wider agent pack. A fixed-`(N,d)`
completeness theorem says nothing about stable compression across `d`. Neither
algorithm supplies a rank-independent statement, and the measured trend below
is evidence, not a theorem, about what would happen if one did.

## The trend that decides the next sprint

🟦 **CONFIRMED within a two-point finite scope.** The oracle ceiling **falls**
from rank 7 to rank 8, 21.41x to 16.38x, while the raw row-reduction ratio
**rises**, 548/4 = 137 to 6,606/31 = 213.

The mechanism is visible in the artifact and is not a resource accident. The
oracle path's own cost is dominated by composition, 0.83 s of its 0.98 s at
rank 8, and composition scales with the retained facet count, which grows fast
across the census: 1, 4, 31, 52, 93 at ranks 6 through 10. So the floor an
external generator cannot remove is growing faster than the work it could take
away.

Two points are a trend, not a theorem, and this is explicitly **not** a claim
about all `d`. `tests/test_bdr_constructor_comparison.py` pins the measured
direction so a rerun cannot lose it silently.

## Scored pre-registration

| id | prediction | verdict | scope | independent scope |
|---|---|---|---|---|
| P1 | every published row screens locally as `certified` | **PASS** | 36 rows | 3 systems |
| P2 | oracle path reproduces the published system exactly | **PASS** | 2 systems | 2 |
| P3 | reduction removes zero rows from each published system | **PASS** | 2 systems | 2 |
| P4 | oracle total for `(3,8)` under 5 s | **PASS** (0.98 s) | 1 | 1 |
| P5 | enumeration is over half the `(3,8)` wall clock | **FAIL** (0.39) | 1 | 1 |
| P6 | `(3,8)` oracle ceiling at least 5 | **PASS** (16.38x) | 1 | 1 |

### P5 FAILED, and it matters

The prediction came straight from `docs/orbit_native_constructor.md`, which
says the orbit-native path "still pays for the orbit ENUMERATION, which is now
the dominant stage".

On this machine that does not reproduce. At rank 8 the split is enumeration
6.25 s, screening 7.64 s, composition 2.11 s: **screening dominates**, and
enumeration is 39 percent, stable to within 1.4 percentage points across three
runs.

This is not asserted as a refutation of that document. Its measurement was
taken on a different machine, where the same construction took 66.0 s against
16.00 s here, a 4x difference, and multi-machine benchmark results are reported
per machine and never averaged. The correct status is 🔬 **OPEN: the
"enumeration dominates" claim does not reproduce on this machine and needs
re-measuring on the original one before it is quoted again.**

The failure strengthens rather than weakens the main conclusion. If screening,
not enumeration, is the dominant stage, then an external tool that replaces
candidate generation is attacking the smaller half of the work.

## What this does and does not establish

**Establishes.**

- The exact convention map between the two representations, exercised on 36
  published rows with round trips that never use a sign flip.
- That every published row of `(3,6)`, `(3,7)` and `(3,8)` recertifies locally
  with a bound certificate, so external output is usable as candidates here.
- A measured upper bound on any external candidate generator under the local
  certification rule: 6.27x at rank 8 non-circular, 1.82x for the birationality
  step alone.
- That the trace stage, where the two algorithms genuinely agree
  mathematically, is already the local pipeline's fastest stage and offers BDR
  nothing to remove.

**Does not establish.**

- Anything at all about BDR's actual runtime, candidate counts, or memory. This
  benchmark did not run the backend. 🔬 OPEN, and now cheap to close: the
  backend has since been fetched, installed and run for
  `docs/bdr_linear_triangular_gate.md`, so the missing piece is a timing
  harness, not access.
- That BDR's candidate enumeration is or is not cheaper than the local
  closed-flat orbit enumeration. Priced at zero in the bound, unmeasured in
  fact.
- That the pinned revision's birationality filter is or is not probabilistic.
  That characterisation is from the in-repo secondary audit and needs a source
  read to confirm.
- Any statement uniform in `d`. Two ranks are two ranks.

## Recommendation

**Do not integrate, and do not spend another sprint on this route at ranks 7
and 8.** Even a free and perfect version of BDR's one distinctive step buys
1.82x at rank 8, and the full non-circular bound is 6.27x with everything
external priced at zero. That is the load-bearing reason and it is a
measurement, not an availability complaint.

The original recommendation also leaned on "the external code cannot be pinned
and reproduced here, so criterion 5 fails before any mathematics". **That leg is
withdrawn:** criterion 5 is satisfiable, the code pins and reproduces here, and
`docs/bdr_linear_triangular_gate.md` does exactly that for one filter. The
recommendation is unchanged because the 1.82x and 6.27x bounds are unchanged.

The route is worth reopening under a checkable condition: if the composition and
screening floor measured here stops dominating, or if the rank-11 bottleneck is
shown to sit in candidate generation rather than in screening and reduction. On
this machine at rank 8 it does not.

## Smallest next experiment justified by this result

Measure the same three-path split at `(3,9)`, where the retained count jumps to
52. The single question is whether the oracle ceiling continues to fall. A third
point turning 21.41, 16.38 into a monotone sequence would make the "external
candidate generation cannot reach the rank-11 bottleneck" statement a
finite-scope confirmed finding rather than a two-point trend, and it costs one
run of an already-written script.

## Reproducing

```sh
uv run python scripts/bdr_constructor_comparison.py
uv run pytest tests/test_bdr_constructor_comparison.py
```

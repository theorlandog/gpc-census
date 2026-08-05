# PRE-REGISTRATION: generator ladder gate 1, exhaustive (3,7)

Date: 2026-08-05. Standing rules in `docs/prereg_template.md`.

Committed BEFORE the measuring script is run, per R1. The commit that adds this
file also adds `scripts/generator_gate1_3_7.py`; the run happens after.

## What is being attempted

Promote `(3,7)` from RELATIVE to EXACT on the claim ladder in
`docs/fixed_n3_generator_methodology.md`. Level 4 (known-system match) is
already held: the Stage-1 manifest's four retained rows agree as an oriented
set with the stored Altunbulak-Klyachko table. Level 5 needs an
exhaustive-enumeration proof OR the independent exact inner/outer closure.

This gate attempts the closure route, not the enumeration route. Its three
obligations, all of which must be exact:

1. every H-row has a validity proof;
2. the H-polyhedron is vertex-enumerated without omission;
3. every enumerated H-vertex has an exact attainable-state certificate.

Validity gives `P subset H`. Exact attainability of every H-vertex gives
`H subset P` by convexity. Hence `P = H`.

## Frozen inputs

- H-rows: the four retained rows of `results/data/stage1_ressayre_3_7.json`,
  each carrying an exact Ressayre certificate, plus the structural rows of the
  ordered trace-three slice (`1 >= l1 >= ... >= l7 >= 0`, `sum l = 3`).
- Attainability: `results/data/states.jsonl`, records with `system == "(3,7)"`.
- Reference vertex list: `results/data/vertices/vertices_3_7.json`.

## Scored predictions

**G1. The H-system has exactly 10 vertices.**
PASS if exact rational vertex enumeration of the frozen H-rows plus structural
rows returns exactly 10 vertices. FAIL at any other count. A count above 10
means the four rows are not enough to cut the slice down to the polytope; a
count below 10 means the reference vertex list contains a non-vertex.

**G2. The enumerated vertex set equals the reference set exactly.**
PASS only on exact set equality of the 10 enumerated vertices with
`vertices_3_7.json`, compared as exact rationals. Any one-sided difference is
a FAIL and the difference is recorded.

**G3. Every enumerated vertex carries an exact attainable state.**
PASS if each of the 10 has a `states.jsonl` record with
`tierB in {EXACT, EXACT-CONSTR}` AND that state's 1-RDM spectrum is re-derived
here, independently of the record's own `integer_form`, and matches. FAIL if
any vertex lacks a state or if any re-derivation disagrees.

**G4. Every H-row is certified valid, replayed here.**
PASS if all four retained rows carry a Ressayre certificate whose trace
condition and nonzero integer determinant are re-checked in this run, not read
off the manifest. FAIL otherwise.

**G5. The closure is therefore exact: `P = H` for (3,7).**
PASS only if G1 through G4 all PASS. This is a composition prediction and has
no independent evidence of its own; it is listed so that a partial result
cannot be written up as a closure.

## Predictions that would falsify the approach, not just the run

**G6. The structural rows are not silently doing the work.**
Re-running the enumeration with the four Ressayre rows REMOVED must give a
strictly larger vertex set. PASS if the ordered slice alone has more than 10
vertices, so the four rows are load bearing. FAIL, and the gate is
uninformative about them, if the ordered slice alone already has exactly the
10 vertices.

**G7. Each of the four rows is individually load bearing.**
Dropping any single one of the four and re-enumerating must give a strictly
larger vertex set. PASS if all four drops increase the count. A row whose
removal changes nothing is redundant on this slice and is recorded as such,
which would not break the closure but would contradict the manifest's
irredundancy claim.

## Scope accounting, per R2

The ten `(3,7)` states are NOT independent out-of-sample observations of any
fiber law; several are transported or padded from `(3,6)`. This gate makes no
fiber-law claim, so R2's independence accounting does not gate it. It is
recorded here so the artifact cannot later be cited as ten independent
attainability observations. The artifact will carry, per vertex, whether its
state is direct or derived, and the count of each.

## What a PASS does NOT establish

- Nothing about `(3,8)` or higher. Gate 2 is a separate prereg.
- Not that the candidate ENUMERATION at `(3,7)` is exhaustive. The closure
  route bypasses that question rather than answering it; the manifest's
  `system_completeness` field stays as it is unless the enumeration route is
  separately completed.
- Not facetness of any row beyond what the existing facet witnesses show.

## Stopping rule

If G1 or G2 fails, stop and report the difference; do not adjust the H-system
to make the counts agree. If G3 fails for some vertex, that vertex is the
finding and the closure stays open.

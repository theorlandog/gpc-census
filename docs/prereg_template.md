# Pre-registration template and standing rules

Every new fiber law gets a pre-registered scored prediction before the
narrative is written, and a FAIL is a publishable finding (RESEARCH.md, "The
falsification protocol"). This file holds the shared rules so each
pre-registration does not restate them, and the template to copy.

## Standing rules

**R1. Predictions are committed before the measuring script is run.** Commit
the prediction file and the script together, then run. The commit hash is the
evidence that the predictions are unedited.

**R2. Derived measurements are not independent observations.** A measurement
obtained from a lower-rank or already-studied state by transport, padding,
frozen-core lift, face embedding, or particle-hole duality is NOT an
independent out-of-sample observation. Its structural fields are the donor's
fields, so it can confirm only the invariance of the derivation itself, which
must be tested as its own separate prediction.

Every scorecard MUST therefore report, per prediction:

- `scope`: total rows scored;
- `independent_scope`: rows whose root origin is not already in the corpus;
- `independent_distinct_states`: distinct root states behind those rows.

and the artifact MUST carry a per-measurement provenance record naming the
root origin. A scorecard that reports only row counts overstates its evidence
by the derivation factor, silently.

This rule was added after the (3,11)/(3,12) bracket holdout was found to have
34 of 38 rows derived from 17 rank-10 donors, with an effective independent
sample size of 2. See the PROVENANCE CORRECTION section of
`docs/prereg_bracket_3_11_3_12.md`.

**R3. An empty antecedent is UNINFORMATIVE, never PASS.** If a prediction says
"no member of class X appears" and no member of X could have appeared given
how the holdout was built, the verdict is UNINFORMATIVE. Record why the
antecedent was empty by construction.

**R4. A prediction already refuted by data in hand is not a discovery.**
Before scoring, check each prediction against the existing corpus. If a
census-wide fact already falsifies it, the verdict is ALREADY FALSIFIED, the
prediction was misformulated, and the write-up says so. Only genuinely new
information counts as a finding.

**R5. Corrections are visible, not silent.** When a scorecard is re-scored,
retain the original verdicts beside the corrected ones, in both the document
and the artifact. Never edit a scored verdict in place.

**R6. Numerical is labeled numerical.** Every measurement carries an
`evidence` field, and exact and numerical results are never totaled together.

## Template

```markdown
# PRE-REGISTRATION: <what is being tested> (<date>)

Written BEFORE the measuring script was run. Standing rules:
`docs/prereg_template.md`. Certifying artifact: `<script>`, `<artifact>`.

## What is being tested

<the law, stated precisely enough to fail>

## Holdout and its provenance

<what is being measured, and for each source: is it independent of the
existing corpus, or derived from it by transport / padding / embedding?
State the expected independent sample size HERE, before scoring.>

## Power caveat, stated before scoring

<what a PASS would and would not establish, given the above>

## Predictions

P1. <statement>. Expectation: <PASS/FAIL> because <rationale>.
...

## Scoring rules

PASS / FAIL / UNINFORMATIVE / ALREADY-FALSIFIED per prediction, decided by
the artifact alone. Each verdict reports scope, independent scope, and
independent distinct states. No prediction is edited after the run.

## SCORED (<date>, after running; predictions above are unedited)

<table with verdict, scope, independent scope, and the measured numbers for
every FAIL>
```

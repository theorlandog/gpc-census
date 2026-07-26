# PRE-REGISTRATION: the natural-orbital support ratio as a trichotomy diagnostic (2026-07)

Standing rules: `docs/prereg_template.md`. Certifying artifact:
`scripts/no_support_ratio.py`, `results/data/no_support_ratio.json`.

## DISCLOSURE, which must be read before the verdict

This is **NOT a blind pre-registration**, and saying so is the point of
writing it down. The aggregate counts this prediction turns on were already
computed and committed in HANDOFF 6, in
`results/data/natural_orbital_summary.json`:

    support grew      142 of 156
    support unchanged  14 of 156
    support shrank      0 of 156

Two of the 14 unchanged records are the cancellation states v89 and v103. So
before scoring anything it is already visible that roughly twelve generic
INTERFERENCE records have ratio 1, and the prediction below is therefore
expected to FAIL with about that many exceptions.

Under standing rule R1 a prediction committed after the relevant data exists
is not evidence. It is registered anyway, for two honest reasons: the
prediction was proposed independently in the handoff that opened this item,
so it should be scored rather than silently dropped; and the identity of the
exceptions is genuinely unknown and is the informative part.

What is NOT known in advance, and is the actual content of the run:

- which twelve or so records they are;
- whether they share a structural property (channel count, loop count,
  transport class, block structure) that explains ratio 1 without
  cancellation;
- the shape of the ratio distribution above 1.

## Definitions, fixed before the run

For a census record, let `s_sparse` be the support in `states.jsonl` and
`s_NO` the support of its natural-orbital representative in
`states_natural_orbital.jsonl`. For DESIGN records the two ledgers coincide,
since a design state already has a diagonal 1-RDM, so `s_NO = s_sparse` by
definition and the record is not an independent test of anything.

    no_support_ratio = s_NO / s_sparse

Amplitudes below 1e-11 are treated as zero, the same threshold the
natural-orbital ledger was built with.

## Prediction

**P-NSR-1 (as proposed in the handoff).** `no_support_ratio` equals 1
exactly on DESIGN and CANCELLATION states, and is strictly greater than 1 on
every generic INTERFERENCE state. If it held, the trichotomy would have a
one-number diagnostic.

Expectation, given the disclosure above: **FAIL**, with roughly twelve
generic INTERFERENCE records at ratio 1.

**P-NSR-2 (the half that may survive).** `no_support_ratio >= 1` on every
record, with no record whose natural-orbital representative is sparser than
its sparse one. Expectation: PASS. Registered separately because a ratio
below 1 would mean the sparse ledger was not sparse, which would be a defect
rather than a finding.

**P-NSR-3 (structural, the informative one).** The generic INTERFERENCE
records with ratio 1 are not arbitrary: they share a structural property
distinguishing them from the 142 that grow. No specific property is named in
advance, so this is scored as a QUALITATIVE FINDING, not PASS or FAIL, and
whatever property is found is POST-HOC and needs its own out-of-sample test
before it can be believed.

## Scoring rules

Verdicts from `results/data/no_support_ratio.json` alone. Scope reports raw
records and transport classes separately, per standing rule R2. The DESIGN
records are reported but excluded from the discriminating counts, since their
ratio is 1 by construction and not by measurement.

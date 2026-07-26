# The trichotomy, restated (2026-07)

Every convention claim below cites the audit artifact that verifies it. No
statement here rests on memory of how the library was built.

The 154-record diagonality audit
(`results/data/attainer_gate_audit.json`) did not weaken the
DESIGN / INTERFERENCE / CANCELLATION classification. It sharpened it, by
showing that the classes differ in whether SPARSITY and DIAGONALITY can be
achieved in the same orbital basis.

## The three classes

**DESIGN** (643 records: 631 DESIGN-INT, 12 DESIGN-REAL). Sparsity and
diagonality coexist in one basis. The support is one-hop free, the weights are
nonnegative, and no phases are needed. Because no two support determinants
differ in a single mode, every off-diagonal 1-RDM entry vanishes TERMWISE, so
the representative is simultaneously sparse and natural-orbital. This is a
basis-free existence statement about lambda, decided combinatorially by
`classify.py`. All 643 pass the diagonality gate
(`results/data/attainer_gate_audit.json`).

**INTERFERENCE, generic** (154 records). The two virtues provably split. The
shipped sparse representative has active channels and a non-diagonal 1-RDM, so
its support bounds s_Q^free; the natural-orbital representative is diagonal
but support-bloated, with no surviving closed form, and its support bounds
s_Q^NO (`results/data/states_natural_orbital.jsonl`). Neither ledger has both
properties, and that is the content of the class rather than an artifact of
how either was computed.

**CANCELLATION** (2 records: (3,10) v89 and v103). The virtues rejoin.
Sparse, diagonal, and interfering at once: phases are present and every
one-hop channel is exactly silent, so the off-diagonal entries cancel rather
than vanish termwise. These are the only two interference records that pass
the diagonality gate. Stated structurally, this is why the two rank-10
closures were the hard ones: they are the vertices where a representative had
to satisfy both demands simultaneously.

## The one-number diagnostic, scored

`no_support_ratio = s_NO / s_sparse` (`results/data/no_support_ratio.json`,
pre-registered with disclosure in `docs/prereg_no_support_ratio.md`).

| regime | records | ratio 1 | ratio > 1 |
|---|---|---|---|
| DESIGN | 643 | 643 | 0 |
| CANCELLATION | 2 | 2 | 0 |
| INTERFERENCE generic | 154 | 13 | 141 |

Generic interference ratio distribution: min 1.00, median 1.33, max 8.14.

**P-NSR-1 FAILED**, as the pre-registration disclosed it would: the ratio is
NOT strictly greater than 1 on every generic interference state. Thirteen
generic records, in twelve transport classes, sit at exactly 1.

**P-NSR-2 PASSED**: no record has ratio below 1.

**What survives is a one-directional diagnostic**, and it is clean:

> ratio > 1 implies generic INTERFERENCE. 141 records, zero exceptions.

The converse fails. Ratio 1 does not imply design or cancellation, because
thirteen generic interference states happen to admit a natural-orbital
representative no larger than their sparse one. So the support blowup IS a
measure of interference when it occurs, but its absence measures nothing.

**Structure of the thirteen exceptions (POST-HOC, no out-of-sample test).**
All thirteen have incidence kernel dimension at least 1 and a one-hop class of
size at least 2, while 93 of the 141 that grow have kernel dimension 0 and
maximum class size 1. But 48 records that grow also have kernel dimension at
least 1, so kernel freedom is NECESSARY for a generic state to stay flat and
is not sufficient. No law is claimed.

## A caveat that limits every s_Q^NO number

When lambda is degenerate the natural-orbital basis is NOT unique: any unitary
acting within an eigenspace of rho preserves diagonality. The supports in the
natural-orbital ledger come from one arbitrary choice of eigenbasis, so each
is an UPPER bound on s_Q^NO for that choice and is not minimized over the
eigenspace freedom.

This is not hypothetical. An earlier version of the ledger rotated v103 even
though it already satisfied the gate, and an arbitrary basis of its degenerate
eigenspaces reported support 79 in place of the true 14. Records that already
pass the gate are now passed through unrotated, and the ratio table above is
what results.

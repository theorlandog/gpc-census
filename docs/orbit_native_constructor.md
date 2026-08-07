# The orbit-native constructor

Status: **IMPLEMENTED and GATED against the materialising reference at ranks 7
and 8.** The fast components were each verified in isolation and none of them
was on the path that actually builds a system. They are now.

Certifying artifacts: `src/gpc_census/generation/orbit_native.py`,
`tests/test_orbit_native.py`.

## The finding this answers

Orbit enumeration, direct survivor generation and the Hall structural-zero
filter were all proved and measured, and the production constructor called none
of them. `generate_fixed_n3_reference_rank` went straight to
`enumerate_admissible_hyperplanes_by_flats`, materialised every hyperplane,
oriented each into two taus, and screened the whole list; and the screening
fallback ran the exponential determinant dynamic program with no matching
pre-test. Every measured speedup was sitting beside the pipeline rather than in
it.

## What the new path does

    orbit-canonical closed-flat enumeration
      -> one exact (h, z) per TOP-level orbit
      -> direct inversion-number survivor generation
      -> existing screening, now Hall-accelerated
      -> existing exact redundancy reduction

It never constructs the full hyperplane list, the full oriented-tau list, or a
single trace-rejected row. Downstream is untouched, because screening consumes
taus: the only job of the new module is to emit exactly the taus that survive
the trace test.

## Why coverage still holds

Two facts, both proved elsewhere in the package and neither an empirical
pattern:

1. the admissible hyperplane set is `S_d` **stable**, so orbit representatives
   plus expansion reach every hyperplane;
2. the trace test is `inv(sigma h) == B` with `B` an orbit invariant, so an
   orbit's survivors are exactly the arrangements of the multiset `h` with `B`
   inversions, and those are generated directly.

Coverage is then **checked, not asserted**, by count conservation. The left side
comes from orbit-stabilizer and the Gaussian multinomial, the right side from
what the run actually did, so a lost orbit or a dropped survivor cannot balance
it:

| system | oriented rows | never built | survivors generated | predicted in closed form |
|---|---|---|---|---|
| (3,7) | 10,682 | 10,133 | 549 | 549 |
| (3,8) | 332,840 | **326,233** | 6,607 | 6,607 |

## The trap: orientation self-pairing

Each unoriented hyperplane yields two taus, so the naive version processes both
orientations of every orbit. That double counts, because the oriented orbit of
`(h, z)` already contains `(-h, -z)` exactly when `z == 0` and the multiset of
`-h` equals that of `h`.

This is not hypothetical. There are **1, 3 and 3** self-paired orbits at ranks
6, 7 and 8, which is precisely why the oriented orbit counts are 15, 35 and 109
rather than 16, 38 and 112. The constructor detects self-pairing and processes
one orientation for those orbits, and a test pins the counts.

## The acceptance gate

Reproducing the reference is the gate, not running. A constructor that skips 98
percent of the rows is only trustworthy if it lands on the same system.

**Rank 7.** The emitted taus equal, as a SET, the 549 rows the reference keeps
after its trace test. The retained constraints, the reduction certificate and
the known-system regression are identical. The screening partition reproduces
the published numbers exactly: **1 structural, 49 certified, 499 determinant
zeros, 0 unresolved**.

**Rank 8.** 166,420 hyperplanes reconstructed from 56 orbits, 6,607 survivors
against 6,607 predicted, and the screening partition reproduces the materialised
run: **1 structural, 137 certified, 6,469 determinant zeros, 0 unresolved**.

At rank 8 every comparison passes: retained rows, reduction certificate,
known-system regression, the certified tau SET and the symbolic-zero tau SET,
with both paths reporting the system proved complete.

One field differs by design and only one: `screening_index`, a row's position in
the screened list, which cannot agree when one path screens 549 rows and the
other 10,682. It is bookkeeping, not mathematics, and the tests compare with it
removed rather than quietly passing.

**Cost, both paths measured back to back on an otherwise idle machine:**

| system | reference | orbit-native | speedup |
|---|---|---|---|
| (3,7) | 7.1 s | 3.9 s | 1.8x |
| (3,8) | 251.9 s | **66.0 s** | **3.8x** |

The speedup is far smaller than the 719x that survivor generation shows in
isolation, and the reason is worth stating rather than hiding: the orbit-native
path still pays for the orbit ENUMERATION, which is now the dominant stage. The
row scan is gone; what replaced it as the bottleneck is enumerating the orbits
themselves. That is exactly where the remaining work is.

## The Hall filter, on the production path

The structural test now runs inside the symbolic fallback, before any
determinant algebra. Its outcome is deliberately **indistinguishable** from the
one the dynamic program returns, same status and same reason, so every stored
manifest stays byte identical and this is a pure accelerator rather than a
change of verdict. That was confirmed directly: the rank-7 reference manifest
regenerates byte-for-byte, and the rank-8 screening counts match the stored
artifact exactly.

The consequence for work done, at rank 8: of the 6,469 rows that reach the
fallback, **6,414 are decided by matching** and only **55** reach the
exponential dynamic program.

The Hall violator is recorded by the orbit-native manifest rather than the
legacy row payload. Adding a field to the row would rewrite every published
manifest for no mathematical gain.

## What this does and does not license

**Does.** A complete system built from orbit representatives alone, gated
against the materialising path at both ranks where that path can still run. The
coverage argument is a theorem plus a conserved count, not a resource claim.

**Does not.** It does not enumerate faster than the orbit enumerator can go, so
rank 11 still waits on that. It does not touch the residual determinant
algebra beyond removing the structural cases. It does not make screening an
`S_d` invariant, which is refuted and stays refuted.

## Reproducing

```sh
uv run pytest tests/test_orbit_native.py
```

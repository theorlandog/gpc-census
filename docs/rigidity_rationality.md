# Rigidity and rationality (2026-07)

Status: the design case is PROVED below and covers 643 of 799 states. The
implication "fixed-rho rigid implies rational weights" is CONFIRMED BY CENSUS
at 740 of 740 with a power warning that materially limits how it may be
quoted. The general theorem is a TARGET with two stated gaps. The v_B Galois
embedding test is CLOSED, its prediction CONFIRMED, and its earlier naive
form is retracted.

Certifying artifacts: `scripts/attainer_gate_audit.py`,
`results/data/attainer_gate_audit.json`,
`scripts/rigidity_rationality_census.py`,
`results/data/rigidity_rationality.json`, `tests/test_rigidity_rationality.py`.

## Correction first: the diagonality gate fails on 154 records, not one

The motivating handoff reports that the shipped v_B record at (4,9) index 65
violates the attainer gate (1-RDM not diagonal, max off-diagonal 3/23) and
instructs: "Expect exactly one; verify that expectation." The expectation is
wrong, and the verification is the point of stating it that way.

Running `gpc_census.validate.check_vertex_attainer` over all 799 records:

    pass 645    DESIGN-INT 631, DESIGN-REAL 12, INTERFERENCE 2
    fail 154    INTERFERENCE 154

The 154 failures are every INTERFERENCE record except (3,10) v89 and v103.
v_B's 3/23 is real and reproduces exactly, but v_B is an ordinary member of
that population, not a singleton defect. Verified independently in exact
sympy arithmetic at (3,10) v17, whose 1-RDM has four nonzero off-diagonal
entries and whose diagonal is not even a permutation of lambda, while its
eigenvalues are exactly lambda.

Three facts keep this in proportion:

1. **Every one of the 154 still attains the spectrum**, to worst error
   1.1e-15, at unit norm. They are correct states in the wrong basis, not
   wrong states.
2. **The library never certified diagonality for them.** The shipped
   certificate is the characteristic-polynomial identity
   (`scripts/verify_states_standalone.py`), which certifies
   spec(rho) = lambda; diagonality is asserted only for design states, where
   it follows from one-hop freedom by inspection.
3. **The repository already records the same fact from the other side.**
   RESEARCH.md states that the cancellation regime, defined as an exactly
   diagonal 1-RDM with silent channels, "has exactly TWO members across all
   799 certified representatives, v89 and v103", and `cascade.jsonl` carries
   a per-state `rdm_max_offdiagonal` that is nonzero throughout the magnitude
   regime.

So this is a property of how the interference library was built, not a
regression, and not something the gate "was bypassed" to allow: the gate
postdates the library and encodes a stricter invariant than the library was
ever built to satisfy.

**What genuinely needs repair is a documentation claim.** RESEARCH.md's
support taxonomy says s_Q^NO is "the census quantity, the one the library
column bounds". For the 154 non-diagonal records the support column bounds
s_Q^free and not s_Q^NO, which is exactly the conflation the taxonomy was
written to prevent.

### The repair: all 156 interference records, as an additive ledger

Decision taken 2026-07: repair all of them, not just v_B. Repairing one record
while 153 identical cases remain would imply a library-wide fix that had not
happened.

`scripts/natural_orbital_ledger.py` produces
`results/data/states_natural_orbital.jsonl`: for every interference record,
diagonalize its 1-RDM as rho = U diag(w) U^H with w ordered to match lambda,
and push the state through the N-th compound of U^T. Every one of the 156
resulting states has 1-RDM exactly diag(lambda) (worst off-diagonal 2.8e-16,
worst diagonal error 7.8e-16) and passes the gate, re-verified from the
shipped amplitudes in `tests/test_natural_orbital_ledger.py` rather than
trusted from the generator. The 2-RDM spectrum, which is basis-free, is
unchanged to 1.6e-15, so the rotation moved the basis and not the state.

A convention note worth recording, because getting it wrong is silent: the
correct compound is of the TRANSPOSE U^T, not the conjugate transpose. The two
agree for real rotations, so an incorrect convention still appears to work on
the real records and fails only on the complex ones.

**The repair is not free, and the cost is why it is additive rather than a
replacement.**

| | states.jsonl | states_natural_orbital.jsonl |
|---|---|---|
| amplitudes | EXACT symbolic closed forms | NUMERICAL |
| support, mean / max | 6.8 / 14 | 11.7 / 57 |
| 1-RDM | not diagonal on 154 records | diag(lambda) on all 156 |
| support bounds | s_Q^free | s_Q^NO |

Rotating into natural orbitals destroys sparsity: support grows on 141 of the
156 records, is unchanged on 15, and never shrinks. (Records that already pass
the gate are passed through unrotated: with a degenerate lambda, eigh returns
an arbitrary basis of each eigenspace rather than the identity, so rotating an
already-diagonal record scrambles it for nothing. An earlier version of this
ledger did exactly that and reported v103 at support 79 instead of 14. The
same non-uniqueness means every s_Q^NO number here is an upper bound for ONE
choice of eigenbasis, not a minimum over the eigenspace freedom.) The exact closed forms do
not survive the rotation, and re-recognizing them is a research task rather
than a mechanical one, so the new ledger is numerical.

So `states.jsonl` is NOT replaced. The two ledgers answer different questions
and both are needed: the exact sparse library for everything that depends on
closed forms (holonomy fields, denominator arithmetic, the whole arithmetic
layer), and the natural-orbital ledger for any claim about s_Q^NO. This also
honours the originating instruction to keep the existing record as a
documented spectrum-attainer rather than overwrite it.

**Consequence for support claims.** Every vertex now has a natural-orbital
support bound, but for the 154 it is LARGER than the number in the census
column. Any s_Q^NO claim must quote the natural-orbital ledger; the census
column remains correct for s_Q^free.

## The census result

`results/data/rigidity_rationality.json`. Rigidity is NUMERICAL, inherited
from `cascade.jsonl`; the weight field is EXACT, recomputed from the closed
form rather than read off the integer-versus-string encoding.

| population | states | rigid | deformable | rigid and rational | counterexamples |
|---|---|---|---|---|---|
| all shipped records | 799 | 740 | 59 | **740** | **0** |
| gate-passing only | 645 | 645 | 0 | 645 | 0 |

Every fixed-rho first-order rigid record has rational squared weights, with
zero counterexamples. The single irrational record is the current v_B, whose
squared weights lie in Q(sqrt 35), and it is deformable (fixed-rho tangent 1).
Direction check: 58 of the 59 deformable records are rational anyway, so
deformability permits irrationality without forcing it and the implication is
one-way.

### The power warning, which is not optional

**Every one of the 645 natural-orbital representatives is fixed-rho rigid,
and all 59 deformable records are among the 154 that fail the gate.**

So the implication's discriminating power lives entirely in the records that
are not natural-orbital representatives. On the gate-passing subset the
antecedent is always true, and there is no rigid-versus-deformable contrast
for the implication to be confirmed against. Quoting "740 of 740" as evidence
about VERTICES overstates it in exactly the way the (3,11) bracket holdout's
38 rows did: the count is correct and the independent content is smaller
(`docs/prereg_template.md`, standing rule R2).

This is expected rather than surprising. The fixed-rho fiber fixes rho itself,
so a record whose rho is not diag(lambda) is constrained differently from one
whose rho is; the two populations are not comparable observations of a single
law.

## The design case, PROVED

This covers 643 of the 799 states outright, by linear algebra.

**Theorem (design case).** Let S be a ONE-HOP FREE support, meaning no two
determinants of S differ in exactly one mode, and let psi = sum_t c_t |t> be a
state on S whose 1-RDM equals diag(lambda) with lambda rational. Then the
squared weights p_t = |c_t|^2 are rational, provided the incidence system has
a unique solution.

*Proof.* One-hop freedom means that for p != q no pair of support
determinants is related by replacing p with q, so every off-diagonal entry
rho_pq receives no contribution at all and vanishes termwise; rho is diagonal
regardless of the phases. Its diagonal is
rho_mm = sum_{t contains m} |c_t|^2 = (A_S p)_m, where A_S is the 0/1
incidence matrix of the support. So rho = diag(lambda) is exactly the linear
system A_S p = lambda, together with normalization sum_t p_t = 1. Both have
rational coefficients and a rational right-hand side. The solution set is
therefore a rational affine subspace of Q^S; when it is a single point that
point is rational. ∎

**Census check.** All 643 design states have incidence kernel dimension 0 and
zero one-hop channels (`cascade.jsonl`), so the uniqueness hypothesis holds
throughout and the theorem applies to every one of them. All 643 do have
rational weights.

Note what the theorem does and does not need: it needs no rigidity input and
no Galois argument, only uniqueness of a rational linear solve. The
interference case cannot be reached this way, because there the phases enter
and the defining system is no longer linear in p. Those 156 states are the
open part.

## The theorem target

> **RIGIDITY-RATIONALITY (target, not proved).** If the support-restricted
> fixed-rho fiber of a census vertex, modulo gauge, is a single point (true
> isolation, not merely first-order rigidity), then that point is defined
> over Q, because the defining equations have rational coefficients and a
> Galois-stable singleton is Galois-fixed.

Two gaps, and they are the work.

**(a) First-order rigidity is not isolation.** A vanishing tangent space is
necessary for isolation, not sufficient; the fiber could be a higher-order
tangency. Upgrade path: stress-matrix or second-order certificates per the
rigidity-theory entries in `docs/prior_art_roadmap.md`, or exact isolation
certificates on the polynomial system for a sample and then census-wide.
Until then every "rigid" in the table above means first-order rigid.

**(b) Isolation gives finiteness, not uniqueness.** MEASURED 2026-07,
`results/data/uniqueness_census.json` and `results/data/v103_fiber_count.json`.
Multi-starting the fixed-rho realization system on twelve rigid interference
states: every one has exactly ONE solution modulo gauge and conjugation, and
exactly one weight vector. The pre-registered uniqueness claim PASSES 12 of 12.

At v103 this is certified rather than sampled. Once the weight vector is known
to be common to every solution, the unknowns are phases alone, and modulo the
9-dimensional orbital-phase gauge they leave exactly 5 degrees of freedom, one
per loop. That search space is a compact 5-torus and can be exhausted: 16807
grid starts, all converging, yield ONE point, with holonomy (0, pi, pi, 0, 0),
so the state is real up to gauge exactly as the library record says.

CORRECTION, recorded because it was published in the other direction first.
This census originally scored v103 as an exception with dozens of points, and
that was an artifact of the measurement. The loop basis had been taken from an
SVD, so its vectors had irrational entries. Amplitude phases are read wrapped
into (-pi, pi], and a gauge change followed by that wrap shifts theta by
2*pi*m for an integer m; for an INTEGER kernel vector the holonomy then moves
by a multiple of 2*pi and its cosine is unchanged, but for a real one it moves
arbitrarily. One physical point was being split into dozens by branch choice.
So gap (b) has no counterexample anywhere it has been measured.

The theoretical statement of the gap is unchanged: Galois may permute several
conjugate isolated points, giving weights algebraic of degree equal to the
orbit size. So the honest statement is:

> rigid implies the weights are algebraic of degree at most the number of
> isolated fiber points modulo gauge, with rationality when that number is 1.

The 740/740 then says that uniqueness-modulo-gauge is empirically the norm,
which is a SEPARATE claim needing its own test (the uniqueness-modulo-gauge
census, not yet run).

v_B is the instructive case for (b): its weights have degree exactly 2 over Q,
generated by sqrt(35), and it sits at the endpoint of a wall family cut out by
the quadratic 85 t^2 - 70 t - 119, whose two roots are a conjugate pair with
one physical. Orbit size 2, degree 2. That is the degree-equals-orbit-size
mechanism visible in the one record that has it
(`docs/holonomy_rationality.md`, "The one non-rational record").

## The v_B Galois embedding test: earlier form RETRACTED

A naive substitution sqrt(35) to minus sqrt(35) on the v_B record fails
attainment with spectrum error 0.068. That failure is diagnostic, not
refuting, and the earlier reading of it is withdrawn.

The amplitudes contain NESTED surds such as sqrt(120 - 18 sqrt(35)), and
independent principal-branch substitutions on the visible radicals do not
constitute a field embedding when the surds satisfy multiplicative relations.
Since the fiber equations are polynomial over Q, any genuine embedding carries
exact solutions to exact solutions; a substitution that breaks attainment is
therefore proof that it was not an embedding, and says nothing about the
geometry.

The proper test is: construct the number field K generated by all amplitude
surds, compute its real embeddings, apply each embedding consistently to the
whole amplitude vector, and test attainment exactly, with the expected orbit
size pre-registered before step 3.

**RUN and CLOSED**, `results/data/vb_galois_orbit.json` via
`scripts/vb_galois_orbit.py`. The pre-registered orbit size 2 is CONFIRMED:
both roots of the wall quadratic `85 t^2 - 70 t - 119` give all-positive
weights and attain exactly, the shipped root `7/17 - 18 sqrt(35)/85` with 64 of
128 sign patterns attaining. This is the one record where degree equals orbit
size is visible, so it is the worked example for the whole
degree-equals-orbit-size statement. (This section previously read "specified
but NOT YET RUN" after the test had been executed and scored; the ledger entry
was correct and this document had drifted.)

## Cross-references

- `docs/holonomy_rationality.md` takes weight rationality as an explicit
  HYPOTHESIS (H1). This document is where that hypothesis gets its grounding:
  proved for the 643 design states, census-confirmed but power-limited for the
  rigid ones, open in general.
- `docs/prereg_denominator_arithmetic.md` has the same dependency and excludes
  v_B for exactly this reason: with irrational weights it has no state
  denominator.
- `docs/two_block_family.md`: v89 and v103 are the rigid-and-rational
  flagships, and they are also the only two interference records that pass the
  diagonality gate, which is not a coincidence but the definition of the
  cancellation regime.

## Open

1. The repair decision for the 154 non-diagonal records.
2. Gap (a), isolation certificates.
3. Gap (b), the uniqueness-modulo-gauge census.
4. Whether the interference case admits any analogue of the design-case
   proof, given that phases make the system nonlinear in p.

Item 4 of the earlier list, the v_B embedding test, is CLOSED and confirmed;
see the section above.

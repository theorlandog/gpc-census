# HANDOFF 4: the gauge session (2026-07)

Supersedes HANDOFF 3 in two places and confirms it in one. Read this before
acting on anything HANDOFF 3 told you to do; two of its instructions were built
on a premise that turns out to be false, and following them would put a
retracted claim back into the paper.

Everything below is re-derivable from the repo. Where I say "checked", the
artifact is named.

## The one-line version

HANDOFF 3 was right that the natural-orbital gauge had never been explored, and
right that the cascade is mostly an orbital rotation engine. It was wrong that
the v103 support-10 endpoint is diagonal, and therefore wrong to reinstate
s_Q^NO(v103) <= 10. Exploring the gauge properly, the library support column
does not move anywhere in the census: 0 of 799 states sparsify and 629 of 799
are certified gauge-minimal outright.

## Correction A: the v103 endpoint is NOT diagonal, and the bound does not stand

HANDOFF 3 says "The v103 endpoint is diagonal to 1.3e-17". It is not. The
shipped support-10 cascade endpoint has

    max |rho - diag(lambda)| = 0.0735294... = 5/68,

exactly, and its diagonal is not lambda either (it reads 9/17, 233/442, 9/17,
35/68, ...). The 1.3e-17 figure is the diagonality of the **library support-14
state**, which I measured at 1.328147661294743e-17. The two got transposed
somewhere between sessions.

The consequence is the opposite of HANDOFF 3's: the support-10 state is not a
natural-orbital representative, so **s_Q^NO(v103) <= 10 does not stand**. Its
eigenvalues are lambda exactly, so it does attain the vertex; it attains it in
another basis. The bound it certifies is on the free-basis minimum.

The certificate itself was never broken. HANDOFF 2 item 4 said the recognized
weights "force an off-diagonal of 5/68 under every one of the 2^10 sign
patterns, so they are not the endpoint's weights". They are the endpoint's
weights. The endpoint really does carry that off-diagonal; the recognition was
correct and the arithmetic is sound. What was wrong was the LABEL. So do not
"redo the recognition with diagonality imposed" as HANDOFF 3 instructs: there is
nothing to fix there, and imposing diagonality would be imposing a condition the
state does not satisfy. I relabelled it instead, and added the exact
non-diagonality as a checked field so the label cannot silently revert
(`scripts/v103_support10_certificate.py`, `tests/test_v103_support10.py`).

This generalizes, and that is the more useful finding: **none of the 71 cascade
endpoints has a diagonal 1-RDM. Zero of 71.** The fixed-spectrum continuation
preserves spec(rho), which is all it ever promised, and does not preserve
diagonality. So every improved cascade bound in the repo was a free-basis bound
wearing a natural-orbital label. All of them are relabelled
(`results/data/orbit_check.json`, `docs/cascade_census.md`).

## Correction B: 8 of the 71 endpoints ARE new states

HANDOFF 3 tested 26 endpoints, found every one matched its library state's
2-RDM spectrum, and concluded "Not one is a new state", instructing: "Do not
describe any of the 70 as new extremal states."

Over all 71, 63 match and **8 do not**, by 3.4e-3 to 4.7e-2. Differing 2-RDM
spectra prove no one-body unitary connects them, so those 8 are genuinely
distinct extremal states over their vertices, by HANDOFF 3's own gate. They are
(4,10) v93, (5,10) v65, v256, v144, v140, v230, v263, and (4,9) v65, which is
v_B. The split is bimodal with nothing between 1e-12 and 1e-3, and all 8 have
certificate residuals ~1e-17, so this is not a tolerance artifact.

Worth noticing: two of the 8 are (5,10) v144 and v256, exactly the endpoints
that resisted rational recognition in the earlier session. That is now explained
rather than mysterious. They are different states, so there was no reason for
them to sit at the library state's weights.

## Confirmed: the cascade is mostly an orbital rotation engine, and now provably

For the other 63, the connecting one-body unitary is exhibited, not inferred.
HANDOFF 3 asked for this at v103 "before it is quoted as certain"; it is done
for all 63. Worst residual 1.7e-15, worst unitarity error ~4e-15, matrices
shipped in `results/data/orbit_check.json`. At v103 the rotation's off-block
entries reach 0.196, so it is a genuine U(10) rotation and not an element of
U(4) x U(6).

## The new task, done: gauge minimization returns nothing

HANDOFF 3 predicted "corrected s_Q^NO improvements from the gauge
minimization". There are none.

    states searched                        799
    states whose support the gauge reduces   0
    certified gauge-minimal                629
    open in both directions                170

The certification is exact and needs no search. A block-diagonal unitary
preserves each determinant's BLOCK PROFILE (|T n B_1|, ..., |T n B_k|), so its
compound preserves the profile-class decomposition and each class's flattening
ranks. Hence

    support of any gauge representative >= sum over classes of flattening rank,

and where the library support attains this, the state is gauge-minimal by an
integer rank computation. That bound is the most reusable thing this session
produced.

Read the null result at its real strength. The search's power is pinned by a
recovery test rather than asserted: a provably rank-1 class must collapse to one
determinant, and the search finds that collapse exactly, at residual 4e-16, with
a rank-2 control where no collapse exists. Two earlier methods FAILED that test
and were deleted rather than shipped: a reweighted-L1 Jacobi sweep, and a
Riemannian descent whose smoothed-L1 surrogate made growing a zero amplitude
free, so it walked away from the sparsest point it was handed (regression test
kept). If you write a new search, run the power test first. I nearly reported a
null result from the sweep before testing it, and it would have been worthless.

Limits: multi-start greedy over a 27-to-100-dimensional curved group is not
exhaustive, and a reduction needing two simultaneous cancellations would be
missed. The 629 certificates do not depend on the search; the 170 open rows do.

## The taxonomy, which is what actually prevents a repeat

This whole mess is one conflation, so `docs/RESEARCH.md` now carries four
distinct minima and the lower bounds that do and do not apply:

- s_M: majorization baseline.
- s_I: diagonal incidence baseline. s_M <= s_I.
- s_Q^NO: least support with rho EXACTLY diag(lambda). The census quantity, the
  one the library column bounds. **s_I <= s_Q^NO.**
- s_Q^free: least support with spec(rho) = lambda. s_M <= s_Q^free <= s_Q^NO.
  **s_I does NOT bound s_Q^free.** Never put them in one inequality.

Plus a per-state fifth: the orbit-minimal support of a given psi over its U(d)
orbit, which is what the 63 same-state endpoints measure.

The gate is now code, not prose: `gpc_census.validate.check_vertex_attainer`
requires diagonal 1-RDM, diagonal equal to lambda, unit norm, and, for anything
claimed NEW, a 2-RDM spectrum differing from the library state's. It would have
caught the mislabelled v103 bound immediately.

## Pre-registration, scored honestly

`docs/prereg_gauge_minimization.md` was committed before the run (check the git
history; that was the point). Because DROP came out identically zero, P-G1 and
P-G2 are NOT SCOREABLE for want of positives, which the pre-registration named
in advance as the mooting outcome. P-G4 FAILED: 170 states have slack between
support and lower bound and not one sparsifies, so the slack predicts nothing
except that the bound is loose there.

So HANDOFF 3's open question, whether the support-restricted fixed-rho tangent
predicts gauge sparsifiability, is still open and now looks unaskable in this
corpus: there is no variation in the outcome variable to predict.

## What I would do next, in order

1. **Sharpen the lower bound on the 170 open rows.** 51 of them are open only
   because no invariant is implemented for a Lambda^a factor with
   3 <= a <= m-3, where I fall back to the trivial value 1. A Lambda^3
   normal-form invariant is a finite classical computation and would likely
   close most of them. This is the highest value per unit effort in the repo
   right now.
2. **Exactly recognize the 8 new states.** They are new attainers, never
   analyzed, and two of them are the pair that defeated the earlier recognition
   pipeline. Their 2-RDM spectra are new data about those vertices.
3. **Decide whether s_Q^NO can beat the library support at a DIFFERENT fiber
   component.** Everything I certified is minimality within one state's gauge
   family. s_Q^NO is a vertex-level quantity and nothing here bounds it below
   except s_I. That gap is now the honest frontier.
4. The companion note (`docs/no_compactness_note.md`) is written and its data
   exist: 48 certified counterexamples to "natural orbitals give the most
   compact CI expansion", each with its explicit rotation. What it lacks is an
   exact closed form for U; the weights machinery should transfer.

## Things I did not do

- I did not touch the 33 multiplicity-drift stops. That finding is sound and
  unaffected.
- I did not re-run the cascade. Its numbers are unchanged; only their labels
  moved.
- I did not attempt the second-order structure at v103, or anything else in the
  fiber program proper.

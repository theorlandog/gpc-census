---
bibliography: results/report/level5/references.bib
csl: results/report/physics-numeric.csl
link-citations: true
abstract: |
  Altunbulak and Klyachko state, without proof, a family of level-five inequalities extending the second-level quadruples of the three-fermion moment polytope by the last occupation number. We give exact Ressayre certificates for the four instances at eleven orbitals, so those four inequalities are proved valid. Certification is a finite exact computation: a trace condition on weights below an affine hyperplane, and a nonzero integer tangent determinant obtained by fraction-free elimination. The certificates are replayed by a verifier that shares no code with the generator, derives the fermionic sign by a different route, and recomputes each determinant by a different exact algorithm. The screen is calibrated at the two largest ranks whose systems are published, by a complete test rather than a sample: every certified row is maximized by linear programming over the published polytope, and none is false. As a consequence, the four spectra that remained undecided in an outer bracket for three fermions in eleven orbitals are refuted, and that bracket closes at nineteen attained and thirty-one refuted. Validity is what is proved. Facetness, irredundancy, and completeness of the eleven-orbital system are not claimed, and are not needed for the refutations.
keywords:
  - generalized Pauli constraints
  - moment polytopes
  - Ressayre certificates
  - quantum marginal problem
  - fermionic one-body marginals
author:
  - |
    James Orlando\
    Independent researcher, Derry, New Hampshire, USA\
    `jamie@orlandonh.com`
title: |
  Exact Ressayre certificates for the level-five inequalities of the three-fermion moment polytope at eleven orbitals
---

# Introduction

For a normalized $\psi\in\wedge^N\mathcal H_d$ the one-body reduced density matrix is $\rho^{(1)}_{ij}(\psi)=\langle\psi|a_i^\dagger a_j|\psi\rangle$, and its ordered eigenvalues $\lambda_1\ge\cdots\ge\lambda_d$ are the natural occupation numbers. Klyachko's solution of the pure-state one-body $N$-representability problem identifies the attainable spectra with a rational moment polytope [@Klyachko2006; @AK2008], whose nontrivial facets are the generalized Pauli constraints. The three-fermion, six-orbital case is due to Borland and Dennis [@BD1972].

Complete inequality descriptions are published for $N=3$ through $d=10$ [@AK2008; @AltunbulakThesis]. At $d=11$ no complete description is available. What exists is a family of *level-five* inequalities, stated in Remark 4.2.1 of Ref. [@AK2008] and in Ref. [@Klyachko2009], which extend each second-level quadruple by the tail occupation number. In the eleven-orbital case the family has four members,
$$
\lambda_A+\lambda_{11}\le 2,
\qquad
A\in\bigl\{\{2,3,4,5\},\{1,3,4,6\},\{1,2,5,6\},\{1,2,4,7\}\bigr\},
\tag{1}
$$
writing $\lambda_A=\sum_{i\in A}\lambda_i$. Both sources state these without proof. We refer to the $A=\{1,2,4,7\}$ instance as AK-RMK-4.2.1 and to the other three as KLY09-EXT-1, -2 and -3, following the labels used in the accompanying data.

This note proves all four.

**Theorem 1.** *Each of the four inequalities in (1) is valid on the moment polytope of $\wedge^3\mathcal H_{11}$.*

The proof is a finite exact computation, one Ressayre certificate per inequality, described in Sec. 2 and recorded row by row in Sec. 3. Section 4 is the independent replay, Sec. 5 the calibration and its negative control, and Sec. 6 the consequence for an outer bracket at eleven orbitals.

We are deliberate about what is and is not claimed. A Ressayre certificate proves that the single inequality it certifies is valid. It does not establish that the inequality is a facet, that the four together are irredundant, or that any list of inequalities at $d=11$ is complete. None of those is needed for the application in Sec. 6: a spectrum violating a valid inequality lies outside the polytope, and that is the whole of the argument.

# The criterion, and what a certificate proves

Ressayre's theorem [@Ressayre2010], building on Belkale and Kumar [@BelkaleKumar2006] and Berenstein and Sjamaar [@BerensteinSjamaar2000], characterizes the inequalities of a moment cone through a condition on an affine hyperplane in the weight polytope together with a transversality condition on an associated tangent map. We use it in the following concrete form for the multiplicity-free representation $\wedge^N\mathbb C^d$.

Let $W\subset\{0,1\}^d$ be the $\binom{d}{N}$ determinant weights, each the indicator vector of an $N$-subset. Let $(h,z)$ be a primitive integer pair defining the affine hyperplane $\{x: h\cdot x=z\}$, and set
$$
\mathrm{On}=\{w\in W: h\cdot w=z\},\qquad
\mathrm{Below}=\{w\in W: h\cdot w<z\}.
$$
Let $R^-=\{(i,j): i<j,\ h_j-h_i<0\}$ be the negative roots the hyperplane makes negative, with root operators $E_{ji}=a_j^\dagger a_i$ acting on determinant weights with the usual fermionic sign. The **tangent matrix** is the $|\mathrm{Below}|\times|R^-|$ matrix whose $(w,\alpha)$ entry collects, over $u\in\mathrm{On}$, the signed action of the root $\alpha$ carrying $u$ to $w$, weighted by an indeterminate attached to $u$.

The certificate consists of two exact conditions.

1. **Trace condition.** $|\mathrm{Below}|=|R^-|$, so the tangent matrix is square. This is a finite count and is decided exactly.
2. **Nonvanishing.** The determinant of the tangent matrix, evaluated at one integer point, is nonzero. The evaluation is deterministic and the determinant is computed by fraction-free Bareiss elimination on exact integers, so a nonzero value is a proof and never a rounding artifact. A single nonzero evaluation certifies that the determinant polynomial is not identically zero, which is the condition the criterion requires; exhausting a finite set of trial points without a nonzero value is recorded as inconclusive and never as vanishing.

The convention throughout is that a certificate $(h,z)$ proves $h\cdot\lambda\ge z$. An inequality $\lambda_A\le k$ on the trace-$N$ slice is converted to this form by $\tau_i=N\,[i\in A]-k$, then $S=\sum_i\tau_i$, $h_i=S-d\,\tau_i$, $z=NS$, made primitive. On $\sum_i\lambda_i=N$ the two statements are equivalent, and Sec. 4 checks that equivalence numerically on random rational points of the slice as well as through the closed form.

# The four certificates

Screening the four rows of (1) at $d=11$ returns `certified` for all four. Each has $|\mathrm{On}|=60$ of the $\binom{11}{3}=165$ weights on the hyperplane, and satisfies the trace condition with $|\mathrm{Below}|=|R^-|=10$, so each tangent matrix is $10\times 10$. The recorded determinant values are

| label | $A$ | $\lvert\mathrm{Below}\rvert=\lvert R^-\rvert$ | tangent determinant |
| --- | --- | ---: | ---: |
| AK-RMK-4.2.1 | $\{1,2,4,7\}$ | 10 | $99\,369\,851\,016\,234\,708$ |
| KLY09-EXT-1 | $\{2,3,4,5\}$ | 10 | $-2\,652\,958\,389\,212\,757\,152$ |
| KLY09-EXT-2 | $\{1,3,4,6\}$ | 10 | $40\,995\,952\,218\,993\,520$ |
| KLY09-EXT-3 | $\{1,2,5,6\}$ | 10 | $-62\,503\,393\,702\,267\,752$ |

Every entry is a nonzero integer, so Theorem 1 follows. The sign of a determinant carries no meaning here; only its nonvanishing does.

# Independent replay

A certificate produced and checked by the same code is a weak object. The determinants above were produced by the project's generation package, so they are replayed by a separate verifier that imports nothing from that package.

The verifier rebuilds, from the recorded $(h,z)$ and the recorded evaluation point alone, the 165 weights of $\wedge^3\mathbb C^{11}$, the sets $\mathrm{On}$ and $\mathrm{Below}$, the negative root set, and the tangent matrix, then recomputes the determinant. Two steps use deliberately different derivations, so that agreement is evidence rather than tautology.

- The fermionic sign is obtained by substituting the created orbital into the ordered occupation tuple and counting the inversions that resorting removes, rather than by tracking annihilation and creation positions separately as the generator does.
- The determinant is recomputed by exact Gaussian elimination over $\mathbb Q$, rather than by the generator's fraction-free Bareiss routine.

The recorded index sets are compared against the recomputation but never used in it, so a tampered $\mathrm{On}$ set is detected rather than propagated. Sixty-four checks pass. The verifier also fails when it should: perturbing the recorded determinant, hyperplane, level set, evaluation point, or root set is caught in each case, and two of those mutations are asserted in the test suite, since a verifier that cannot fail proves nothing.

# Calibration, and an honest negative control

Theorem 1 rests on the criterion being correctly implemented at $d=11$. That is calibrated at $d=9$ and $d=10$, the two largest ranks whose $N=3$ systems are published [@AK2008; @AltunbulakThesis], by a test that is complete rather than sampled: every subset row $\lambda_A\le k$ with $2\le|A|\le 7$ and $k\in\{1,2\}$ is screened, and every certified row is then maximized by linear programming over the published polytope. A maximum at or below $k$ proves the certified row valid on the whole polytope, not merely at sampled points.

| $d$ | rows screened | certified | invalid | verbatim in published table |
| ---: | ---: | ---: | ---: | ---: |
| 9 | 984 | 19 | 0 | 4 |
| 10 | 1{,}914 | 31 | 0 | 5 |

Fifty certified rows, none false. That most of them do not appear verbatim in the published tables is expected: they are valid but redundant, which is what a validity screen making no facetness claim should produce.

The same runs supply the negative control. A screened row whose linear-programming maximum exceeds its right-hand side is demonstrably false, and 1,857 of the screened rows are false. None was certified.

That statement should be read with one qualification, which we state rather than omit. Two of the rejection routes are structural shortcuts that never form the tangent determinant: a row can fail to define a codimension-one weight hyperplane, or fail the trace condition. Splitting the false rows by route:

| $d$ | false rows | inadmissible | trace-rejected | reached the determinant | certified |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 9 | 626 | 36 | 587 | 3 | 0 |
| 10 | 1{,}231 | 0 | 1{,}228 | 3 | 0 |

The trace condition does nearly all of the rejecting. Only six false rows in the whole calibration reached the tangent determinant, and all six were rejected there on an exactly zero determinant. A targeted attempt to manufacture more, by lowering the right-hand side of each of the 145 published rank-9 and rank-10 facets by one and by two, produced 290 further false rows, every one of which was rejected before the determinant. A seeded sweep of 6,000 rows at $d=9$ with coefficients drawn from $\{-1,0,1,2\}$ and right-hand sides from $\{0,1,2,3\}$ certified five rows, all valid, among 3,791 demonstrably false ones, and contributed one further false row that reached the determinant and died there.

So the calibration establishes that the screen as a whole certified nothing false across 8,898 screened rows, and that admissibility together with the trace condition already excludes essentially every false row of these shapes. It does not stress the determinant step heavily. We record that asymmetry because Theorem 1 depends on the determinant step, and a reader is entitled to know which part of the pipeline the calibration exercised least.

# Consequence: an outer bracket at eleven orbitals closes

The accompanying data set carries a fifty-point outer bracket for $\wedge^3\mathcal H_{11}$. Forty-six of its points were previously settled without any $d=11$ inequality system: nineteen attained, by transport of certified states from the rank-9 and rank-10 census, by an explicit $\mathbb Z_{11}$ difference design for the uniform point $(3/11)^{11}$, and by a frozen-core lift; twenty-seven refuted, by frozen-core reduction to $N=2$ together with the even-degeneracy pairing property of two-fermion one-body marginals, and by restriction to lower rank against known systems. Four points remained undecided, precisely because each violates one or more of the inequalities in (1) and none violates any previously proved inequality.

With Theorem 1 those four are refuted:

| index | spectrum | violates | value |
| ---: | --- | --- | ---: |
| 23 | $(6/7)^2(1/7)^9$ | AK-RMK-4.2.1, KLY09-EXT-3 | $15/7$ |
| 26 | $(7/9)^2(11/27)^3(1/27)^6$ | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-3 | $55/27$ |
| 34 | $(11/17)^4(1/17)^7$ | AK-RMK-4.2.1, KLY09-EXT-1, KLY09-EXT-2 | $35/17$ |
| 44 | $(1/2)^5(1/12)^6$ | KLY09-EXT-1 | $25/12$ |

Each value exceeds $2$ by exactly one over the spectrum's denominator. The bracket therefore closes at nineteen attained and thirty-one refuted, with none undecided.

Two consistency checks accompany this. None of the nineteen attained spectra violates any of the four certified inequalities, and neither does the Slater point $(1,1,1,0,\dots,0)$ nor the uniform point $(3/11)^{11}$. At a rank whose polytope is unknown this is the only soundness check available, and it is reported as such.

Index 44 deserves a separate remark. It violated exactly one of the four inequalities, which made attaining it the cheapest possible falsification of a published claim in this family, and therefore a natural target. That target is closed from both directions. A support in which no two determinants share exactly two orbitals is a partial Steiner triple system, and for such a support the one-body marginal is diagonal with $\lambda_p=\sum_{S\ni p}w_S$; an exact counting argument shows no such support realizes index 44. Writing $a,b,c,d$ for the total weight on determinants meeting the five orbitals at $1/2$ in three, two, one and zero places, the incidence identities give $3a+2b+c=5/2$, $a+b+c+d=1$ and $b+2c+3d=1/2$, hence $a-c-2d=1/2$ and so $a\ge 1/2$. At most two determinants fit inside a five-set under the Steiner condition, and each of the three resulting cases fails on covering the six orbitals at $1/12$. With KLY09-EXT-1 now certified, no state of any kind attains it.

# What is not proved

The calibration is complete at $d=9$ and $d=10$. Soundness of the implementation at $d=11$ follows from the criterion itself, which is rank-independent, together with the evidence that the implementation realizes it correctly at the ranks where the answer is independently known. It is not an enumeration at $d=11$, and we keep that distinction.

Nothing here enumerates the eleven-orbital system. The bracket closes because valid inequalities refute its last four points, not because a facet list was computed. Cutting these points creates new intersections that a complete description would still have to enumerate and certify, so $\wedge^3\mathcal H_{11}$ remains open as a polytope.

Facetness and irredundancy of the four inequalities are likewise not addressed. They may well be facets; the certificates do not say so.

# Reproducibility

Every claim above maps to one artifact and one generating script in the accompanying data set, and each is independently checked.

| claim | artifact | check |
| --- | --- | --- |
| the four certificates | `level5_ressayre_3_11.json`, block `rank_eleven` | `verify_level5_ressayre_3_11_standalone.py` |
| complete calibration at $d=9,10$ | same, block `calibration` | linear programming against `constraints.json` |
| negative control and route split | same, block `calibration.negative_control` | recomputed in the test suite |
| randomized sweep | same, block `randomized_sweep` | seeded, regenerated by the script |
| the bracket closing at 19 and 31 | `bracket_3_11_settlement.json` | regenerated by `settle_bracket_3_11.py` |
| no design realizes index 44 | `rank11_open_candidate_designs.json` | exact counting, replayed in the test suite |

The standalone verifier imports nothing from the project package. The settlement generator applies the certified rows only after the four earlier methods, so it decides only what they left undecided, and the other forty-six verdicts retain their original certificates.

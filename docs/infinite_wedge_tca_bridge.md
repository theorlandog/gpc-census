# The infinite wedge, the stable semigroup, and why it does not compress the facets (2026-08-08)

Pre-registration: `docs/prereg_infinite_wedge_tca_sanity.md`, committed at
`cac70c5` before the measuring script was run. Artifacts:
`results/data/infinite_wedge_tca_sanity.json` (the measurement),
`results/data/infinite_wedge_tca_dependency.json` (the theorem-dependency
graph, machine readable). Generators: `scripts/infinite_wedge_tca_sanity.py`,
`scripts/infinite_wedge_tca_dependency.py`. Guard test:
`tests/test_infinite_wedge_tca_sanity.py`.

ONE PARAGRAPH ANSWER. All finite-rank fermionic moment polytopes for a fixed
particle number `N` are recoverable from one `GL_infinity` object, and the
recovery is completely explicit: the multiplicity of `V_lambda` in
`Sym^k(wedge^N C^d)` is the plethysm coefficient `<s_lambda, h_k[e_N]>`, which
does not mention `d`, so the finite-rank semigroups are the length truncations
`Sigma_{N,d} = Sigma_N^stab intersect {ell(lambda) <= d}` of a single
`d`-independent monoid, and the moment polytope is its normalized degree-one
slice. That is where the good news stops. The stable object is not the hard
part and never was; the hard part is that truncation by length is benign on
cones and hostile on FACETS. Every noetherianity theorem in the intended toolkit
(topological noetherianity of polynomial functors, bounded Plucker varieties and
the infinite wedge, FI/OI noetherianity, tca noetherianity) bounds EQUATIONS or
MODULES up to symmetry, and a facet of a moment cone is neither. Worse, the
covariant functoriality that any FI/OI argument would need does not exist: this
campaign refutes it exactly, by showing that a published GPC row generally
becomes INVALID at the next rank under zero padding, and also under every other
order-preserving embedding. The one direction that does work is contravariant
and is now proved: restriction to `lambda_{d+1} = 0` carries the rank `d+1`
system onto the rank `d` system exactly, at all six consecutive census pairs.

## 1. Exact finite-d semigroup convention

Fix `N` and `d`. Write `Lambda_d^+` for dominant integral weights of `GL_d`.
Define

    Sigma_{N,d} = { (lambda, k) in Lambda_d^+ x Z_{>0} :
                    mult(V_lambda, Sym^k(wedge^N C^d)) > 0 }.

CONVENTION, fixed and not assumed. `lambda` is a PARTITION with nonnegative
parts and `|lambda| = N k`, taken with NO dual, NO transpose and NO complement;
the normalization is `lambda / k`, read directly as the decreasingly ordered
occupation spectrum. Every part satisfies `lambda_1 <= k`, which is the
occupation bound `<= 1`.

✅ `Sigma_{N,d}` is a monoid. `Sym(wedge^N C^d)` is an integral domain and its
`U`-invariants form a subalgebra, so the product of a highest-weight vector of
weight `lambda` in degree `k` with one of weight `mu` in degree `l` is a nonzero
highest-weight vector of weight `lambda + mu` in degree `k + l`. Concretely
`Sigma_{N,d}` is the weight monoid of `Sym(wedge^N C^d)^{U_d}` graded by degree.

✅ `Sigma_{N,d}` is FINITELY GENERATED, by Hadziev-Grosshans: `C[X]^U` is
finitely generated for any affine `G`-variety `X` with `G` reductive and `U`
maximal unipotent. This is finiteness AT EACH FIXED `d` and carries no
uniformity whatsoever.

🟦 THE CONVENTION IS PINNED, 93 of 93, and the pin is not vacuous. Three
independently implemented routes to `mult(V_{k lambda}, Sym^k(wedge^N C^d))`
agree exactly at every one of the 93 (vertex, `k`) pairs that carry a measured
multiplicity sequence: the Murnaghan-Nakayama pairing `<s_{k lambda}, h_k[e_N]>`
in `gpc_census.plethysm`, the Weyl alternating sum over exact integer weight
counts in `scripts/quantization_multiplicity.py` (recomputed here at 45 of the
93 pairs, the rest inherited from the artifact that route generated), and the
stored integers of `results/data/quantization_multiplicities.json`. The control
matters more than the agreement: substituting the conjugate partition or the
reversed complement `(k - lambda_d, ..., k - lambda_1)` changes the answer at 77
of the 93 pairs, so the gate distinguishes conventions rather than passing
everything. [P1 PASS, P2 PASS]

## 2. Stable object and specialization

Let `A_N = Sym(wedge^N C^infinity)` as a polynomial functor, and

    Sigma_N^stab = { (lambda, k) : <s_lambda, h_k[e_N]> > 0 },

a sub-monoid of {partitions} x `Z_{>0}` in which `d` does not appear.

🟨 SPECIALIZATION LEMMA (proved here; elementary, and it is the whole content
of "one stable object"). For every `d`,

    Sigma_{N,d} = Sigma_N^stab intersect { ell(lambda) <= d }.

Proof. `Sym^k(wedge^N C^d) = sum_lambda S_lambda(C^d)^{a_lambda}` with
`a_lambda = <s_lambda, h_k[e_N]>` independent of `d`, and `S_lambda(C^d) = 0`
exactly when `ell(lambda) > d`. So membership at rank `d` is membership in the
stable monoid plus a length bound, nothing else. ∎

🟨 CONE COROLLARY. `{ell(lambda) <= d} = {lambda_{d+1} = 0}` is a FACE of the
dominant cone, and if a positive combination of partitions has `lambda_{d+1}=0`
then every summand does. Hence

    cone(Sigma_{N,d}) = cone(Sigma_N^stab) intersect { lambda_{d+1} = 0 },

so the whole family of moment cones consists of face-slices of ONE cone.

❗ THE POINT. Both statements are easy, and that is the finding. Every
difficulty the compression program hoped to import from stable representation
theory is already dissolved at this step: there is nothing subtle about the
stable object, its specializations, or the passage to cones. What survives is a
purely convex-geometric question, namely how the FACETS of a face-slice depend
on `d`, and no theorem in the intended toolkit speaks about that.

## 3. Moment-polytope recovery theorem

✅ STANDARD (Mumford, Ness, Brion; used in this project as Klyachko's practical
algorithm). The fermionic moment polytope is the normalized degree-one slice,

    Delta(N,d) = closure { lambda / k : (lambda, k) in Sigma_{N,d} },

and since `X = P(wedge^N C^d)` is the FULL projective space, every point of `X`
is an `N`-fermion pure state and no Zariski-closure subtlety about a proper
subvariety enters.

🟦 VERIFIED EXACTLY, and with both inclusions, not just the easy one.
[recovery block; P-level gate]

| system | degrees | generated points | every point satisfies every certified row | every certified vertex occurs | conv = Delta |
|---|---|---|---|---|---|
| (3,6) | m <= 6 | 20 | yes, 0 violations | 4 of 4 | **yes** |
| (3,7) | m <= 10 | 348 | yes, 0 violations | 10 of 10 | **yes** |

The argument closes because both halves are certified: every generated point is
an attainable spectrum AND satisfies every published row, so
`conv(points) subset Delta`; every certified vertex occurs among the points, so
`Delta subset conv(points)`. The second half leans on `P = O` from
`docs/census_inner_sandwich.md`, which is what makes the certified vertex list a
complete vertex list. The largest first degree needed at `(3,7)` is 8, so
`m <= 7` is NOT enough, which reproduces the note already in
`scripts/plethysm_inner_hull.py`.

❌ AND THE SEMIGROUP IS NOT RECOVERED. `Sigma_{N,d}` is NOT saturated, so the
polytope determines `cone(Sigma_{N,d})` and does NOT determine `Sigma_{N,d}`.
Four exact witnesses, all in the measured window, each a dominant lattice point
of the cone with multiplicity zero:

| vertex | lattice point `(lambda, k)` | multiplicity | quantization period |
|---|---|---|---|
| `(3,6)#3` | `((1,1,1,1,1,1), 2)` | **0** | 4 |
| `(3,7)#5` | `((2,2,1,1,1,1,1), 3)` | **0** | 6 |
| `(3,7)#7` | `((1,1,1,1,1,1,0), 2)` | **0** | 4 |
| `(3,7)#8` | `((2,2,2,2,2,1,1), 4)` | **0** | 8 |

Each `(lambda_int, q)` is in the cone because `lambda_int / q` is a certified
vertex of `Delta`, and it is an integral dominant weight; it is outside the
semigroup exactly when the quantization period exceeds `q`. So this is the same
`Z/2` obstruction that `docs/quantization_character.md` identified as a finite
isotropy character, read one level up as a saturation defect. The character side
predicts 54 such vertices census wide; ⚠️ those 54 are PREDICTED, not measured,
and the honest count of measured witnesses is 4 (3 independent, since `(3,7)#7`
is `(3,6)#3` padded). [P3 PASS, flagged ALREADY-KNOWN under R4: the numbers
were already in the artifact, only the semigroup reading is new]

CONSEQUENCE FOR THE PROGRAM, and it is a real one. "Finite generation of the
stable semigroup" and "finitely many facet orbits of the stable cone" are
DIFFERENT statements, and only the second is what a constructor needs. Any
compression claim must say which one it is about. The example above is the
proof that the distinction is not pedantic.

## 4. Candidate noetherianity theorems

⚠️ CITATION CAVEAT, recorded before the audit. The network policy of the
drafting environment blocked every primary host (arXiv, AMS, Springer, HAL,
De Gruyter), so each entry below was verified to ABSTRACT level only. A human
must confirm the precise statements before this section is quoted outside the
repository. The scope conclusions do not depend on fine print, and are stated so
that they survive any reasonable sharpening of the statements, but the
attributions do.

| theorem | object covered | finiteness proved | covers `A_N` | covers the weight semigroup | facet-orbit finiteness follows |
|---|---|---|---|---|---|
| Draisma, topological noetherianity of polynomial functors, JAMS 32 (2019) | `wedge^N C^infinity` as a polynomial functor | descending chains of `GL`-stable Zariski-closed subsets stabilize; equations up to orbit | yes | no | **no** |
| Draisma-Eggermont, Plucker varieties and higher secants of Sato's Grassmannian, Crelle 737 (2018); Nekrasov, dual infinite wedge equivariantly noetherian (2020) | bounded Plucker varieties, the infinite wedge | ideals generated in bounded degree up to `GL_infinity` | yes | no | **no** |
| FI/OI noetherianity: Church-Ellenberg-Farb; Aschenbrenner-Hillar; Hillar-Sullivant; Draisma, Noetherianity up to symmetry | finitely generated FI- and OI-modules | submodules of f.g. modules are f.g. | not directly | no | **no** |
| Nagpal-Sam-Snowden, degree two tcas, Selecta Math. 22 (2016) and sequel | unbounded tcas; known ONLY for `Sym(Sym^2)` and `Sym(wedge^2)` | module and ideal noetherianity | `A_2` yes, `A_N` for `N >= 3` **open** | no | **no** |
| Hadziev-Grosshans | `C[X]^U` for affine `G`-varieties | finite generation of the invariant algebra, hence of `Sigma_{N,d}` | yes, one `d` at a time | **yes** | no (no uniformity in `d`) |
| Vergne-Walter, J. Symplectic Geom. 15 (2017); Bulois-Denis-Ressayre (2025) | the moment cone of one FIXED representation, fermionic cone implemented explicitly | a finite, and in BDR minimal, inequality list at each fixed `d` | one `d` at a time | the cone, yes | no (no uniformity in `d`) |

❗ A COINCIDENCE WORTH RECORDING. The one fermionic case where tca noetherianity
IS known, `N = 2`, is exactly the case where the GPC answer is classically
trivial (`lambda_{2i-1} = lambda_{2i}`, one symmetric template of arity 2 at
every `d`). This project's main series, `N = 3`, sits precisely in the open
range. That is not evidence for anything, but it does mean the tca literature
has never been tested against a fermionic case where the answer is interesting.

## 5. Why each theorem does or does not imply facet compression

Every failure above has one of two causes, and both are structural rather than
technical.

**Cause 1: equations are not inequalities.** Topological noetherianity, bounded
Plucker varieties and the dual infinite wedge all bound the EQUATIONS of
Zariski-closed `GL`-stable sets, up to symmetry. A facet of a moment cone is a
real convex-geometric datum: an extreme ray of a dual cone. The family
`d -> Delta(N,d)` is not a descending chain of closed subsets of anything, and
"finitely many equations up to symmetry" carries no information about "finitely
many inequalities up to symmetry". This is the gap the mission asked to be made
explicit, and it is total, not partial. It is sharpened by an observation
specific to this problem: the GPC polytope is the moment image of the FULL
projective space `P(wedge^N C^d)`, not of a proper subvariety, so there is no
interesting Plucker variety in play at all. The equations these theorems bound
are equations of objects the GPC problem does not use.

**Cause 2: modules are not cones, and the module does not exist anyway.** FI/OI
and tca noetherianity bound submodules of finitely generated modules. Two
independent failures follow. First, module generation allows arbitrary
`R`-linear combinations while a cone needs NONNEGATIVE ones, so finite
generation of the span of the facet normals says nothing about the extreme rays
of the dual cone. Second, and fatally, the covariant action such an argument
needs does not exist, which this campaign establishes exactly rather than
suspecting.

❌ ZERO PADDING DOES NOT PRESERVE VALIDITY. Padding a published row with a zero
coefficient and testing it at every certified vertex of the next rank, which is
a test on the whole polytope because `P = O` is proved at both ranks:

| pair | rows still valid | rows tested | also literally published higher |
|---|---|---|---|
| (3,6) -> (3,7) | 0 | 4 | 0 |
| (3,7) -> (3,8) | 4 | 4 | 4 |
| (3,8) -> (3,9) | **10** | 31 | 10 |
| (3,9) -> (3,10) | **44** | 52 | 37 |
| (4,8) -> (4,9) | **8** | 15 | 8 |
| (4,9) -> (4,10) | **19** | 60 | 19 |

[P4 FAIL, and the FAIL is the finding.] The pre-registered expectation was PASS
above rank 6, on the strength of the disclosed pilot that the `(3,8)` table
visibly contains the four `(3,7)` rows padded. That pair is the exception, not
the rule. A named witness: `lambda_1 + lambda_8 <= 1` is a `(3,8)` row, and at
`(3,9)` the padded row is violated at the certified vertex
`(1,1/4,1/4,1/4,1/4,1/4,1/4,1/4,1/4)` with slack `-1/4`. The mechanism is plain
once seen: a published row names orbitals by POSITION relative to `d`, so
padding silently re-points it at a different orbital.

❌ AND NO OTHER ORDER-PRESERVING EMBEDDING SAVES IT. 🔬 POST-HOC block, added
after P4 failed and not pre-registered. For each row, all `d+1` order-preserving
injections `[d] -> [d+1]` were tried, which is exactly the OI action:

| pair | rows with SOME valid embedding | rows tested | higher rows that are OI images of a lower row |
|---|---|---|---|
| (3,6) -> (3,7) | 0 | 1 | 0 of 4 |
| (3,7) -> (3,8) | 4 | 4 | 4 of 31 |
| (3,8) -> (3,9) | 16 | 31 | 10 of 52 |
| (3,9) -> (3,10) | 47 | 52 | 39 of 93 |
| (4,8) -> (4,9) | 8 | 15 | 8 of 60 |
| (4,9) -> (4,10) | 35 | 60 | 19 of 125 |

So the inequality side of the census carries NO covariant OI structure: rows
with no valid embedding at all exist at five of the six pairs, and most of each
higher system is not an OI image of anything below it. There is no OI-module
here to which an OI-noetherianity theorem could be applied. That is the precise
sense in which the FI/OI route is not merely insufficient but inapplicable.

🟦 WHAT DOES WORK, AND IT IS A THEOREM. Restriction is complete downward.
`Delta(N,d) = Delta(N,d+1) intersect {lambda_{d+1} = 0}`, because an orbital of
occupation zero is annihilated by its own annihilation operator, so a state with
`lambda_{d+1} = 0` lies in `wedge^N C^d`; padding gives the converse. Hence
restricting the rank `d+1` system to `lambda_{d+1} = 0` cuts out exactly
`Delta(N,d)`. Verified by exact affine Farkas over the rationals at all six
consecutive census pairs, with every published lower row implied and none left
over. [P5 PASS, a forced gate on the tables, and passing it upgrades a citation
into a replayed proof.] Chaining the pairs: the published `(3,10)` system
implies the `(3,9)`, `(3,8)`, `(3,7)` and `(3,6)` systems, including the three
Borland-Dennis EQUALITIES, which are recovered as consequences rather than
assumed; and `(4,10)` implies `(4,9)` and `(4,8)`.

So the census inequality system is a purely CONTRAVARIANT functor on the rank
ladder. Compression down the ladder is free and now certified. Compression up
the ladder, which is the only direction that would help a constructor, is
refuted.

## 6. Smallest missing theorem

Stated no stronger than necessary, on the AMBIENT Weyl-invariant cone, because
`docs/screening_orbit_structure.md` already established that a chamber-facing
list is not the orbit-invariant object.

> **(T1) Bounded-arity facet templating.** Fix `N`. There are an integer
> `a = a(N)` and a finite set `T` of pairs `(c, b)`, with `c` a multiset of at
> most `a` nonzero rationals, such that for every `d` every facet of the ambient
> moment cone `C(N,d)` is, up to `S_d`, of the form
> `sum_j c_j lambda_{i_j} <= (b/N) sum_i lambda_i`.

> **(T2) Decidable instance rule.** There is an algorithm deciding, in time
> polynomial in `d`, whether a given instance of a template in `T` is a facet at
> rank `d`.

🟨 LEMMA (proved here). T1 implies that `C(N,d)` has `O(d^a)` facets, because an
`S_d`-orbit of an arity-`a` template has at most `binom(d,a) * a!` members.
Contrapositive: superpolynomial facet growth in `d` REFUTES T1. This converts
the compression question into a facet-counting growth question, which is the
smallest decisive experiment the program could run.

❗ T2, NOT T1, IS THE LOAD-BEARING HALF, and the measurement is what shows it.
The symmetric template stock is already far smaller than the row count:

| system | published rows | distinct symmetric templates | new templates vs all lower ranks | max symmetric arity |
|---|---|---|---|---|
| (3,6) | 1 | 1 | 1 | 3 |
| (3,7) | 4 | **1** | 1 | 4 |
| (3,8) | 31 | **9** | 8 | 7 |
| (3,9) | 52 | **8** | 6 | 7 |
| (3,10) | 93 | **16** | 9 | 8 |
| (4,8) | 15 | 3 | 3 | 4 |
| (4,9) | 60 | 13 | 11 | 7 |
| (4,10) | 125 | 15 | 10 | 9 |
| (5,10) | 161 | 21 | 21 | 9 |

All four `(3,7)` rows are one template, `lambda_a + lambda_b + lambda_c +
lambda_d <= 2`. Ninety-three `(3,10)` rows are 16 templates. So the ORBIT count
is small and the ORDERED count is large, which is the same shape the repository
already measured for flats (362 hyperplanes in 8 orbits at rank 6, 5,341 in 19
at rank 7). Knowing the templates therefore does NOT give a constructor; knowing
which ordered instances are facets does. That is T2.

❌ AND THE EVIDENCE AGAINST T1 IS ALREADY BAD. Two independent symptoms, both
measured, neither a proof:

1. The template stock does not saturate. New symmetric templates appear at every
   rank measured: 1, 1, 8, 6, 9 at `N = 3` and 3, 11, 10 at `N = 4`, cumulative
   1, 2, 10, 16, 25. [P6 PASS]
2. The maximum symmetric arity tracks the rank rather than a constant: 3, 4, 7,
   7, 8 at `d = 6..10` for `N = 3`, which is roughly `0.8 d` at the top, and 4,
   7, 9 for `N = 4`. It is not strictly increasing, since it repeats at 7 between
   ranks 8 and 9, so the pre-registered form FAILED. [P7 FAIL, and the FAIL is
   informative: arity plateaus locally while growing overall, so the sharp
   question is the trend, not monotonicity.]

If arity is genuinely unbounded then T1 is FALSE, and the whole templating idea
dies with it by the Lemma.

IS THAT THEOREM KNOWN OPEN OR FALSE? Neither, as stated for the fermionic cone:
it is not in the literature in this form. But the closest classical analogue is
decisive in the wrong direction. For the Horn cone of sums of Hermitian
matrices, the facet description IS completely known (Klyachko;
Knutson-Tao-Woodward), and the essential inequalities are indexed by subsets of
every size up to `n-1`, so their arity is UNBOUNDED in `n`. A moment cone whose
facets are completely understood fails bounded-arity templating. Absent a reason
why the fermionic cone should be better behaved than the Horn cone, T1 should be
expected to be false.

EASIER OR HARDER THAN OPEN TCA NOETHERIANITY? Logically INDEPENDENT, in both
directions. Section 5, Cause 2 shows tca noetherianity does not imply T1; and
T1, being a statement about extreme rays of a dual cone, says nothing about
ideals or modules, so it does not imply tca noetherianity either. Solving
`Sym(wedge^3 C^infinity)` noetherianity would therefore not advance this program
by one step. That is worth stating plainly, because it is the single most
likely misdirection available here: the open problem nearest to hand is not on
the path.

## 7. Finite sanity check

Scorecard, decided by `results/data/infinite_wedge_tca_sanity.json` alone.
Predictions are unedited from `docs/prereg_infinite_wedge_tca_sanity.md`
at `cac70c5`.

| prediction | verdict | scope | independent scope |
|---|---|---|---|
| P1 three routes to the multiplicity agree | PASS, 93 of 93 | 93 pairs | 13 vertices |
| P2 the convention is load bearing | PASS, 77 of 93 pairs discriminate | 93 pairs | 13 vertices |
| P3 the semigroup is not saturated | PASS, 4 witnesses | 14 vertices | 13 vertices |
| P4 zero padding preserves validity above rank 6 | **FAIL** | 162 rows | 5 pairs |
| P5 restriction is complete | PASS, 6 of 6 pairs | 169 rows | 6 pairs |
| P6 no template saturation in the window | PASS | 181 rows | 4 pairs |
| P7 symmetric arity strictly increases | **FAIL**, 3,4,7,7,8 | 181 rows | 5 systems |

P1 and P5 are forced gates on the code and the tables and are not quoted as
evidence for the science. P3 carries the ALREADY-KNOWN flag under R4: the
numbers were already in `results/data/quantization_multiplicities.json` and only
the semigroup reading is new. P4 is the campaign's result.

⚠️ SCOPE. Ranks 6 to 10, `N` in {3,4,5}, nine systems that are not independent of
each other (the wider ones reach the census partly through particle-hole duality
and lifts). Four consecutive ranks constrain no asymptotics. The multiplicity
side is permanently capped at rank 7 by the weight-count budget recorded in
`docs/symplectic_invariants.md`, which is why P1 to P3 rest on 14 vertices.

## 8. Rank-20 and rank-40 implication

WHAT T1 AND T2 WOULD BUY. Together they replace flat enumeration by template
instantiation: emit the `O(d^a)` instances of the finitely many templates and
decide each in polynomial time, instead of enumerating admissible hyperplanes.
The current constructor reaches 10,004,154 admissible hyperplanes in 231 orbits
at `(3,9)` and needed streaming before `(3,11)` fit in memory, so the comparison
is against a real and measured bottleneck.

🟦 AND THE ARITHMETIC SAYS IT IS NOT ENOUGH. Taking the 16 measured `(3,10)`
symmetric templates AS IF they were the complete list, which is an underestimate
by construction since new templates appear at every rank, the ordered instance
counts are:

| rank | candidate rows |
|---|---|
| 11 | 8,153,332 |
| 12 | 23,961,465 |
| 20 | 5,834,590,015 |
| 40 | 3,514,539,839,590 |

So even under the most favourable possible assumption, a template-instantiating
constructor faces about `8e6` candidates at rank 11, which is the same order as
the `1e7` hyperplanes the existing orbit enumerator already handles at rank 9;
about `6e9` at rank 20; and about `3.5e12` at rank 40. RANK 20 IS ARGUABLY
REACHABLE, RANK 40 IS NOT, and neither follows from T1 alone: the `d^a` factor
is the ordered instantiation, which is precisely what T2 would have to prune.

❗ AND THERE IS A CHEAPER ROUTE TO THE PHYSICS QUESTION, which the program should
decide about explicitly. Membership in a moment polytope is in NP and coNP
(Burgisser-Christandl-Mulmuley-Walter), with a polynomial-time weak membership
oracle by tensor scaling (Burgisser-Franks-Garg-Oliveira-Walter-Wigderson). If
the goal is "is this occupation spectrum admissible at rank 40", no facet list
is needed at all. If the goal is the exact facet system as a mathematical
object, the templating route is the one under discussion and the arithmetic
above applies. These are different projects and the compression program has been
conflating them. Relevant to the choice: Reuvers (J. Math. Phys. 62, 032204
(2021)) shows the GPC polytope volume approaches the Pauli hypersimplex volume as
`d` grows, so the generalized constraints are volumetrically negligible at large
`d` except in few-fermion settings, which weakens the physics case for rank 40
specifically.

## 9. Nonclaims

- ⚠️ Nothing here is a statement about all `d`. Every measurement lives in ranks
  6 to 10.
- ⚠️ T1 is NOT disproved. Two measured symptoms point against it and the Horn
  analogy points against it, and none of that is a proof. The Lemma says what
  would settle it: a superpolynomial lower bound on the facet count.
- ⚠️ The literature audit was verified to abstract level only, because the
  drafting environment blocked the primary hosts. The scope conclusions are
  robust to sharpening; the attributions need a human check before external use.
- ⚠️ P4's refutation is about the PUBLISHED rows. It shows that these particular
  valid rows do not transport up, which is enough to kill the naive covariant
  structure, and it is not a proof that no covariant structure on some other
  generating set exists.
- ⚠️ P5 proves that the restricted system cuts out the same SET. It proves
  nothing about facetness or irredundancy of any published row, both of which
  the repository already records as not established.
- ⚠️ The 54 census-wide non-saturation witnesses are PREDICTED, inherited from
  the character side, not measured. Four are measured.
- ⚠️ No facet of a rank-11 or higher system is computed here, and nothing here
  bears on candidate exhaustiveness at any rank.
- ⚠️ The full FI instantiation search, meaning all placements of a symmetric
  template's coefficients rather than the `d+1` order-preserving ones, was NOT
  run; it is combinatorially large at the measured arities. The refutation in
  section 5 is therefore stated for the OI action, which is the ordered
  embedding category the mission named.
- ⚠️ Nothing here computes a symplectic slice, a reduced space, a
  Duistermaat-Heckman density, or a fiber dimension. Every semigroup statement
  is about the moment polytope, which is FIXED-SPECTRUM data, and no stabilizer
  or fiber claim is made.

## MERGE_NOTES

For `docs/RESEARCH.md`, Part I, after the cross-realization block:

> ❌ THE STABLE-RANK COMPRESSION ROUTE IS REFUTED IN ITS COVARIANT FORM
> (2026-08-08). All finite-rank moment polytopes at fixed `N` are face-slices of
> one `d`-independent cone, by the elementary specialization
> `Sigma_{N,d} = Sigma_N^stab intersect {ell <= d}` with `Sigma_N^stab` the
> plethysm support. That step is easy and is not where the difficulty lives.
> ❌ A published GPC row generally becomes INVALID one rank up: zero padding
> keeps only 10 of 31 rows at (3,8) to (3,9), 44 of 52 at (3,9) to (3,10), 8 of
> 15 at (4,8) to (4,9), 19 of 60 at (4,9) to (4,10), and allowing ANY
> order-preserving embedding still leaves rows with no valid image at five of the
> six census pairs. So the inequality side carries no covariant FI/OI structure
> and no noetherianity theorem of that family can apply. 🟦 The contravariant
> direction is a theorem and is now replayed: restriction to `lambda_{d+1}=0` is
> complete at all six pairs by exact Farkas, so the published (3,10) system
> implies (3,9), (3,8), (3,7) and (3,6) including the Borland-Dennis equalities,
> and (4,10) implies (4,9) and (4,8). ❌ `Sigma_{N,d}` is NOT saturated: four
> exact witnesses, which are exactly the four period-`2q` vertices, so the
> quantization character is a saturation defect one level up. 🔬 The smallest
> missing theorem is bounded-arity facet templating plus a polynomial instance
> rule; measured arity tracks the rank (3,4,7,7,8 at `N=3`) and the template
> stock does not saturate, and the Horn cone, whose facets ARE known, fails
> bounded-arity templating, so the prior is bad. Even granting it, instantiating
> the 16 measured (3,10) templates costs 6e9 candidate rows at rank 20 and 3.5e12
> at rank 40.
> [docs/infinite_wedge_tca_bridge.md; docs/prereg_infinite_wedge_tca_sanity.md;
> results/data/infinite_wedge_tca_sanity.json;
> results/data/infinite_wedge_tca_dependency.json]

## Final answer

```text
FINITE-D SEMIGROUP:
  Sigma_{N,d} = {(lambda,k) : ell(lambda) <= d,
                 mult(V_lambda, Sym^k(wedge^N C^d)) > 0},
  the degree-graded weight monoid of Sym(wedge^N C^d)^{U_d}. Partition
  convention, no dual, no transpose, no complement; normalization lambda/k is
  the decreasing occupation spectrum. Pinned at 93 of 93 measured pairs by
  three independent routes, with 77 of 93 discriminating against the wrong
  conventions.

STABLE OBJECT:
  Sigma_N^stab = {(lambda,k) : <s_lambda, h_k[e_N]> > 0}, the plethysm support
  of A_N = Sym(wedge^N C^infinity). Specialization is literal truncation:
  Sigma_{N,d} = Sigma_N^stab intersect {ell(lambda) <= d}, and after passing to
  cones, cone(Sigma_{N,d}) = cone(Sigma_N^stab) intersect {lambda_{d+1} = 0},
  a face-slice of one d-independent cone.

KNOWN THEOREM THAT APPLIES:
  Hadziev-Grosshans (Sigma_{N,d} finitely generated, each fixed d);
  Mumford-Ness-Brion (the normalized slice IS the moment polytope, verified
  exactly at (3,6) and (3,7)); Vergne-Walter and Bulois-Denis-Ressayre (a
  finite, and minimal, facet list at each fixed d). All three are
  fixed-representation statements with no uniformity in d.

KNOWN THEOREMS THAT DO NOT APPLY:
  Draisma topological noetherianity of polynomial functors; Draisma-Eggermont
  bounded Plucker varieties and Nekrasov's dual infinite wedge; FI/OI
  noetherianity; Nagpal-Sam-Snowden tca noetherianity. Two causes: they bound
  equations of Zariski-closed sets or submodules of modules, and a facet is
  neither; and the covariant OI action their hypotheses need does not exist on
  the inequality side, which this campaign refutes by exact counterexample.

SMALLEST MISSING THEOREM:
  (T1) bounded-arity facet templating for the ambient fermionic moment cone,
  plus (T2) a polynomial-time rule deciding which instances are facets. T2 is
  the load-bearing half: the orbit count is already small (16 templates for 93
  rows at (3,10)) and the ordered instantiation is where the cost is.

IS THAT THEOREM KNOWN OPEN?:
  Not in the literature in this form, and logically independent of open tca
  noetherianity in both directions. The closest solved analogue, the Horn cone,
  FAILS bounded-arity templating, so T1 should be expected false. Measured
  symptoms agree: arity 3,4,7,7,8 over d = 6..10, template stock not
  saturating.

RANK-20/RANK-40 CONSEQUENCE:
  Granting T1 and T2 with the 16 measured (3,10) templates, which underestimates
  the template list, instantiation costs 8.2e6 candidate rows at rank 11,
  2.4e7 at rank 12, 5.8e9 at rank 20 and 3.5e12 at rank 40. Rank 20 is
  arguably reachable; rank 40 is not. For the physics membership question a
  polynomial-time weak membership oracle already exists and needs no facets.

GO/NO-GO:
  NO-GO on stable-algebra compression of the facet system. The stable object
  exists and is trivial; the noetherianity toolkit is inapplicable for
  structural reasons; the covariant functoriality it needs is refuted exactly;
  and the arithmetic denies rank 40 even if the missing theorem were granted.
  GO, already banked, on the contravariant direction: the highest-rank system
  implies every lower one, proved and replayed. CONDITIONAL GO on one narrow
  successor: measure facet count and maximum symmetric arity at (3,11) and
  (3,12) as the frontier reaches them. If arity keeps tracking the rank, T1 is
  dead and this line closes for good.
```

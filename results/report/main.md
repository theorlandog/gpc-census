---
bibliography: results/report/references.bib
csl: results/report/american-physics-society.csl
link-citations: true
nocite: "@*"
abstract: |
  We present an algorithmic pipeline that constructs exact extremal states for fermionic natural-occupation-number (moment) polytopes, and use it to complete the $\wedge^4\mathcal{H}_9$ polytope of Altunbulak and Klyachko \[Comm. Math. Phys. 282, 287 (2008)\]: the two vertices verified in that work "only numerically"---$(16,16,16,6,6,6,6,6,6)/21$ and $(20,14,14,14,14,4,4,4,4)/23$--- now have explicit extremal states with elementary, hand-checkable proofs, so the completeness argument for the published $\wedge^4\mathcal{H}_{10}$ constraint list is fully analytic. The pipeline decides, vertex by vertex, whether a one-hop-independent state (equivalently, a root-distinct support in the sense of Wildberger--Sjamaar and Maciazek--Tsanov) exists---by integer and mixed-integer programming with solver-independent certificates---and, when it does not, repairs minimal exchange-coupled ansatz classes by an exact phase-feasibility criterion, followed by verification in exact arithmetic. The resulting census shows all vertices of $\wedge^3\mathcal{H}_6$ and $\wedge^3\mathcal{H}_7$, and $27$ of $38$ vertices of $\wedge^3\mathcal{H}_8$ ($22$ integer at the natural denominator), are attainable by such states, while the second open $\wedge^4\mathcal{H}_9$ vertex is provably not: certified infeasibility, reproduced on two independent solvers, shows its extremal states cannot be one-hop independent in any orbital basis. The explicit state we construct for it carries a single forced phase, $\cos\gamma = 3/(4\sqrt{14})$; exhaustive enumeration excludes real amplitudes within its ansatz class, and a numerical orbit search indicates the state admits no real form under any orbital rotation (maximal antiunitary-symmetry overlap $\approx 0.9332$), while the existence of some other real extremal state for this vertex remains open. We describe how the same pipeline is being scaled toward the one-particle dimensions ($d\sim20$--$40$) required for quantum-chemical applications.
author:
  - James Orlando[^1]
date: |
  July 2026 --- DRAFT v0.7\
  Preprint: [doi:10.5281/zenodo.21313834](https://doi.org/10.5281/zenodo.21313834)
title: |
  Algorithmic construction of exact extremal states for fermionic moment\
  polytopes: completing $\wedge^4\mathcal{H}_9$
---

# Introduction

The Pauli exclusion principle bounds fermionic natural occupation numbers (NONs) by $0 \le \lambda_i \le 1$. Klyachko's solution of the one-body pure $N$-representability problem [@Klyachko2006; @AK2008] showed that the antisymmetry of $N$-fermion wave functions imposes far stronger conditions: the achievable ordered spectra $\lambda_1 \ge \cdots \ge \lambda_d$ of the one-body reduced density matrix (1-RDM) of a pure state in $\wedge^{N}\mathcal{H}_{d}$ form a convex polytope $\Pi_{N,d}$, cut out by finitely many _generalized Pauli constraints_ (GPCs). These constraints have measurable consequences: (quasi-)pinning of occupation spectra to facets of $\Pi_{N,d}$ constrains the structure of the wave function itself [@SGC2013; @Liebert2025], and the polytopes govern applications from reduced-density-matrix functional theory to quantum information [@Castillo2021; @Liebert2025].

Altunbulak and Klyachko computed $\Pi_{N,d}$ completely for all systems with $d \le 8$ (and for $\wedge^3\mathcal{H}_9$, $\wedge^3\mathcal{H}_{10}$, $\wedge^4\mathcal{H}_9$, $\wedge^4\mathcal{H}_{10}$, $\wedge^5\mathcal{H}_{10}$), supplying for each vertex of the small systems an explicit extremal state "for those who don't trust computer assisted proofs" [@AK2008]. For $\wedge^4\mathcal{H}_9$, however, two vertices resisted their representation-theoretic constructions: $$
v_A = \tfrac{1}{21}(16,16,16,6,6,6,6,6,6),
\qquad
v_B = \tfrac{1}{23}(20,14,14,14,14,4,4,4,4),$$ {#eq:AB} which were, in their words, "checked only numerically"---and on which the completeness of the published $\wedge^4\mathcal{H}_{10}$ constraint list depends [@AK2008]. To our knowledge these verifications have remained numerical for the intervening eighteen years.

#### Contributions.

1.  **An explicit extremal state for $v_A$** (Theorem [1](#thm:psiA)), with an elementary proof: seven Slater determinants, pairwise sharing at most two orbitals, with squared amplitudes $k_T/21$ for integer weights $k_T$ summing to $21$. The state would occupy a single line of Table 6 of Ref. [@AK2008].

2.  **A design/interference dichotomy.** Abstracting the structure of $\psi_A$, we define _weighted-design_ states (Definition [1](#def:design)) and show that design-attainability of a vertex is a pure integer-programming question. An exhaustive census (Section [-@sec:census]) over all vertices of $\Pi_{3,6}$, $\Pi_{3,7}$, $\Pi_{3,8}$ finds designs for $4/4$, $10/10$, and $22/38$ vertices respectively; for $v_B$ we prove by mixed-integer infeasibility that _no_ design exists, so every extremal state of $v_B$ requires phase cancellation between one-hop-connected configurations. The vertices that resisted the constructions of Ref. [@AK2008] are precisely of this interference type.

3.  **Canonical forms of GPC polytopes.** We compute $f$-vector data and the degrees of the adjoint hypersurfaces (equivalently, of the numerators of the Arkani-Hamed--Bai--Lam canonical forms [@ABL2018]) for all rigorously solved fermionic moment polytopes (Table [-@tbl:census]), extending to the pure-state GPC polytopes a bridge between $N$-representability and positive geometry known previously at the level of the hypersimplex [@Castillo2021; @LPW2020].

4.  **Methods.** An alternating-projection solver for "pure state with prescribed occupation spectrum," immune to the second-order degeneracy flatness that stalls gradient and moment-matching methods; a sparsifying concentration flow, deployed as a jittered multi-start swarm, which discovered $\psi_A$; and a fully independent verification pipeline, including a constructive re-verification of all $38$ published extremal states of $\Pi_{3,8}$.

# An explicit extremal state for $v_A$ {#sec:theorem}

We write $|i_1 i_2 i_3 i_4\rangle$, $i_1<i_2<i_3<i_4$, for Slater determinants in $\wedge^{4}\mathcal{H}_{9}$ over an orthonormal basis of $\mathcal{H}_9$, and call two determinants _one-hop connected_ if they share exactly three orbitals.

::: {#thm:psiA .theorem}
**Theorem 1**. \*The normalized state $$\psi_A \;=\; \frac{1}{\sqrt{21}}\Bigl(\sqrt{2}\,|1257\rangle + \sqrt{3}\,|1347\rangle + \sqrt{2}\,|1369\rangle - \sqrt{3}\,|1459\rangle + \sqrt{6}\,|1789\rangle + 2\,|2679\rangle + |3579\rangle \Bigr)$$ {#eq:psiA} has natural occupation numbers exactly $v_A = \tfrac1{21}(16,16,16,6,6,6,6,6,6)$, with the three heavy natural orbitals being the basis orbitals $\{1,7,9\}$. Hence the vertex $v_A$ of $\Pi_{4,9}$ is attained.\*
:::

::: proof
_Proof._ No two of the seven determinants in [@eq:psiA] are one-hop connected (each pair shares at most two orbitals), so every off-diagonal matrix element $\langle a_i^\dagger a_j\rangle_{\psi_A}$, $i\neq j$, vanishes identically: the 1-RDM is diagonal in the given basis, and the amplitudes' phases are pure gauge. The diagonal entries are the incidence sums $\lambda_m = \sum_{T \ni m} k_T/21$ with squared weights $k = (2,3,2,3,6,4,1)$ on the respective determinants. Orbital $1$ lies in the determinants of weight $2,3,2,3,6$, giving $16/21$; likewise orbitals $7$ and $9$; each remaining orbital collects weight $6$. Since $\sum_T k_T = 21$, the state is normalized, and its ordered spectrum is exactly $v_A$. $\square$
:::

::: remark
**Remark 1**. \*Relabeling orbitals so that the heavy modes are $\{1,2,3\}$ gives the equivalent canonical form $\psi_A = \tfrac1{\sqrt{21}}(\sqrt6|1239\rangle + \sqrt2|1247\rangle

- \sqrt3|1256\rangle + \sqrt2|1358\rangle + \sqrt3|1367\rangle
- 2|2348\rangle + |2357\rangle)$.\*
  :::

::: remark
**Remark 2**. _The proof is checkable by hand in minutes; the discovery was not. The state was found by the numerical pipeline of Section [-@sec:methods] and only then recognized to have exact structure. We regard the division of labor--- stochastic geometric search for discovery, elementary combinatorics for proof---as a portable pattern for moment-polytope problems._
:::

# Weighted designs and the dichotomy {#sec:census}

_Terminology and prior art._ The notion below is not new. A support in which no two determinants differ by a single orbital is, in representation-theoretic language, a _root-distinct_ set of weights: the weights of the supporting basis vectors differ pairwise by no root of the acting group. The idea is due to Wildberger [@W92] and was developed by Sjamaar [@S98]; in the quantum-marginal setting it appears as Definitions 3--4 of Maciazek and Tsanov [@MT17], stated there for spinful and spinless fermions, bosons and distinguishable particles alike, with applications to qubit systems [@M15] and to pinned fermionic occupations [@MS20]. Ref. [@MT17] further records that maximal root-distinct supports were already used to determine the $\wedge^3\mathcal{H}_6$ and $\wedge^3\mathcal{H}_7$ polytopes---with the support $\{123,145,246,356\}$ that our census recovers---so our rank-6 and rank-7 results are a rediscovery, reported here only for completeness. We retain the term _one-hop independent_ for the same condition, since our arguments are phrased in terms of single-orbital exchanges rather than roots. What is new here is the algorithmic layer: the exhaustive vertex-by-vertex decision procedure with solver-independent certificates, the exact repair step for vertices admitting no such state, and the resulting classification.

::: {#def:design .definition}
**Definition 1**. *Let $v = n/D$ be a vertex of $\Pi_{N,d}$ with integer numerator vector $n$ and $D = \tfrac1N\sum_m n_m$. A *weighted design* for $v$ is a family of $N$-subsets $\{T\}$ of $\{1,\dots,d\}$, pairwise sharing at most $N-2$ elements, together with positive real weights $w_T$ with $\sum_{T\ni m} w_T = n_{\sigma(m)}/D$ for some permutation $\sigma$. The associated state $\psi = \sum_T \sqrt{w_T}\,|T\rangle$ has exactly diagonal 1-RDM and ordered spectrum $v$; the design is *integer* if $w_T = k_T/D$ with $k_T \in \mathbb{Z}_{>0}$.*
:::

Design-attainability is thus a feasibility question for a system of integer incidence equations under an independence (anti-adjacency) condition---no quantum mechanics remains. We decided it exactly for every vertex of the rigorously solved small systems, by depth-first search with exact-cover pivoting and, where search was inconclusive, by mixed-integer linear programming with binary support variables and pairwise independence constraints (HiGHS; certificates of infeasibility retained).

::: {#prop:census .proposition}
**Proposition 1** (Census). _All $4$ vertices of $\Pi_{3,6}$ and all $10$ vertices of $\Pi_{3,7}$ admit integer weighted designs (reproducing, in the latter case, the extremal states of Ref. [@AK2008]). Of the $38$ vertices of $\Pi_{3,8}$, exactly $27$ admit weighted designs ($22$ of them integer at the natural denominator); the remaining $11$---the exotic vertex $(5,5,5,5,2,2,2,2)/28$ and its relatives---admit no weighted design with any positive real weights (mixed-integer certificates with continuous weight variables and binary support, each independently reproduced with a second solver of a different lineage, alongside feasibility controls on design-class vertices). The vertex $v_A$ of $\Pi_{4,9}$ admits the integer design of Theorem [1](#thm:psiA).\_
:::

::: {#thm:B .theorem}
**Theorem 2** ($v_B$ is interference-type). *No weighted design with any positive real weights exists for $v_B = (20,14,14,14,14,4,4,4,4)/23$. Since the family of Slater determinants and its one-hop structure are identical in every orthonormal orbital basis, the statement is basis-free: no extremal state of $v_B$ is one-hop independent in *any* orbital basis. Equivalently, in any basis in which the 1-RDM of an extremal state is diagonal, the supporting determinants necessarily fail one-hop independence; the vertex is attained through interference. (We thank T. Maciazek for prompting this basis-free formulation: a $2\times2$ rotation inside the mixing block of the state of Theorem [3](#thm:psiB) diagonalizes its 1-RDM while destroying independence, illustrating both halves of the statement.)*
:::

::: proof
_Proof (computer-assisted)._ A state with one-hop-independent support has exactly diagonal 1-RDM, so its squared amplitudes would give nonnegative real weights solving the incidence system on an independent support. The corresponding feasibility problem---$126$ continuous weight variables, $126$ binary support variables, $9$ incidence equalities, and $1260$ pairwise independence inequalities---is infeasible. The certificate was produced by the HiGHS solver, independently reproduced with the COIN-OR CBC solver, and is archived with the code; the same certificates were obtained for the integer-restricted problem at denominators $23$, $46$, and $69$. $\square$
:::

::: {#thm:psiB .theorem}
**Theorem 3** (Explicit extremal state for $v_B$). \*Let $\gamma = \arccos\bigl(3/(4\sqrt{14})\bigr)$. The normalized state $$
\psi_B = \tfrac{1}{\sqrt{23}}\Bigl(
2|1236\rangle + \sqrt2\,|1248\rangle + \sqrt7\,e^{i\gamma}|1249\rangle

- \sqrt3\,|1345\rangle + \sqrt2\,|1378\rangle + \sqrt2\,|1379\rangle
- |2358\rangle + \sqrt2\,|3489\rangle\Bigr)$$ {#eq:psiB} has natural occupation numbers exactly $v_B = \tfrac1{23}(20,14,14,14,14,4,4,4,4)$. Hence both numerically-verified vertices of Ref. [@AK2008] are attained exactly, completing the verification of the $\wedge^4\mathcal{H}_9$ moment polytope.\*
  :::

::: proof
_Proof._ The eight determinants pairwise share at most two orbitals except the two pairs $\{1248,1249\}$ and $\{1378,1379\}$, each differing by the exchange $8 \leftrightarrow 9$; hence every off-diagonal 1-RDM element vanishes identically except $\rho_{89}$. The diagonal entries are integer incidence sums of the squared weights $(4,2,7,3,2,2,1,2)$: orbitals $1$--$7$ receive $20,14,14,14,4,4,4$ respectively, and orbitals $8,9$ receive $7$ and $11$ (all divided by $23$). The exchange pairs contribute $\rho_{89} = \bigl(\sqrt{14}\,e^{i\gamma} + 2\bigr)/23
= (11 + i\sqrt{215})/92$ (both fermionic signs are $+1$), so $|\rho_{89}|^2 = (121+215)/92^2 = 21/23^2$. The remaining $2\times2$ block $\bigl(\begin{smallmatrix} 7 & z \\ \bar z & 11\end{smallmatrix}\bigr)/23$ has characteristic polynomial $\lambda^2 - 18\lambda + (77 - 21)
= (\lambda-14)(\lambda-4)$, giving eigenvalues $14/23$ and $4/23$ exactly. The ordered spectrum is $v_B$. $\square$
:::

::: {#thm:realexclusion .theorem}
**Theorem 4** (Real single-block exclusion). _No state of the form $\psi = \sum_T \varepsilon_T\sqrt{k_T/23}\,|T\rangle$ with $\varepsilon_T \in \{\pm1\}$, positive integers $k_T$, and one-hop connections confined to a single orbital pair attains $v_B$. (Exhaustive enumeration: $16$ ansatz equivalence classes, over $4.5\times10^6$ weighted supports, certified optimal terminations.) Consequently the complex phase in Theorem [3](#thm:psiB) is not an artifact: within the single-block class, $v_B$ is attainable only with nontrivial relative phase._
:::

::: {#rem:orbit .remark}
**Remark 3** (Non-realifiability under orbital rotations (numerical)). _Coefficient complexity is basis-dependent, so a sharper question is whether any orbital rotation makes the state of Theorem [3](#thm:psiB) real. A state $\psi$ admits a real form iff there exists $V \in U(9)$ with $\Lambda^4 V\,\bar\psi = e^{i\alpha}\psi$; such $V$ must map the eigenspaces of $\bar\rho$ to those of $\rho$, reducing the search to $U(1)\times U(4)\times U(4)$ aligned with the natural eigenspaces. Diagonal (phase-only) rotations are excluded analytically: the two exchange pairs force contradictory phase differences $2\gamma$ and $0$ on the same mode pair. Numerically, maximizing $|\langle\psi,\Lambda^4 V\bar\psi\rangle|$ over the full reduced group converges from independent random starts to the same value $\approx 0.933184 < 1$, indicating that $\psi_B$ admits no real form under any orbital rotation, with the deficit an apparent invariant of its orbit. This is a statement about the constructed state, not the vertex: whether $v_B$ admits some other, real extremal state remains open (exhaustive enumeration excludes integer-weight real single-exchange-block states; the continuous-weight real case on larger supports, and deeper block structures, remain open)._
:::

::: remark
**Remark 4** (Relation to the superselection rule of Ref. [@Liebert2025]). *Liebert *et al._ prove that saturation of a (spin-adapted) GPC restricts which configurations may contribute to the wave function. A vertex saturates many GPCs at once; the dichotomy sharpens the picture at this extreme point: for design vertices the intersected selection rules admit a diagonal, effectively classical weighted configuration mixture, while for interference vertices they force coherent superpositions of connected configurations. The historically difficult vertices of Ref. [@AK2008] are exactly the interference ones, which we propose as the structural reason they evaded sparse representation-theoretic constructions._
:::

# GPC polytopes as positive geometries {#sec:posgeo}

Every convex polytope carries a canonical form in the sense of Arkani-Hamed--Bai--Lam [@ABL2018]: a rational top-form with simple poles exactly on facets, whose numerator is the adjoint polynomial [@KohnRanestad]. A connection between $N$-representability and positive geometry was observed at the level of the hypersimplex [@Castillo2021; @LPW2020]; the pure-state GPC polytopes themselves appear not to have been examined in this light. For the rigorously solved systems the invariants are: $\Pi_{3,6}$ ($4$ vertices, $4$ facets, adjoint degree $0$; a simplex), $\Pi_{3,7}$ ($10$, $10$, $3$), $\Pi_{3,8}$ ($38$, $39$, $31$), and $\Pi_{4,8}$ ($22$, $23$, $15$), with adjoint degree $=F-\dim-1$. Table [-@tbl:census] reports the basic invariants for all rigorously solved systems, computed from exact rational vertex/facet enumeration (dual-volume evaluations of the canonical forms and numerical residue tests on facets---including the exotic facets of $\Pi_{3,8}$---are described in the repository).

Two structural observations from the same computations: positivity of the smallest occupation number is _implied_ by the GPCs at ranks $6$ and $7$ and returns as a genuine facet at rank $8$; and both rank-$8$ systems satisfy $F = V+1$, which we flag without explanation.

# Methods {#sec:methods}

#### Alternating projections for prescribed spectra.

Optimizing any functional of the eigenvalues of the 1-RDM stalls near degenerate targets: eigenvalues respond only at second order to intra-cluster perturbations, so matching moments to machine precision $\varepsilon$ constrains the spectrum only to $O(\sqrt{\varepsilon})$---which quantitatively reproduces the $\sim\!10^{-8}$ floors we observed for gradient, moment-matching (L-BFGS with verified analytic gradients), and support-function methods. The cure is to keep eigenvalues out of the inner loop: alternately (i) fix the natural orbitals $U$ of the current state and set the target matrix $T = U\,\mathrm{diag}(v)\,U^\dagger$, then (ii) minimize the smooth quartic $\|\rho(\psi) - T\|_F^2$ by L-BFGS with exact gradients. The resulting solver reaches prescribed spectra to $10^{-9}$ from random initializations, and its calibration on vertices with known extremal states gated every claim in this paper.

#### Concentration flow and swarm discovery.

On the solution manifold $\{\psi : \operatorname{spec}\rho(\psi) = v\}$ (dimension $\approx 242$ for $v_A$) we ran a sparsifying flow: ascent of the inverse participation ratio $\sum_T |c_T|^4$, alternated with full-power reprojection, guarded by rejecting any step whose reprojected spectral distance exceeds $10^{-8}$, with hard removal of sub-threshold amplitudes under the same guard. Deployed as a multi-start swarm (12 jittered walkers, consumer hardware), one walker cascaded from $126$ to $7$ determinants at spectral distance $2.3\times10^{-14}$ within minutes; the resulting support and rational squared amplitudes were then recognized and proved exactly (Theorem [1](#thm:psiA)). We note for honesty that single-trajectory runs of the same flow plateau near support $120$: conclusions about $242$-dimensional solution manifolds require swarm-scale exploration, and an early version of this project wrongly concluded from local searches that the extremal states were irreducibly dense.

#### Verification pipeline.

All claims passed: (i) independent re-derivation of the 1-RDM with separate fermionic-sign bookkeeping; (ii) exact rational arithmetic for the final proofs (no floating point); (iii) constructive re-verification of all $38$ published extremal states of $\Pi_{3,8}$ from Table 6 of Ref. [@AK2008], each reproducing its claimed vertex to machine precision; and (iv) calibration of every attainability method on vertices with known extremal states before its application to open ones.

# Scaling toward chemically relevant dimensions {#sec:scaling}

The gap this paper addresses at rank $9$ is an instance of a larger one: constraint tables exist only for $d \le 10$ (with $(3,d)$ systems reported to $d = 12$), while quantum-chemical applications need $d \sim 20$--$40$; the polytope program has largely paused at this frontier, in part because extremal-state construction has been manual. We outline a four-instrument program toward that regime, distinguishing carefully between rigorous and heuristic layers.

_(i) Exact inner structure (rigorous; this work, scaled)._ Root-distinct enumeration produces certified rational polytope points at any rank: the CP-SAT/MILP feasibility problems grow with $\binom{d}{N}$---$1140$ determinants at $(3,20)$, $9880$ at $(3,40)$---within modern solver reach, and every constructed state is an exact extreme-point candidate with a hand-checkable proof. The exchange-block repair step is a per-candidate polygon test whose cost is independent of $d$. This yields exact inner approximations (design hulls and their repairs) at all chemically relevant dimensions on commodity hardware.

_(ii) Rigorous boundary anchors from known ranks._ The polytope $\Delta(N,d)$ is precisely the face $\{\lambda_{d+1}=0\}$ of $\Delta(N,d+1)$, so the rigorously known tables at $d \le 12$ supply exact face data---boundary anchors---for every higher rank, and the closed-form Grassmann inequality families remain valid at all $d$. (We note that valid inequalities do _not_ in general lift from $\Delta(N,d)$ to $\Delta(N,d+1)$; the face relation, not naive lifting, is the correct rigorous statement.)

_(iii) Pointwise membership oracles (rigorous, pointwise)._ Membership in moment polytopes of reductive group representations---including the fermionic case explicitly---is in $\mathrm{NP}\cap\mathrm{coNP}$ [@BCMW17], and the tensor-scaling algorithms of Ref. [@BFGOWW18] provide an efficient weak membership oracle. Since the oracle iterates on state vectors of dimension $\binom{d}{N}$, pointwise boundary location (bisection along rays) is computationally cheap even at $(3,40)$. Our alternating-projection solver is a primitive relative of these methods; the scaling literature supplies convergence guarantees we lacked.

_(iv) Toward full outer certification._ The representation-theoretic computation of complete constraint lists becomes prohibitive beyond $d \approx 12$ along the original route. A recent algorithmic line [@vdBerg25] computes moment polytopes for general reductive-group representations via Franz's description, with a current frontier of $27$--$64$-dimensional representations (all of $\mathbb{C}^3{\otimes}\mathbb{C}^3{\otimes}\mathbb{C}^3$, and $\mathbb{C}^4{\otimes}\mathbb{C}^4{\otimes}\mathbb{C}^4$ with high probability). This does not yet reach the fermionic tables (already at representation dimensions $126$--$252$ in 2008), but it scales along a different complexity axis than the classical approach, and to our knowledge has not been applied to $\Lambda^N\mathbb{C}^d$. Validating it against the known fermionic tables, and attempting $(3,13)$---which would be the first new generalized Pauli constraint system since 2008---is a concrete experiment we highlight.

_Resources._ Instruments (i)--(iii) run at $(3,d)$ for $d \le 20$ on a single workstation in weeks; extending the inner atlas to $(3,40)$ is a modest cloud campaign (order $10^2$ vCPU-months). Instrument (iv) is gated on algorithm engineering rather than hardware. The binding constraint throughout is machine-readable vertex/facet data for the known systems (ranks $9$--$12$), for validation; we would welcome pointers to surviving copies of the original datasets.

# The census at ranks 9 and 10 {#sec:census910}

polytope $\dim$ vertices constraints design (int) design (real) interference

---

$\Pi_{3,6}$ 3 4 $1{+}3$eq all design (doubly-excited regime [@MT17])  
$\Pi_{3,7}$ 6 10 4 all design (doubly-excited regime [@MT17])  
$\Pi_{3,8}$ 7 38 31 26 1 11
$\Pi_{4,8}$ 7 22 15 22 0 0
$\Pi_{3,9}$ 8 58 52 37 1 20
$\Pi_{4,9}$ 8 103 60 85 2 16
$\Pi_{3,10}$ 9 113 93 69 2 42
$\Pi_{4,10}$ 9 159 125 132 2 25
$\Pi_{5,10}$ 9 292 161 246 4 42

: Every determinate fermionic moment polytope, with the design/interference census. Ranks $\le 8$ are rigorously complete; the rank-$10$ constraint lists are conjecturally complete as in Ref. [@AK2008], and all rank-$10$ rows are conditional on them. Verdicts transport under particle--hole duality, extending the classification to all $(N,d)$ with $d\le 10$. For $\Pi_{4,8}$ (half filling) the particle--hole involution maps vertex and facet sets to themselves; by uniqueness of the canonical form, $\Omega$ is exactly PH-invariant. Interference is absent at rank $8$ in the $N{=}4$ series and first appears at rank $9$---the system of this paper. The $16$ interference vertices of $\wedge^4\mathcal{H}_9$ (integer forms at natural denominator): $(12{,}9{,}5{,}5{,}5{,}3{,}3{,}3{,}3)$, $(10{,}7{,}7{,}4{,}4{,}4{,}2{,}1{,}1)$, $(10{,}7{,}5{,}5{,}5{,}2{,}2{,}2{,}2)$, $(9{,}6{,}6{,}4{,}4{,}4{,}1{,}1{,}1)$, $(18{,}12{,}12{,}7{,}7{,}4{,}4{,}4{,}4)$, $(9{,}6{,}5{,}5{,}5{,}2{,}2{,}1{,}1)$, $(15{,}9{,}8{,}8{,}8{,}3{,}3{,}3{,}3)$, $(7{,}4{,}4{,}4{,}4{,}2{,}1{,}1{,}1)$, $(16{,}9{,}9{,}9{,}9{,}4{,}4{,}2{,}2)$, $(18{,}10{,}10{,}10{,}10{,}4{,}4{,}3{,}3)$, $(28{,}15{,}15{,}15{,}15{,}6{,}6{,}6{,}6)$, $(20{,}12{,}12{,}12{,}12{,}4{,}4{,}4{,}4)$, $(14{,}11{,}11{,}6{,}6{,}6{,}2{,}2{,}2)$, $(14{,}9{,}9{,}9{,}9{,}3{,}3{,}2{,}2)$, $(11{,}8{,}7{,}7{,}7{,}2{,}2{,}2{,}2)$, and $v_B=(20{,}14{,}14{,}14{,}14{,}4{,}4{,}4{,}4)$. Real-not-integer vertices: $(6{,}3{,}3{,}3{,}3{,}3{,}1{,}1{,}1)$ and $(10{,}10{,}10{,}7{,}7{,}2{,}2{,}2{,}2)$. Full lists, verdicts, and validation code accompany the data bundle. {#tbl:census}

After the results above were obtained, we carried out the program proposed in the Discussion: extending the design/interference census to every rank-$9$ and rank-$10$ system. The machine-readable constraint data of Ref. [@AK2008] (hosted at Bilkent) is no longer accessible; we recovered the complete inequality lists for $\wedge^3\mathcal{H}_9$, $\wedge^4\mathcal{H}_9$, $\wedge^3\mathcal{H}_{10}$, $\wedge^4\mathcal{H}_{10}$ and $\wedge^5\mathcal{H}_{10}$ from the text of Altunbulak's thesis [@AltunbulakThesis], with parsed counts matching the published totals ($52$, $60$, $93$, $125$, $161$). Exact vertex enumeration (lrs, rational arithmetic) and the MILP census of Section 5 were then applied to all five systems, with every step cross-validated by independent structural invariants: embedding coherence ($\Delta(N,d)$ is the face $\{\lambda_{d+1}=0\}$ of $\Delta(N,d+1)$), frozen-core lifts $\lambda \mapsto (1,\lambda)$, particle--hole self-duality of $\wedge^5\mathcal{H}_{10}$ at the level of inequalities, vertices, and verdicts, and agreement with the published vertex tables at rank $8$.

(One inequality of the $\wedge^3\mathcal{H}_9$ table of Ref. [@AltunbulakThesis] is lost at a page boundary; we work with the repaired $52$-inequality system, which yields $58$ vertices consistent with all cross-checks. Details accompany the data bundle.)

**The census.** Table [-@tbl:census] summarizes the classification. For the system of this paper, $\wedge^4\mathcal{H}_9$: of $103$ vertices (matching the count of Ref. [@AK2008]), $85$ admit integer weighted designs at the natural denominator, $2$ admit real but not integer designs, and $16$ require interference. The vertex $v_A$ is classified DESIGN (consistent with Theorem 1) and $v_B$ INTERFERENCE (consistent with Theorem 2), each verdict landing blind inside a uniform sweep. Notably, $v_B$ is not isolated: the vertex $(20{:}12{:}12{:}12{:}12{:}4{:}4{:}4{:}4)/21$ shares its architecture (head $20$, degenerate quadruple, tail $4^4$) at the neighbouring natural denominator and likewise requires interference, so $v_B$ is the first-found member of a family. Across ranks $9\to 10$ the interference fraction is stable within each particle-number series ($34\to 37\%$ at $N{=}3$, $16\%$ at $N{=}4$, $14\%$ at $N{=}5$; the rank-$8$ systems sit lower, $29\%$ and $0\%$) while padding and frozen-core lifting generate the majority of each rank's interference vertices from lower-rank originals; at $\wedge^4\mathcal{H}_{10}$ they generate _all_ of them. Verdicts transport under particle--hole duality (the complement bijection on determinants maps designs to designs at the same denominator), so the classification extends to all systems $(N,d)$, $d \le 10$. The rank-$10$ constraint lists remain conjecturally complete as in Ref. [@AK2008]; all rank-$10$ statements here are conditional on those lists, which our tests probe but cannot fully certify.

# Discussion

The remaining open object is an explicit extremal state for $v_B$, now known to require interference; a natural attack is exact phase-cancellation solving over small connected supports, guided by the MILP census, or the facet conditions of Ressayre-type inequalities [@Ressayre]. The dichotomy itself invites theory: which incidence data of moment-polytope vertices are weighted-independent-set realizable is a clean combinatorial question, and the census suggests the answer encodes the boundary of "classically representable" vertices. The entire toolkit transfers to other pure-state quantum marginal settings---most directly the entanglement polytopes of Ref. [@WDGC2013]---where the same questions (explicit vertex states, design versus interference, canonical forms) appear open. Section [-@sec:census910] carries out this extension: the census now covers every rank-$9$ and rank-$10$ system, cross-validated by embedding, lifting, and duality invariants.

## Three questions this raises {#three-questions-this-raises .unnumbered}

\(i\) _Does the complex class exist, and persist?_ Theorem [3](#thm:psiB) exhibits the first genuinely complex extremal state known for these polytopes, Theorem [4](#thm:realexclusion) shows real amplitudes cannot be repaired within its mechanism class, and Remark [3](#rem:orbit) indicates the state admits no real form in any orbital basis; whether $v_B$ admits _any_ real extremal state is open, as is the same question for the rank-$9$ and rank-$10$ polytopes, whose interference vertices Section [-@sec:census910] now enumerates ($16$ at $\wedge^4\mathcal{H}_9$ alone, including a structural sibling of $v_B$); whether each admits a real extremal state, and whether field of definition is decidable vertex-by-vertex the way design-attainability is, remain open. We note that numerical verification cannot distinguish the classes: real-restricted optimization approaches $v_B$ to working precision even though no real state attains it. (ii) _Is the design/real/complex trichotomy an instance of a broader hierarchy?_ The three classes are naturally reminiscent of the positive/signed/complex gradation familiar from classical simulability and sign-problem considerations; whether that resemblance is structural or superficial we leave open. (iii) _What is the physics of a pinned phase?_ Occupation spectra are phase-blind, yet pinning to $v_B$ forces the relative phase $\cos\gamma = 3/(4\sqrt{14})$ exactly; the dynamical and metrological consequences of a geometrically pinned algebraic phase appear unexplored.

# Acknowledgments {#acknowledgments .unnumbered}

We are grateful to Tomasz Maciazek for generous and rapid correspondence on a first version of this manuscript: for pointing out the root-distinct literature and lineage now cited in the weighted-designs section, for the basis-free reformulation prompt absorbed into Theorem [2](#thm:B), and for emphasizing that the algorithmic construction, rather than any particular state, is the contribution most likely to be useful---advice that reshaped this draft. The computations, literature analysis, and drafting in this work were carried out in close collaboration with Claude (Anthropic); all results were verified through the independent pipeline described in Section [-@sec:methods], and responsibility for the claims rests with the author. The multi-start swarm computations were performed on the author's personal hardware.

# References {.unnumbered}

::: {#refs}
:::

[^1]: Independent researcher, Derry, NH, USA. `jamie@orlandonh.com`. ORCID: 0009-0008-3158-771X.

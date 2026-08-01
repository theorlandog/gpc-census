---
bibliography: results/report/references.bib
csl: results/report/physics-numeric.csl
link-citations: true
abstract: |
  The pure-state one-body fermionic marginal problem asks which ordered spectra occur as one-body reduced density matrices of states in an exterior power. A facet description of the resulting moment polytope does not itself provide pure-state preimages of its vertices. We give symbolic preimages for all 799 spectra appearing as vertices of nine tabulated rational descriptions with 6 through 10 orbitals. Two constructions replace numerical attainability checks left open for four fermions in nine orbitals. The first is root-distinct. The second uses two coherent exchange channels, and a finite rational certificate proves that no root-distinct preimage of that spectrum exists in any orthonormal orbital basis. The atlas contains 643 positive root-distinct witnesses, and 12 negative cases have finite no-root-distinct proofs. A separately implemented rational half-space calculation recovers all 799 tabulated vertices. Together with a row-by-row source concordance, these results complete the published four-fermion, nine-orbital description. The corresponding ten-orbital conclusion is conditional on validity of the tabulated inequalities.
keywords:
  - generalized Pauli constraints
  - quantum marginal problem
  - pure-state N-representability
  - moment polytopes
  - symbolic state preimages
author:
  - |
    James Orlando\
    Independent researcher, Derry, New Hampshire, USA\
    `jamie@orlandonh.com`
title: |
  Symbolic pure-state preimages for tabulated vertex spectra of fermionic one-body marginals through ten orbitals
---

# Introduction

For a normalized state $\psi\in\wedge^N\mathcal H_d$, define the one-body reduced density matrix (1-RDM) by
$$
\rho^{(1)}_{ij}(\psi)=\langle\psi|a_i^\dagger a_j|\psi\rangle,
\qquad \operatorname{Tr}\rho^{(1)}=N.\tag{1}
$$
Its ordered eigenvalues $\lambda_1\ge\cdots\ge\lambda_d$ are the natural occupation numbers. They obey Pauli's bounds $0\le\lambda_i\le1$, while antisymmetry imposes further linear inequalities. The three-fermion, six-orbital case was identified by Borland and Dennis [@BD1972]. Klyachko's general solution of the pure-state one-body $N$-representability problem identifies the attainable spectra with a rational moment polytope $\Pi_{N,d}$ [@Klyachko2006; @AK2008]. The additional facets are commonly called generalized Pauli constraints. They underlie work on pinning and quasipinning [@SGC2013; @Schilling2018; @Liebert2025], large-system structure [@Reuvers2021], and effective algorithms for moment polytopes [@Castillo2021; @vdBergSTOC2025].

An inequality description answers a membership question. It does not automatically provide a state with a prescribed admissible spectrum. At a vertex this distinction is especially useful: a pure-state preimage is a direct sufficiency witness, while a vertex obtained from a candidate outer description still requires an attainability argument. In this paper an *extremal state* means a normalized pure fermionic state whose 1-RDM spectrum is a vertex of the relevant polytope. It does not mean an extreme point of the convex set of many-body density operators.

Altunbulak and Klyachko supplied explicit states for the small systems in their Tables 5 and 6. In Sec. 6.2.2 of Ref. [@AK2008], however, two vertices of the 103-vertex description for $\wedge^4\mathcal H_9$ were checked only numerically:
$$
v_A=\frac1{21}(16,16,16,6,6,6,6,6,6),\qquad
v_B=\frac1{23}(20,14,14,14,14,4,4,4,4).\tag{2}
$$
The same source observed that attaining these spectra would finish the vertex argument for its published four-fermion, nine-orbital inequalities. Our first result replaces both numerical checks by short symbolic constructions.

The two states exhibit different mechanisms. If no two determinants in a support differ by one orbital, every off-diagonal 1-RDM contribution vanishes termwise. Such supports are root-distinct in moment-map language [@W92; @S98; @MT17]. They are also the fermionic analogue of free tensor supports, for which reduced density matrices are diagonal [@vdBergNonFree2025]. The state for $v_A$ has this form. The spectrum $v_B$ admits no root-distinct preimage after any unitary change of orbitals, but it has a sparse preimage with two coherently contributing exchange channels. Here and below, orbital bases are orthonormal and basis changes lie in $U(d)$; no claim is made about arbitrary $\mathrm{GL}(d)$ tensor transformations.

The second contribution is a uniform state atlas. It contains one symbolic determinant expansion for each of 799 tabulated vertices in nine systems with $6\le d\le10$, together with finite data that permit the states, classifications, and H-to-V calculations to be checked separately. A published algorithm computes moment polytopes from representation-theoretic input [@vdBergSTOC2025], while its expanded account develops the tensor and representation-theoretic setting in greater detail [@vdBergMomentPolytopes2025]. Our task is different: the rational H-descriptions are inputs, and the output is a preimage for each of their vertices.

The main conclusions are as follows.

1. Theorems [1](#thm:psiA) and [3](#thm:psiB) give hand-checkable states for $v_A$ and $v_B$.
2. Theorem [2](#thm:no-design) gives a finite proof that $v_B$ has no root-distinct preimage in any orthonormal orbital basis.
3. Theorem [4](#thm:atlas) gives symbolic preimages for all 799 tabulated vertices; Corollary [1](#cor:38) determines the root-distinct split for $\Pi_{3,8}$.
4. Proposition [3](#prop:exhaustion) reconstructs all nine tabulated V-representations from their rational H-representations.
5. Corollary [2](#cor:49complete) completes the published four-fermion, nine-orbital description. Corollary [3](#cor:complete) states the precise conditional result at ten orbitals.

Constraint generation beyond ten orbitals and fixed-spectrum fiber geometry are not addressed here.

# Root-distinct preimages and finite obstruction certificates {#sec:design}

Write $|T\rangle=|i_1\cdots i_N\rangle$ for the Slater determinant indexed by the $N$-subset $T=\{i_1<\cdots<i_N\}\subseteq[d]$. Two determinants are *one-hop connected* when $|T\cap U|=N-1$. A support is *one-hop independent* when it contains no connected pair. This is the root-distinct condition for the weights of the exterior-power representation.

::: {#lem:incidence .lemma}
**Lemma 1** (Incidence construction). *Let $S$ be a one-hop-independent family of $N$-subsets of $[d]$, and let $p_T>0$ with $\sum_{T\in S}p_T=1$. Then*
$$
\psi=\sum_{T\in S}\sqrt{p_T}e^{i\theta_T}|T\rangle,\tag{3}
$$
*has diagonal 1-RDM in the displayed orbital basis, with*
$$
\rho^{(1)}_{ii}=\sum_{T\in S:i\in T}p_T.\tag{4}
$$
*Consequently, any nonnegative solution of these incidence equations gives a pure-state preimage of the resulting occupation vector.*
:::

::: proof
*Proof.* A matrix element $\langle T|a_i^\dagger a_j|U\rangle$ with $i\ne j$ can be nonzero only when $T$ and $U$ differ by the exchange $j\leftrightarrow i$. This is excluded by the hypothesis. The diagonal formula follows from $a_i^\dagger a_i|T\rangle=\mathbf 1_{i\in T}|T\rangle$. $\square$
:::

For a rational target $v$, let $D$ be the least positive integer such that $Dv\in\mathbb Z^d$. A *weighted design* is a one-hop-independent support with positive rational weights whose incidence vector equals $v$. Requiring rational weights loses no real solution: a nonempty feasible set on a fixed support is a rational polyhedron and contains a rational point, after which zero weights may be removed. Designs related by a permutation preserving $v$ are identified.

Root-distinct attainability is a finite problem. Let $G_{N,d}$ be the graph whose vertices are the $N$-subsets of $[d]$ and whose edges join one-hop-connected determinants. Set
$$
a_T=(\mathbf 1_T,1)\in\mathbb Z^{d+1},\qquad
b=(v,1)\in\mathbb Q^{d+1}.\tag{5}
$$
A weighted design exists if and only if some independent set $S$ of $G_{N,d}$ supports a nonnegative solution of $\sum_{T\in S}p_Ta_T=b$.

For $y\in\mathbb Z^{d+1}$ with $y\mathbin{\cdot}b>0$, define
$$
P_y=\{T:y\mathbin{\cdot}a_T>0\}.\tag{6}
$$
Any nonnegative incidence solution must use a determinant in $P_y$. Thus each $y$ gives a hitting clause, and a finite branch proof may show that no independent set hits all clauses.

::: {#prop:certificate .proposition}
**Proposition 1** (Soundness of a no-design certificate). *Suppose a finite set $Y\subset\mathbb Z^{d+1}$ satisfies $y\mathbin{\cdot}b>0$ for all $y\in Y$, and a complete branch proof shows that no independent set of $G_{N,d}$ intersects every $P_y$. Then $v$ has no one-hop-independent pure-state preimage in any orthonormal orbital basis.*
:::

::: proof
*Proof.* If such a preimage existed, Lemma [1](#lem:incidence) would make its 1-RDM diagonal in that basis. Its diagonal is a permutation of $v$, and its squared coefficients give a nonnegative incidence solution on an independent support. After the corresponding orbital permutation, that support must hit every $P_y$, contrary to the branch proof. $\square$
:::

::: {#prop:certificate-complete .proposition}
**Proposition 2** (Existence of a finite certificate). *If no weighted design for $v$ exists, then there is a finite certificate of the form in Proposition [1](#prop:certificate).*
:::

::: proof
*Proof.* There are finitely many independent sets $S$ in $G_{N,d}$. For each $S$, infeasibility of $\sum_{T\in S}p_Ta_T=b$ with $p_T\ge0$ gives, by the rational Farkas lemma, a rational vector $y_S$ such that $y_S\mathbin{\cdot}b>0$ and $y_S\mathbin{\cdot}a_T\le0$ for every $T\in S$. Clearing denominators makes $y_S$ integral, and $S$ misses $P_{y_S}$. Taking these clauses over the finite family of independent sets gives a finite obstruction. A finite search tree that branches on an unhit clause supplies the branch proof. $\square$
:::

The supplementary certificates store integral Farkas vectors and a branch directed acyclic graph. At a node, the selected determinants are independent and a cited clause remains unhit. Its children cover the members of that clause that are not forbidden by one-hop conflicts, modulo permutations within equal-eigenvalue blocks. A dead leaf has no allowable member. Appendix [B](#app:no-design) gives the node invariant and a worked separator.

# The two four-fermion vertices {#sec:targets}

## A root-distinct state for $v_A$

::: {#thm:psiA .theorem}
**Theorem 1**. *The normalized state*
$$
\begin{aligned}
\psi_A=\frac1{\sqrt{21}}(&\sqrt2|1257\rangle+\sqrt3|1347\rangle
+\sqrt2|1369\rangle+\sqrt3|1459\rangle\\
&+\sqrt6|1789\rangle+2|2679\rangle+|3579\rangle)
\end{aligned},\tag{7}
$$
*has 1-RDM spectrum $v_A$.*
:::

::: proof
*Proof.* Every pair of supporting determinants shares at most two orbitals, so Lemma [1](#lem:incidence) makes the 1-RDM diagonal. The squared weights are $(2,3,2,3,6,4,1)/21$. Orbitals 1, 7, and 9 each collect weight $16/21$; every other orbital collects $6/21$. The weights sum to one, and sorting the diagonal gives $v_A$. $\square$
:::

## An obstruction for $v_B$

::: {#thm:no-design .theorem}
**Theorem 2**. *No normalized pure state with 1-RDM spectrum*
$$
v_B=\frac1{23}(20,14,14,14,14,4,4,4,4),\tag{8}
$$
*has one-hop-independent support in any orthonormal orbital basis.*
:::

::: proof
*Proof.* The certificate for $v_B$ contains 20 primitive integral separators. Closing their hitting sets under the $S_4\times S_4$ symmetry of the two fourfold-degenerate spectral blocks gives 2736 clauses. A 22-node branch directed acyclic graph of maximum depth 8 proves that no support can obey the 1260 one-hop exclusions while hitting every clause. Proposition [1](#prop:certificate) then applies. The finite objects are checked by direct integer arithmetic according to Appendix [B](#app:no-design). $\square$
:::

The theorem excludes a root-distinct expansion, not a sparse expansion of another kind. It does not say that complex phases or exactly two exchange channels are unavoidable.

## A coherent two-channel state for $v_B$

Define the algebraic unit phase
$$
\eta=\frac{3+i\sqrt{215}}{4\sqrt{14}},\qquad |\eta|=1.\tag{9}
$$

::: {#thm:psiB .theorem}
**Theorem 3**. *The normalized state*
$$
\begin{aligned}
\psi_B=\frac1{\sqrt{23}}(&2|1236\rangle+\sqrt2|1248\rangle
+\sqrt7\eta|1249\rangle-\sqrt3|1345\rangle\\
&+\sqrt2|1378\rangle+\sqrt2|1379\rangle
-|2358\rangle+\sqrt2|3489\rangle)
\end{aligned},\tag{10}
$$
*has 1-RDM spectrum $v_B$.*
:::

::: proof
*Proof.* All pairs of supporting determinants share at most two orbitals except $\{1248,1249\}$ and $\{1378,1379\}$, both associated with the exchange $8\leftrightarrow9$. With the convention $\rho^{(1)}_{ij}=\langle a_i^\dagger a_j\rangle$, their two contributions have the same sign and
$$
\rho^{(1)}_{89}=\frac{\sqrt{14}\eta+2}{23}
=\frac{11+i\sqrt{215}}{92}.\tag{11}
$$
The diagonal numerator is $(20,14,14,14,4,4,4,7,11)$. Put
$$
z=\sqrt{14}\eta+2=\frac{11+i\sqrt{215}}4,
\qquad |z|^2=21.\tag{12}
$$
The principal block on orbitals 8 and 9 is
$$
B=\frac1{23}\begin{pmatrix}7&z\\\bar z&11\end{pmatrix},\tag{13}
$$
and
$$
\det(xI-B)=\left(x-\frac{14}{23}\right)
\left(x-\frac4{23}\right).\tag{14}
$$
The other seven eigenvalues are the displayed diagonal entries, so the ordered spectrum is $v_B$. The squared amplitudes sum to $23/23$. $\square$
:::

The two exchange channels therefore add coherently to produce the required nonzero principal block. Theorem [2](#thm:no-design) explains why termwise diagonality cannot attain the same spectrum.

# State atlas through ten orbitals {#sec:atlas}

The census covers the nine systems in Table I. Particle-hole dual systems follow by the Hodge-star map and are not counted again. Each rational H-representation includes Pauli bounds, ordering, normalization, and its tabulated inequalities.

The rank-9 and rank-10 rows were transcribed from Appendix A of Ref. [@AltunbulakThesis]. A fresh ordered extraction matches all 60 rows at $(4,9)$ and all 379 rows at $d=10$ coefficient by coefficient. The appendix states that $(3,9)$ has 52 independent inequalities but prints 51. The final row in our normalized $(3,9)$ table equals printed $(3,10)$ row 82 restricted to $\lambda_{10}=0$; this identity does not establish that it was the source's intended omitted row. Appendix [D](#app:source-map) and the supplementary concordance give the page-level details.

Write $P^{\mathrm{tab}}_{N,d}$ for the rational polytope defined by the corresponding tabulated rows together with ordering, normalization, and Pauli bounds. This notation does not assume that a tabulated description equals $\Pi_{N,d}$.

\newpage

**Table I.** Census counts. A positive witness is a rational weighted design; a solver-negative label has no finite global obstruction in the supplement.

<!-- Alt text: Nine fermionic systems are compared by number of tabulated vertices, inequality rows, positive root-distinct witnesses, finite negative proofs, and solver-only negative labels. The totals are 799 vertices, 643 positive witnesses, 12 finite negative proofs, and 144 solver-only labels. -->

| system | vertices | tabulated rows | positive witness | finite no-root proof | solver-negative only |
|:---|---:|---:|---:|---:|---:|
| $P^{\mathrm{tab}}_{3,6}$ | 4 | 1 ineq. + 3 eqs. | 4 | 0 | 0 |
| $P^{\mathrm{tab}}_{3,7}$ | 10 | 4 | 10 | 0 | 0 |
| $P^{\mathrm{tab}}_{3,8}$ | 38 | 31 | 27 | 11 | 0 |
| $P^{\mathrm{tab}}_{4,8}$ | 22 | 15 | 22 | 0 | 0 |
| $P^{\mathrm{tab}}_{3,9}$ | 58 | 52 | 38 | 0 | 20 |
| $P^{\mathrm{tab}}_{4,9}$ | 103 | 60 | 87 | 1 | 15 |
| $P^{\mathrm{tab}}_{3,10}$ | 113 | 93 | 71 | 0 | 42 |
| $P^{\mathrm{tab}}_{4,10}$ | 159 | 125 | 134 | 0 | 25 |
| $P^{\mathrm{tab}}_{5,10}$ | 292 | 161 | 250 | 0 | 42 |
| **total** | **799** |  | **643** | **12** | **144** |

The 643 positive records comprise 631 witnesses on the target's least-denominator grid and 12 rational witnesses off that grid. Both types give symbolic incidence proofs. The assertion that an off-grid record has no on-grid witness is based on a finite integer search and is not used as a theorem. The search reported no weighted design for 156 targets. Twelve have the finite certificates of Proposition [1](#prop:certificate); the other 144 are solver classifications and are not rerun by the proof supplement. Symbolic 1-RDM reconstruction verifies the preimage in all 799 records, including the 156 states outside the positive design class.

::: {#cor:38 .theorem}
**Corollary 1** (Root-distinct split for three fermions in eight orbitals). *Exactly 27 of the 38 vertices of $\Pi_{3,8}$ admit a root-distinct preimage; the other 11 do not.*
:::

::: proof
*Proof.* The 27 positive records contain weighted-design witnesses, so Lemma [1](#lem:incidence) gives root-distinct preimages. The other 11 targets carry finite certificates satisfying Proposition [1](#prop:certificate). $\square$
:::

Altunbulak and Klyachko printed extremal-state formulas covering 70 system-vertex records represented in the present atlas: all 10 records at $(3,7)$, all 38 at $(3,8)$, and all 22 at $(4,8)$ [@AK2008]. The atlas supplies a common encoding and verification route for those states. With padding assigned before frozen-core descent to make the categories disjoint, 235 records descend by padding, 161 by a frozen-core lift, and 403 are primitive relative to those two operations. The novelty claim is restricted to exactifying the two cases reported numerically in the source; the other records are presented as a uniform symbolic atlas.

::: {#thm:atlas .theorem}
**Theorem 4** (Symbolic state atlas). *For every one of the 799 spectra in the nine supplementary vertex tables, the atlas contains a normalized symbolic state $\psi\in\wedge^N\mathcal H_d$ satisfying*
$$
\det(xI-\rho^{(1)}_\psi)=\prod_{i=1}^d(x-\lambda_i),\tag{15}
$$
*as an algebraic identity. The state records and vertex rows are in bijection.*
:::

::: proof
*Proof.* Each record specifies a determinant support and symbolic amplitudes built from rational numbers, radicals, and algebraic phases. Appendix [A](#app:state-check) gives the record convention and the direct reconstruction of every matrix element $\langle T|a_i^\dagger a_j|U\rangle$. Applied to each record, that calculation proves normalization and the displayed characteristic-polynomial identity. For every positive weighted-design record it also checks one-hop independence and the incidence weights. Matching the system, vertex index, canonical integer spectrum, and denominator establishes a bijection between the 799 records and the nine vertex tables. $\square$
:::

::: {#prop:exhaustion .proposition}
**Proposition 3** (Vertex exhaustion of the tabulated H-polytopes). *For each of the nine supplementary rational H-representations, the corresponding vertex table is exactly its vertex set. The nine sets contain 799 vertices in total.*
:::

::: proof
*Proof.* Start from the ordered Pauli simplex subject to the trace condition and any further affine equalities. Add the tabulated half-spaces one at a time. When a polytope $Q$ is cut by a closed half-space $H$, every vertex of $Q\cap H$ is either a retained vertex of $Q$ or the intersection of $\partial H$ with an edge of $Q$ that crosses $\partial H$. Two vertices of $Q$ form an edge precisely when the affine equalities and their common active inequality normals have rank $d-1$. The calculation in Appendix [C](#app:vertex-check) implements this induction with rational arithmetic. A final feasibility and active-rank check removes no point, and comparison with the nine tables has empty set difference in both directions. $\square$
:::

Theorem [4](#thm:atlas) and Proposition [3](#prop:exhaustion) are statements about the supplied algebraic states and rational coefficient tables. Identifying an H-representation with the true fermionic moment polytope additionally requires validity of its tabulated inequalities. That distinction matters at ten orbitals.

::: {#cor:49complete .theorem}
**Corollary 2** (Four fermions in nine orbitals). *The rational H-description published for $\wedge^4\mathcal H_9$ equals the pure-state moment polytope $\Pi_{4,9}$.*
:::

::: proof
*Proof.* The source inequalities are valid on $\Pi_{4,9}$ [@AK2008; @AltunbulakThesis], and the row concordance identifies all 60 tabulated rows in the supplementary table with the source. Let $P$ be their intersection with the Pauli, ordering, and trace conditions. Proposition [3](#prop:exhaustion) gives the complete 103-vertex set of $P$. Theorem [4](#thm:atlas), including the two formerly numerical cases proved in Theorems [1](#thm:psiA) and [3](#thm:psiB), gives a pure-state preimage for every vertex. Thus $V(P)\subseteq\Pi_{4,9}$, and convexity gives $P=\operatorname{conv}V(P)\subseteq\Pi_{4,9}$. The source inequalities give the reverse inclusion. $\square$
:::

::: {#cor:complete .theorem}
**Corollary 3** (Conditional ten-orbital conclusion). *Let $\mathcal S_{N,10}$ be one of the supplementary lists of 93, 125, or 161 tabulated inequalities, and define*
$$
P_{\mathcal S}=\left\{\lambda\in\mathbb R^{10}:
1\ge\lambda_1\ge\cdots\ge\lambda_{10}\ge0,
\ \sum_i\lambda_i=N,
\ A_{\mathcal S}\lambda\le b_{\mathcal S}\right\}.\tag{16}
$$
*If every inequality in $\mathcal S_{N,10}$ is valid on $\Pi_{N,10}$, then $P_{\mathcal S}=\Pi_{N,10}$.*
:::

::: proof
*Proof.* The hypothesis gives $\Pi_{N,10}\subseteq P_{\mathcal S}$. Proposition [3](#prop:exhaustion) lists every vertex of $P_{\mathcal S}$, and Theorem [4](#thm:atlas) supplies a pure-state preimage for each of the 564 listed vertices. Hence $V(P_{\mathcal S})\subseteq\Pi_{N,10}$. Convexity gives the reverse inclusion $P_{\mathcal S}=\operatorname{conv}V(P_{\mathcal S})\subseteq\Pi_{N,10}$. $\square$
:::

The source concordance verifies transcription, not representation-theoretic validity. No independent derivation of every ten-orbital inequality is claimed.

# Construction and verification {#sec:methods}

## State construction

For a positive weighted-design record, the state is built directly from its incidence witness,
$$
\psi=\sum_{T\in S}\sqrt{p_T}\,|T\rangle.\tag{17}
$$
Lemma [1](#lem:incidence) then proves the target spectrum without numerical optimization. Padding adds an empty orbital and occupation number 0; a frozen-core lift adds an always-occupied orbital and occupation number 1; particle-hole duality acts by the Hodge star. Transported states are checked again at their destination.

For targets not reached by a root-distinct support, sparse numerical states served as discovery data. Their amplitudes were subsequently recognized as rational or algebraic expressions, and only the symbolic expressions enter Theorem [4](#thm:atlas). Discovery convergence is therefore not a premise of the theorem.

The least-denominator integer search used OR-Tools CP-SAT [@PerronDidierGay2023], and the real-weight feasibility search used COIN-OR CBC [@ForrestLougeeHeimer2005]. SymPy performs algebraic simplification and characteristic-polynomial calculations [@Meurer2017]. The original vertex tables were generated with `lrs` [@Avis2000LRS]; Proposition [3](#prop:exhaustion) uses a separate implementation of incremental double description [@FukudaProdon1996].

## Computer-assisted components

\newpage

**Table II.** Computer-assisted components and their mathematical status.

<!-- Alt text: Nine claim classes are distinguished as symbolic identities, direct constructions, finite computer-assisted proofs, source concordances, solver-only classifications, or conditional theorems. -->

| Claim | Status | Basis |
|:---|:---:|:---|
| all 799 preimages | symbolic, computer-checked | determinant expansions and characteristic polynomials |
| all 643 root-distinct witnesses | direct | rational incidence weights |
| 12 negative labels | finite computer-assisted proof | branch certificates |
| 144 negative labels | solver classification | exploratory CP-SAT and CBC outcomes; not used in formal results |
| 490 printed rank-9 and rank-10 rows | source-checked | ordered Appendix A concordance |
| final $(3,9)$ table row | derived identity | restriction of printed $(3,10)$ row 82 |
| all 799 tabulated vertices | finite computer-assisted proof | rational H-to-V reconstruction |
| completeness at $(4,9)$ | theorem | Corollary [2](#cor:49complete) |
| completeness at $d=10$ | conditional theorem | Corollary [3](#cor:complete) |

The proof supplement supports the symbolic preimages, 12 no-design certificates, source concordance, and rational H-to-V calculation. The 144 solver classifications record exploratory CP-SAT and CBC outcomes and are not part of that proof supplement. None enters Theorem [4](#thm:atlas), Proposition [3](#prop:exhaustion), or Corollaries [2](#cor:49complete) and [3](#cor:complete).

**Generative-AI assistance.** OpenAI Codex in ChatGPT Work (GPT-5.6, OpenAI, accessed July--August 2026) and Anthropic Claude (Claude Opus 5, Anthropic, accessed July 2026) assisted with exploratory-code development and inspection, test-case generation, and language editing because the computational study contains 799 state records. The author reviewed the code, proofs, data, and text and is responsible for the work.

# Discussion

The two historical gaps in $\wedge^4\mathcal H_9$ have different resolutions. The first required a missing combinatorial construction: $v_A$ has a seven-term root-distinct preimage. The second lies outside that mechanism: no unitary orbital basis makes a $v_B$ preimage root-distinct, while an eight-determinant state attains it through a nonzero coherence built from two exchange channels. This is a unitary-basis obstruction attached to a marginal-polytope vertex. It is distinct from non-freeness under general local $\mathrm{GL}$ transformations [@vdBergNonFree2025].

The atlas turns isolated existence checks into a common algebraic test set. It can be used to test symbolic 1-RDM calculations, state-construction heuristics, and methods designed for pinned occupations. These are prospective uses; no many-body solver, density-matrix functional, or state-preparation protocol is benchmarked here.

The principal limitation concerns the negative root-distinct classification. State verification covers all 799 records, but global no-design certificates cover 12 of the 156 negative labels. The other 144 remain solver classifications. Extending the finite certificate collection would strengthen that taxonomy but is not needed for the preimage theorem or for the four-fermion, nine-orbital completion. At ten orbitals, the remaining limitation is different: the coefficient tables agree with their source, but completeness remains conditional on the mathematical validity of those inequalities.

# Conclusion

We have replaced two numerical preimages in the four-fermion, nine-orbital analysis by symbolic states. One is root-distinct; the other has a finite obstruction to every root-distinct realization and is attained by two coherent exchange channels. The accompanying atlas supplies symbolic preimages for all 799 vertices of nine tabulated rational H-descriptions, while a rational half-space calculation reconstructs the complete vertex lists. A source concordance closes the provenance of the four-fermion, nine-orbital rows, yielding an unconditional completion of that published description. The ten-orbital implication is stated with its remaining validity hypothesis.

# Supplementary Material {.unnumbered}

See the Supplementary Material for the state-record schema, proof-object formats, source-concordance metadata, archive manifest, and reproduction instructions. The accompanying machine-readable research-data archive contains the 799 symbolic states, nine rational H- and V-tables, 12 no-design certificates, source concordance, and checking programs.

# Author declarations {.unnumbered}

## Conflict of interest {.unnumbered}

The author has no conflicts to disclose.

## Author contributions {.unnumbered}

James Orlando: Conceptualization; Methodology; Software; Validation; Formal analysis; Investigation; Data curation; Writing - original draft; Writing - review and editing; Project administration.

## Data availability {.unnumbered}

The Supplementary Material and machine-readable research-data archive are submitted with this manuscript. Versioned project materials are also available from Zenodo [@OrlandoZenodo2026] and the [public source repository](https://github.com/theorlandog/gpc-census).

# Appendix A: State records and 1-RDM reconstruction {#app:state-check .unnumbered}

Formulas in the manuscript number orbitals from 1; the supplementary determinant arrays number them from 0. Each state record is keyed by its system and vertex index and stores the canonical integer spectrum $n$, its least denominator $D$, determinant arrays, and symbolic amplitudes. The Supplementary Material gives the complete schema.

For a sorted determinant $U=(u_1<\cdots<u_N)$, the action of $a_j$ is zero unless $j=u_p$, in which case removal contributes $(-1)^{p-1}$. If $i$ is absent from the remaining determinant, inserting it in sorted position $q$ contributes $(-1)^{q-1}$. Thus the verifier evaluates
$$
\rho^{(1)}_{ij}=\sum_{T,U}\overline{c_T}c_U
\langle T|a_i^\dagger a_j|U\rangle,\tag{A1}
$$
directly from the serialized determinant arrays. It checks $\sum_T|c_T|^2=1$ and compares every coefficient of
$$
\det(xI-\rho^{(1)})-\prod_{k=1}^d\left(x-\frac{n_k}{D}\right),\tag{A2}
$$
as an algebraic number. The cross-link checker separately confirms that each `(system, index, n, D)` agrees with exactly one vertex-table row. Positive weighted-design records receive the additional one-hop and incidence checks of Lemma [1](#lem:incidence).

# Appendix B: No-design certificates {#app:no-design .unnumbered}

For each integral separator $y$, the checker recomputes $y\mathbin{\cdot}b$ and the clause $P_y$ from the determinant incidence columns. As a simple example from the $v_B$ certificate, take
$$
y=(2,2,2,2,2,-3,-3,-3,-3;-3).\tag{B1}
$$
The last entry pairs with the normalization coordinate. Direct calculation gives $y\mathbin{\cdot}(v_B,1)=35/23>0$. For a four-subset $T$, the inequality $y\mathbin{\cdot}a_T>0$ holds precisely when all four orbitals of $T$ lie among the first five. Hence $P_y$ is the five-element family of four-subsets of $\{1,2,3,4,5\}$. Any weighted design must select at least one of them.

At every branch node the certificate records an independent selected set $S$ and an unhit clause $P_y$. For every $T\in P_y$ that is compatible with $S$, a child records $S\cup\{T\}$, with orbit-equivalent children represented once under permutations inside equal-eigenvalue blocks. The checker expands those orbits, verifies that their union is exactly the compatible part of $P_y$, and verifies the child invariant. A leaf is valid only when that compatible part is empty. Acyclicity, reachability, and complete child coverage prove that no independent support hits all clauses.

# Appendix C: Rational H-to-V reconstruction {#app:vertex-check .unnumbered}

Let $E$ contain the affine equalities and let $R$ be the inequalities already imposed. For a vertex $u$, write $A(u)\subseteq R$ for its active inequalities. The update used in Proposition [3](#prop:exhaustion) is:

1. Start with the vertices of the ordered Pauli simplex in the affine space $E$.
2. For the next half-space $h(x)\le c$, retain every old vertex satisfying it.
3. For each pair $u,v$ on opposite sides, test
   $$
   \operatorname{rank}\bigl(E\cup A(u)\cap A(v)\bigr)=d-1.\tag{C1}
   $$
   If the equality holds, add the unique intersection of the segment $[u,v]$ with $h(x)=c$.
4. Merge coincident points, recompute their complete active sets, and continue.

The rank condition is the standard adjacency criterion in the affine hull. The vertices of a half-space cut are exactly the retained vertices and intersections of the cutting hyperplane with crossed edges, so induction proves that every intermediate vertex set is complete. The implementation uses `Fraction` arithmetic for solves, slacks, intersections, and ranks. After the final cut it verifies every inequality and checks that the equality normals plus active inequality normals have full rank $d$ at each point. Exact set comparison with the nine supplied tables then proves Proposition [3](#prop:exhaustion).

# Appendix D: Constraint-source concordance {#app:source-map .unnumbered}

The source checker identifies the dissertation PDF by SHA-256, extracts Appendix A with layout-preserving text conversion, and parses each displayed $\lambda_i$ coefficient and right-hand side. It compares the rows in order with both the transcription and the normalized H-table. All 490 printed rows match: 51 at $(3,9)$, 60 at $(4,9)$, and 93, 125, and 161 at $(3,10)$, $(4,10)$, and $(5,10)$. The unprinted 52nd row in our $(3,9)$ table equals the restriction of printed $(3,10)$ row 82. This establishes the algebraic identity but not the historical intent of the omitted source row. The PDF itself is not redistributed.

# References {.unnumbered}

::: {#refs}
:::

# Prior Art and Cross-Field Roadmap for the Inverse-Fiber Program (2026-07)

Purpose: every object in RESEARCH_NEXT.md has a mature sibling in some other
field. This file names the sibling, cites the canonical sources, and states
the concrete import. Ordered by expected payoff. GPC-internal prior art
(pinning selection rules, UDA/UDP, entanglement polytopes, stratified
symplectic geometry) is covered in the paper's related-work section and the
2026-07 literature review; this file is the OUTSIDE-GPC map.

---

## 1. Phase retrieval: the methodological twin

THE MATCH. The coherence equations C_S(p,θ) = 0 are a phaseless-measurement
recovery problem: quadratic magnitude/interference data, unknown phases,
gauge group U(1)^d × U(1), conjugation ambiguity. The census's time-reversal
sheets (±Θ at v_B) are phase retrieval's conjugate-reflection ambiguity; the
"complex representative vs real fiber point" confusion the closure resolved
is a standard PR phenomenon.

Sources:
- A. Conca, D. Edidin, M. Hering, C. Vinzant, "An algebraic characterization
  of injectivity in phase retrieval," Appl. Comput. Harmon. Anal. 38 (2015)
  346–356. [Generic frames of 4n−4 vectors give injectivity in C^n; algebraic
  variety methods for uniqueness thresholds.]
- P. Grohs, S. Koppensteiner, M. Rathmair, "Phase Retrieval: Uniqueness and
  Stability," SIAM Review 62 (2020) 301–350. [Survey; note it reaches back to
  the PAULI PROBLEM literature, wavefunction determination from marginal
  distributions, so the bridge to quantum marginals is already in print.]
- D. Edidin, "The geometry of ambiguity in one-dimensional phase retrieval,"
  SIAM J. Appl. Algebra Geom. 3 (2019) 644–660. [A worked example of
  "geometry of the ambiguity set" = fiber geometry, exactly the program's
  genre.]
- T. Bendory, D. Edidin et al., "Toward a mathematical theory of the
  crystallographic phase retrieval problem," SIAM J. Math. Data Sci. (2020);
  arXiv:2002.10081. [SPARSE phase retrieval: uniqueness from sparsity, the
  direct analog of sparse realization / the defect program.]

IMPORTS:
(a) Injectivity thresholds as a function of measurement count ↔ fiber
    triviality as a function of support size: state and prove a census
    analog of the 4n−4 theorem for support-restricted spectral maps.
(b) The ambiguity-group formalism: classify the census fiber's discrete
    symmetries (conjugation, sign flips) before computing moduli, as PR does.
(c) Crystallographic PR's sparsity-uniqueness results are the template for
    Problems 5–6 (infinite defect / rigid families).

## 2. Rigidity theory of frameworks: the silence-rank lemma's home field

THE MATCH. Determinants = joints, channels = bars, the silence Jacobian =
the rigidity matrix, first-order fixed-spectrum rigidity = infinitesimal
rigidity, and the open "is within-block-rank = kernel-dim forced?" question
is a GENERICITY question, the exact kind rigidity theory answers.

Sources:
- L. Asimow, B. Roth, "The rigidity of graphs," Trans. AMS 245 (1978)
  279–289. [Infinitesimal vs continuous rigidity; generic rigidity is a
  property of the graph alone.]
- R. Connelly, "Generic global rigidity," Discrete Comput. Geom. 33 (2005)
  549–563; S. Gortler, A. Healy, D. Thurston, "Characterizing generic
  global rigidity," Amer. J. Math. 132 (2010). [Stress matrices certify
  GLOBAL rigidity, the analog of proving a fiber point is the UNIQUE state
  over its support, not just isolated.]
- J. Graver, B. Servatius, H. Servatius, "Combinatorial Rigidity," AMS
  Grad. Studies 2 (1993). [The rigidity matroid.]

IMPORTS:
(a) Define the COHERENCE RIGIDITY MATROID on the channel set; the lemma's
    rank conditions become matroid statements, and "forced vs accidental"
    (sample size 2) becomes a generic-rank theorem or a matroid
    counterexample search.
(b) Stress-matrix certificates: a second-order tool for upgrading
    first-order rigidity (v89) to global-over-support rigidity, and for the
    germ-smoothness gap in Remark 2 of the lemma.

## 3. Structured inverse eigenvalue problems: the realization problem is a PEIEP

THE MATCH. "Find a Hermitian matrix with prescribed spectrum and prescribed
zero pattern (plus the ψ-realizability structure)" is a Prescribed-Entries
Inverse Eigenvalue Problem, a named problem class with fifty years of
solvability theory and algorithms. The fixed-spectrum fiber is an
ISOSPECTRAL MANIFOLD intersected with the support constraint; Chu's
isospectral flows move along exactly this manifold.

Sources:
- M.T. Chu, G.H. Golub, "Structured inverse eigenvalue problems," Acta
  Numerica 11 (2002) 1–71.
- M.T. Chu, G.H. Golub, "Inverse Eigenvalue Problems: Theory, Algorithms,
  and Applications," Oxford UP (2005). [Chapters on Partially Described
  IEPs and Structured Low-Rank Approximation are the closest.]
- M.T. Chu, "Inverse eigenvalue problems," SIAM Review 40 (1998) 1–39.

IMPORTS:
(a) Isospectral flow algorithms (double-bracket / projected gradient on the
    isospectral orbit) as production fiber explorers, a principled
    replacement for the ad-hoc continuation used in the sparsification
    cascade.
(b) Solvability results for zero-pattern IEPs: candidate exact statements
    for which supports admit ANY matrix with the target spectrum, a
    matrix-level relaxation of quantum feasibility sitting between s_I and
    s_Q, i.e. a NEW rung for the defect ladder (call it s_H, the Hermitian
    pattern bound: matrix exists but maybe no ψ).

## 4. Numerical algebraic geometry: the professionalized version of this session

THE MATCH. Witness sets, monodromy, trace tests, and α-certification are the
industrial-strength versions of what the fiber probes did by hand, plus the
one thing the probes lacked: RIGOROUS certification of numerical points.

Sources:
- A.J. Sommese, C.W. Wampler, "The Numerical Solution of Systems of
  Polynomials Arising in Engineering and Science," World Scientific (2005).
- J.D. Hauenstein, F. Sottile, "Algorithm 921: alphaCertified: certifying
  solutions to polynomial systems," ACM TOMS 38 (2012) Art. 28. [Smale
  α-theory in exact rational arithmetic: turns "residual 1e-16" into a
  THEOREM that a true solution exists nearby.]
- A. Leykin, J.I. Rodriguez, F. Sottile, "Trace test," Arnold Math. J. 4
  (2018) 113–125. [Certifies completeness of a witness set.]
- P. Breiding, S. Timme, "HomotopyContinuation.jl," (software; the modern
  workhorse).
- Model workflow: "96120: The degree of the linear orbit of a cubic
  surface," arXiv:1909.06620, monodromy + α-certification + trace test
  assembled into a "numerical theorem," with exactly the epistemic honesty
  register this repo's house rules require.

IMPORTS (highest immediate payoff of the whole file):
(a) α-CERTIFY the v103 support-13/-12 states and the v89 (0,3,6) family
    points: converts the cascade's numerical claims (s_Q ≤ 12, deformable
    fiber) into certified theorems without finding closed forms first.
(b) Witness set + numerical irreducible decomposition of the v103 fiber
    variety: degree, number of components, dimension per component, the
    complete answer to "what is this 6-manifold" at numerical-theorem
    strength, then exactify per component. Same for the v_B surface.
(c) Monodromy on the fiber variety computes its GALOIS GROUP, directly
    feeding the holonomy/abelianness program (Conjecture 2) with an
    independent method.

## 5. Tensor rank vs border rank: the border-support question, already a field

THE MATCH. Minimal support of an exact realization = a fermionic RANK;
minimal support achievable in the closure = fermionic BORDER RANK. For
tensors these famously differ, and their gap is the source of ill-posedness.

Sources:
- V. de Silva, L.-H. Lim, "Tensor rank and the ill-posedness of the best
  low-rank approximation problem," SIAM J. Matrix Anal. Appl. 30 (2008)
  1084–1127. [Border rank < rank on positive-measure sets; no best low-rank
  approximant.]
- J.M. Landsberg, "Tensors: Geometry and Applications," AMS GSM 128 (2012).
- L. Chiantini, G. Ottaviani, N. Vannieuwenhoven, "An algorithm for generic
  and low-rank specific identifiability of complex tensors," SIAM J. Matrix
  Anal. Appl. 35 (2014). [Identifiability = fiber-of-decomposition
  triviality.]
- For two fermions, "Slater rank" is classical: J. Schliemann, J.I. Cirac,
  M. Kuś, M. Lewenstein, D. Loss, Phys. Rev. A 64, 022303 (2001), the
  program's min-support is its N-fermion, spectrum-constrained descendant.

IMPORTS:
(a) Formal definitions: support-rank and border-support-rank of a spectrum;
    the RESEARCH_NEXT "border support" claim becomes a crisp question with
    tensor-land intuition predicting the OPPOSITE of the v89 finding, so
    either fermionic spectral fibers are better behaved (a theorem worth
    having) or a gap instance exists at higher rank (a hunt worth running).
(b) Secant-variety language for the support layers Δ^(k) (they are images
    of determinantal secant constructions); known secant machinery
    (dimension counts, Terracini) may compute dim Δ^(k) for free.

## 6. Euclidean distance degree and algebraic matroids: explain the floors

THE MATCH. The contraction attack and every "floor" measured this session
are critical points of a distance function to the fiber variety. The ED
degree counts exactly these critical points; algebraic matroids say which
coordinate subsets determine fiber points.

Sources:
- J. Draisma, E. Horobeț, G. Ottaviani, B. Sturmfels, R.R. Thomas, "The
  Euclidean distance degree of an algebraic variety," Found. Comput. Math.
  16 (2016) 99–149.
- F. Király, L. Theran, R. Tomioka, "The algebraic combinatorial approach
  for low-rank matrix completion," J. Mach. Learn. Res. 16 (2015). [Algebraic
  matroids for completability, which data subsets pin the object.]

IMPORTS:
(a) ED degree of the fiber variety = principled multi-start budgets and an
    a-priori census of the floor values the attacks can land on (turns the
    "thin missed basin" caveat from the (3,11) audit into arithmetic).
(b) The algebraic matroid of the fiber variety answers "which pair
    occupations / coherences determine the fiber point", the exact
    formalization of the 1-RDM → 2-RDM information-gap program (T7/γ²).

## 7. Sums of squares / Positivstellensatz: industrialize the kills

THE MATCH. The ladder's survivor kills and the witness infeasibility proofs
are polynomial-system infeasibility certificates; SOS gives these
systematically with RATIONAL certificates.

Sources:
- J.B. Lasserre, "Global optimization with polynomials and the problem of
  moments," SIAM J. Optim. 11 (2001) 796–817.
- H. Peyrl, P.A. Parrilo, "Computing sum of squares decompositions with
  rational coefficients," Theor. Comput. Sci. 409 (2008) 269–281.

IMPORT: replace per-support hand branch analyses with an automated pipeline:
LP/SDP infeasibility → rational SOS certificate → exact verification in
sympy. The repo's sos_*.py scripts are the seed. This is what completes the
defect ladder to a theorem.

## 8. GIT / scaling algorithms: the image side, for completeness

- P. Bürgisser, C. Franks, A. Garg, R. Oliveira, M. Walter, A. Wigderson,
  "Efficient algorithms for tensor scaling, quantum marginals, and moment
  polytopes," FOCS 2018; "Towards a theory of non-commutative optimization,"
  FOCS 2019. [Membership and norm-minimization on the image; the fiber
  program is the complementary problem these papers do not treat.]
- M. van den Berg et al., "Computing moment polytopes, with a focus on
  tensors, entanglement and matrix multiplication," arXiv:2510.08336.
  [State of the art on computing the polytopes themselves; confirms the
  fiber side is open territory.]

---

## Ordered action list (mapped to RESEARCH_NEXT problems)

1. α-certify the cascade endpoints and the (0,3,6) family (Tool 4a),
   upgrades this month's numerical findings to certified theorems. [feeds
   the s_Q corrections + Problem 2]
2. Witness set + irreducible decomposition of the v103 fiber and the v_B
   surface via HomotopyContinuation.jl; monodromy Galois groups (4b, 4c).
   [Problems 2, 2b; Conjecture 2]
3. Coherence rigidity matroid + generic-rank theorem for the silence
   Jacobian (Tool 2a). [Problem 2b's open remainder]
4. SOS pipeline for support kills; finish the v89 size-7..9 ladder (Tool
   7). [defect theorem]
5. Add the s_H rung (Hermitian pattern IEP bound) between s_I and s_Q;
   import isospectral-flow explorers (Tools 3a, 3b). [defect hierarchy;
   realization algorithms, Immediate Task 6]
6. State the census analogs of the PR injectivity threshold and of
   rank-vs-border-rank; run the gap hunt at (3,11) when Stage 1 lands
   (Tools 1a, 5a). [Problems 5, 6; border-support question]
7. ED degrees of small fiber varieties to calibrate all future floors
   (Tool 6a). [audit methodology]

## One-line thesis for the paper's related-work section

The inverse-fiber program imports its questions from phase retrieval
(ambiguity geometry), its rigidity from framework theory, its algorithms
from structured inverse eigenvalue problems and numerical algebraic
geometry, and its pathology expectations from tensor rank, applied to a
map, the fermionic spectral moment map, on which none of those fields has
yet worked.

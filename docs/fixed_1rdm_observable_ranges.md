# Certified observable ranges on fixed one-body fibers

## Project objective

For a pure fermionic state \(\Psi\), let

\[
\gamma_\Psi{}_{pq}=\langle\Psi|a_q^\dagger a_p|\Psi\rangle
\]

be its one-particle reduced density matrix. Given a complete target 1-RDM
\(\gamma\) and a Hermitian at-most-two-body observable \(W\), define

\[
F_W^-(\gamma)=
\min_{\Psi\mapsto\gamma}\langle\Psi|W|\Psi\rangle,
\qquad
F_W^+(\gamma)=
\max_{\Psi\mapsto\gamma}\langle\Psi|W|\Psi\rangle.
\]

The project goal is to compute certified intervals

\[
\mathcal I_W(\gamma)=[F_W^-(\gamma),F_W^+(\gamma)]
\]

rather than returning one uncontrolled reconstruction of the missing many-body
information. A narrow interval means the 1-RDM nearly determines the
observable. A wide interval quantifies the correlation information hidden from
all one-body data.

The constrained-search idea is standard in density-matrix functional theory,
and the 1-RDM is the moment map for the orbital-unitary action. The intended
new contribution is an exact, machine-verifiable range calculus: explicit
fiber reductions, exact critical values, observable-identifiability tests, and
software that imports physical Hamiltonians.

## First exact fiber

The first target is a nondegenerate pinned point of the Borland--Dennis system
\(\wedge^3\mathbb C^6\):

\[
\gamma=
\operatorname{diag}\left(
\frac{17}{20},\frac34,\frac35,\frac25,\frac14,\frac{3}{20}
\right).
\]

Its occupations are strictly decreasing and satisfy

\[
\lambda_1+\lambda_6=
\lambda_2+\lambda_5=
\lambda_3+\lambda_4=1
\]

and the saturated generalized Pauli constraint

\[
2-(\lambda_1+\lambda_2+\lambda_4)=0.
\]

The exact selection rules reduce every pure state with this complete 1-RDM, in
its natural-orbital frame, to

\[
|\Psi(\alpha,\beta)\rangle=
\sqrt{\frac35}|012\rangle
+\frac{e^{i\alpha}}2|034\rangle
+\sqrt{\frac{3}{20}}e^{i\beta}|135\rangle.
\]

This is a global statement for the stated pure-state fiber, not an imposed
three-determinant ansatz. The generator verifies it by intersecting:

1. the eight determinants containing exactly one orbital from each pair
   \((0,5),(1,4),(2,3)\); and
2. the nine-determinant kernel of \(2-(n_0+n_1+n_3)\).

The intersection is exactly

\[
S=\{012,034,135\}.
\]

No two determinants in \(S\) differ by one excitation, so every one-body
off-diagonal matrix element vanishes identically. The singleton-incidence
matrix has rank three and a unimodular maximal minor, so the diagonal of
\(\gamma\) uniquely fixes

\[
(w_0,w_1,w_2)=\left(\frac35,\frac14,\frac{3}{20}\right).
\]

After removing global phase, the fixed-1RDM fiber is a two-torus of relative
phases.

## Every real carrier Hamiltonian is physical at order two

For \(p<q\), use the pair-annihilator convention

\[
A_{pq}=a_q a_p.
\]

The signed channel enumerator finds one collision-free two-body channel for
each carrier edge:

| Carrier edge | Channel | Fermionic coefficient |
|---|---|---:|
| \(012\leftrightarrow034\) | \((12),(34)\) | \(+1\) |
| \(012\leftrightarrow135\) | \((02),(35)\) | \(-1\) |
| \(034\leftrightarrow135\) | \((04),(15)\) | \(+1\) |

The singleton-incidence minor on orbitals \(0,1,2\) is invertible, so selected
number operators realize an arbitrary diagonal on the carrier. The three
collision-free pair hoppings independently realize its three real
off-diagonal entries. Consequently every rational real-symmetric matrix

\[
H_S=\begin{pmatrix}
h_{00}&h_{01}&h_{02}\\
h_{01}&h_{11}&h_{12}\\
h_{02}&h_{12}&h_{22}
\end{pmatrix}
\]

is the restriction of an explicitly emitted at-most-two-body observable.
This makes the first range engine a complete solver for the real-symmetric
observable algebra on this fiber, not just a solver for one hand-selected
operator.

## Exact constrained-search theorem

Write

\[
E_H(\alpha,\beta)=
E_0+A\cos\alpha+B\cos\beta+C\cos(\beta-\alpha),
\]

where

\[
\begin{aligned}
E_0&=\sum_a w_a h_{aa},\\
A&=2\sqrt{w_0w_1}\,h_{01},\\
B&=2\sqrt{w_0w_2}\,h_{02},\\
C&=2\sqrt{w_1w_2}\,h_{12}.
\end{aligned}
\]

For fixed \(t=\cos\alpha\), optimization over \(\beta\) is elementary:

\[
E_\pm(t)=E_0+At\pm\sqrt{B^2+C^2+2BCt},
\qquad -1\le t\le1.
\]

The upper branch is concave and the lower branch is convex. Therefore the
global maximum and minimum occur among the four branch endpoints at
\(t=\pm1\) and at most one feasible interior stationary value. When
\(ABC\ne0\), the only possible interior cosine is

\[
t_*=\frac{(BC/A)^2-B^2-C^2}{2BC},
\]

with critical value

\[
E_*=E_0+At_*-\frac{BC}{A}.
\]

If \(ABC>0\), this candidate belongs to the lower branch; if \(ABC<0\), it
belongs to the upper branch. It is retained only when \(-1\le t_*\le1\).
Thus every rational real-symmetric carrier matrix is solved by exact comparison
of at most five values. For this target 1-RDM all values lie in
\(\mathbb Q(\sqrt{15})\), and the implementation compares them without
floating-point arithmetic.

The emitted certificate also includes a realizable triple

\[
\bigl(\cos\alpha,\cos\beta,\cos(\beta-\alpha)\bigr)
\]

for each selected endpoint. Its cosine Gram determinant is checked exactly.

## Example 1: unfrustrated coherence triangle

Absorbing the fermionic channel signs into the operator gives the carrier
matrix

\[
H_+=\begin{pmatrix}0&1&1\\1&0&1\\1&1&0\end{pmatrix}.
\]

Its expectation is

\[
E_+(\alpha,\beta)=
\frac{\sqrt{15}}5\cos\alpha
+\frac35\cos\beta
+\frac{\sqrt{15}}{10}\cos(\beta-\alpha).
\]

It also has the phasor form

\[
E_+=
\left|
\sqrt{\frac35}+\frac{e^{i\alpha}}2
+\sqrt{\frac{3}{20}}e^{i\beta}
\right|^2-1.
\]

The three radii close a triangle. Hence

\[
\boxed{
F_{H_+}^-(\gamma)=-1,
\qquad
F_{H_+}^+(\gamma)=\frac{3}{10}(2+\sqrt{15}).
}
\]

The minimum is the interior stationary point with

\[
\cos\alpha=-\frac{7\sqrt{15}}{30},\qquad
\cos\beta=-\frac56,\qquad
\cos(\beta-\alpha)=\frac{2\sqrt{15}}{15}.
\]

The maximum is the aligned phase point \(\alpha=\beta=0\).

The same interval follows independently from the complete-coherence polygon
identity

\[
2\sum_{a<b}\sqrt{w_aw_b}\cos(\theta_b-\theta_a)
=\left|\sum_a\sqrt{w_a}e^{i\theta_a}\right|^2-1.
\]

For any number of carrier vertices, this gives the exact range

\[
\left(\max\left\{0,2r_{\max}-\sum_a r_a\right\}\right)^2-1
\le \langle W_S\rangle \le
\left(\sum_a r_a\right)^2-1,
\qquad r_a=\sqrt{w_a},
\]

whenever every carrier pair has a signed collision-free observable channel.

## Example 2: a \(\pi\)-flux triangle

Flip one carrier edge:

\[
H_\pi=\begin{pmatrix}0&1&1\\1&0&-1\\1&-1&0\end{pmatrix}.
\]

The product of the three real edge signs is negative and cannot be removed by
rephasing the carrier basis. The exact interval is

\[
\boxed{
F_{H_\pi}^-(\gamma)=-\frac{3}{10}(2+\sqrt{15}),
\qquad
F_{H_\pi}^+(\gamma)=1.
}
\]

Now the interior stationary point supplies the maximum:

\[
\cos\alpha=\frac{7\sqrt{15}}{30},\qquad
\cos\beta=\frac56,\qquad
\cos(\beta-\alpha)=\frac{2\sqrt{15}}{15}.
\]

The minimum is the endpoint \((\alpha,\beta)=(\pi,\pi)\). This second case is
important as a regression test: it exercises the frustrated branch rather
than merely changing a scale or adding a constant.

## Why this is useful

The implementation now certifies the full chain

\[
\text{GPC face}
\longrightarrow
\text{exact carrier}
\longrightarrow
\text{fixed-1RDM fiber}
\longrightarrow
\text{physical operator basis}
\longrightarrow
\text{global exact range}.
\]

The result gives exact regression cases for every later solver. A numerical
optimizer that reports an interval outside the certified endpoints is wrong;
a solver that cannot approach both endpoints is incomplete. Because the target
1-RDM is nondegenerate and the support reduction is global, there is no hidden
ansatz error in this benchmark.

For quantum chemistry, the intended output is the same calculation for an
imported interaction operator: rigorous lower and upper correlation energies
compatible with a computed or measured 1-RDM. For quantum information,
replacing the observable by variance or Fisher-information expressions can
bound resources left undetermined by one-body data.

## Next implementation stages

### Stage 1: complex Hermitian carrier observables

A general Hermitian triangle carries a gauge-invariant cycle phase. Extend the
exact solver from real signs \(0,\pi\) to rational complex matrix elements.
The stationarity equations can be polynomialized with tangent-half-angle
variables and solved by exact elimination, with isolating intervals for any
higher-degree algebraic critical values.

### Stage 2: complete \((3,6)\) stratification

Enumerate the interior and boundary strata of the Borland--Dennis polytope,
including occupation degeneracies and vanishing coefficients. For each stratum
record:

- selection-rule carrier;
- fixed-1RDM fiber dimension;
- stabilizer dimension and orbit type;
- singular charts;
- supported one- and two-body observable spaces;
- exact or certified-numerical range method.

Degenerate occupations require care because the natural-orbital frame is not
unique and the simple pinning selection rule needs a stabilizer-aware
formulation.

### Stage 3: physical Hamiltonians

Add FCIDUMP and OpenFermion import. Project the one- and two-electron integrals
onto each certified carrier, classify which coefficients are constant on the
fiber, and compute rigorous energy intervals. Initial benchmark families:

- three-site or three-orbital Hubbard-type models;
- minimal-basis three-electron atoms or molecules;
- random symmetry-respecting two-body Hamiltonians for coverage tests.

The decisive comparison is against full configuration interaction: how much
exact energy uncertainty remains after fixing the complete 1-RDM, and how much
a GPC face reduces it.

### Stage 4: stability away from pinning

Exact pinning is rare in realistic calculations. For
\(D(\lambda)=\varepsilon>0\), combine quasipinning selection-rule estimates
with operator-norm bounds to enlarge the pinned interval by a certified error
term. The target form is

\[
\mathcal I_W(\gamma)
\subseteq
\mathcal I_W^{\mathrm{pinned}}(\gamma_0)
+[-C_W\sqrt\varepsilon,C_W\sqrt\varepsilon]
\]

under explicit nondegeneracy and gap assumptions. No such bound is claimed by
the present artifact.

## Scope and nonclaims

The current exact engine covers rational real-symmetric carrier matrices on
one nondegenerate pinned Borland--Dennis fiber. It does not yet cover:

- complex cycle phases;
- mixed-state constrained search;
- degenerate natural occupations;
- unpinned or quasipinned targets;
- arbitrary \((N,d)\) fibers.

The constrained-search framework, moment-map interpretation, Borland--Dennis
conditions, and pinning selection rule are established background. A paper
must distinguish those ingredients from the exact range algorithms,
stratification, physical applications, and stability certificates developed by
this project.

## Literature position

The benchmark relies on established ingredients:

- R. E. Borland and K. Dennis, *J. Phys. B* **5**, 7--15 (1972),
  DOI `10.1088/0022-3700/5/1/009`;
- C. Schilling, D. Gross, and M. Christandl, *Phys. Rev. Lett.* **110**,
  040404 (2013), DOI `10.1103/PhysRevLett.110.040404`, for pinning selection
  rules;
- C.-C. Wang, J. Liebert, M. Penz, and C. Schilling, “Unified Framework for
  Functional Theories of Quantum Systems,” arXiv:`2606.06676` (2026), for the
  finite-dimensional constrained-search and moment-map framework;
- C. L. Benavides-Riveros, T. Wasak, and A. Recati, arXiv:`2311.12596`, for an
  example of extracting a many-body resource through 1-RDM constrained search.

## Reproduction

The generator named in the artifact's `generated_by` field,
`scripts/fixed_1rdm_observable_ranges.py`, was not delivered to this
repository, nor were the standalone verifier and test module the upstream
write-up lists. What is here instead is an independent verifier that
recomputes the whole calculation from the target 1-RDM rather than replaying
the artifact:

```sh
uv run python scripts/verify_fixed_1rdm_observable_ranges.py
uv run pytest -q tests/test_fixed_1rdm_observable_ranges.py
```

It rebuilds the carrier from the two selection rules, solves the six
occupation equations for the weights, derives the couplings, re-derives both
ranges by the branch rule and the polygon identity separately, and checks that
every recorded phase triple is realizable. Comparison against the artifact is
by symbolic equality, so an equivalent radical form still passes: the branch
rule returns the maximum as `sqrt(15)/5 + sqrt(12*sqrt(15) + 51)/10`, which
denests to the recorded `3(2 + sqrt(15))/10`.

Artifact:

```text
results/data/fixed_1rdm_observable_ranges.json
```

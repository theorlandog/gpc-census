# Pair coordinates and exact fiber symmetries

This note records two exact results about the inverse fibers. The first is a
support-level separation criterion and a complete certificate over the shipped
atlas. The second is an exact stabilizer-invariant coordinate on the known
one-parameter slice of the \(v_B\) reduced fiber.

The scope is important. Pair occupations are expressed in a fixed orbital
frame, so they are not invariants of the full degeneracy-stabilizer quotient.
The cubic 2-RDM moment used for \(v_B\) is invariant under every one-body
unitary. Neither result globally integrates the second direction of the
two-dimensional phase-gauge-quotiented support tangent at \(v_B\).

## 1. Pair-incidence separation

Let \(S=\{T_1,\ldots,T_m\}\) be a support of \(N\)-element determinants and
write \(w_\alpha=|c_{T_\alpha}|^2\). Define the singleton and pair-incidence
matrices

\[
 (A_S)_{i,\alpha}=\mathbf 1_{i\in T_\alpha},\qquad
 (B_S)_{\{i,j\},\alpha}=\mathbf 1_{\{i,j\}\subset T_\alpha}.
\]

The diagonal one-body and two-body occupations are \(A_Sw\) and \(B_Sw\).
For every orbital \(i\),

\[
 \sum_{j\ne i}(B_Sw)_{\{i,j\}}=(N-1)(A_Sw)_i.
\]

Consequently,

\[
 \ker B_S\subseteq\ker A_S.
\]

This proves the following conditional statement.

**Pair-coordinate proposition.** If \(B_S\) has full column rank, the
diagonal pair occupations determine the complete squared-amplitude vector on
\(S\). In particular, every nonzero incidence-kernel weight deformation
changes at least one pair occupation.

The condition is not automatic. A nonzero vector in \(\ker B_S\) is a signed
combinatorial 2-trade, and such supports exist. The project therefore checks
the rank rather than treating it as a theorem about arbitrary supports.

### Exact atlas result

For every one of the 799 certified supports, \(B_S\) has full column rank. The
artifact contains a nonzero integer maximal minor for each support:

- 799 of 799 records certified;
- 403 transport classes represented;
- 63 records have positive incidence-kernel dimension;
- those 63 records comprise 39 transport classes;
- on all 63, the pair map is injective on the complete incidence kernel.

Thus every weight-sector deformation in the shipped atlas is visible to a
diagonal two-body observable in the chosen frame. This proves the
support-restricted form of the proposed cumulant-separation mechanism for the
atlas. It does not prove that every component or every support over the same
vertex has this property.

The certificate is functorial under the project transports. Orbital
relabeling only permutes rows and columns. Empty-orbital padding and
frozen-core lift retain the old pair rows. For Hodge complement, assuming
both particle and hole number are at least two,

\[
 B^{\mathrm c}_{ij}=1-A_i-A_j+B_{ij}.
\]

If \(B^{\mathrm c}x=0\), the singleton identity for the complemented
determinants first gives \(A^{\mathrm c}x=0\), hence
\(\sum_\alpha x_\alpha=0\) and \(A x=0\). The displayed identity then gives
\(Bx=0\). Applying the same argument twice proves rank preservation in both
directions.

## 2. Time reversal on the exact \(v_B\) slice

Use the repository's zero-based determinant support

\[
\begin{split}
S_B=\{&(0,1,2,4),(0,1,2,8),(0,1,3,5),(0,2,3,6),\\
      &(0,3,4,7),(0,3,7,8),(1,3,6,8),(2,3,4,8)\}.
\end{split}
\]

The exact fixed-spectrum family has squared weights

\[
\frac1{23}(1+t,\,8-t,\,4,\,3,\,2-t,\,2+t,\,1,\,2).
\]

Fixed-spectrum, not fixed-diagonal: every member has 1-RDM spectrum exactly
\(\lambda_{v_B}\), but \(\rho\) is not diagonal in the frame this support is
written in, and its off-diagonal moves along the family. At the interior point
\(t=0\) used in section 4 below,

\[
 \rho_{4,8}=-\frac{5}{92}-\frac{\sqrt{119}}{92}\,i\neq0,
 \qquad
 \operatorname{diag}\rho=\tfrac1{23}(20,14,14,14,5,4,4,4,13),
\]

with spectrum \(\{20/23,\,(14/23)^4,\,(4/23)^4\}\). The shipped `states.jsonl`
representative is likewise non-diagonal (\(\rho_{4,8}=3/23\)); the
natural-orbital frame is the separate ledger
`results/data/states_natural_orbital.jsonl`. Everything below is therefore a
fixed-spectrum statement, which is what the constraint
\(P_\lambda\,d\rho\,P_\lambda=0\) in section 4 imposes.

Put

\[
Q(t)=(1+t)(8-t)(2-t)(2+t),\qquad
P(t)=3+7t-2t^2.
\]

The coherence equation fixes the gauge-invariant phase through

\[
\cos\theta=\frac{P(t)}{2\sqrt{Q(t)}}.
\]

For an interior physical value of \(t\), there are two solutions
\(\theta\) and \(-\theta\). Complex conjugation is therefore the involution

\[
\mathcal T:(t,\theta)\longmapsto(t,-\theta).
\]

Modulo the orbital-phase gauge, its fixed points on this slice are exactly the
two real walls

\[
t_\pm=\frac7{17}\pm\frac{18\sqrt{35}}{85}.
\]

The pair occupation on zero-based modes \((1,4)\) is

\[
\langle n_1n_4\rangle=\frac{1+t}{23}.
\]

It recovers \(t\), but it is frame-dependent and does not by itself separate
points in the reduced quotient. A fixed-frame 2-RDM coherence detects the two
time-reversal sheets:

\[
(\Gamma_2)_{(2,4),(3,5)}
=\frac{2\sqrt{1+t}}{23}e^{-i\theta}.
\]

Its imaginary part changes sign under \(\mathcal T\) and vanishes at the real
walls. This proves sheet separation modulo the diagonal orbital-phase gauge.
Whether the two interior sheets at fixed \(t\) are equivalent under the full
\(U(1)\times U(4)\times U(4)\) degeneracy stabilizer remains open.

An exact enumeration gives only the identity among the
\(S_4\times S_4\) block permutations that preserve this determinant support.
This rules out a simple support-preserving orbital relabeling, not a general
continuous block unitary.

## 3. A full one-body-unitary invariant

Let \(\Gamma_2\) use the pair-vector Gram convention in
`gpc_census.gauge`, so \(\operatorname{tr}\Gamma_2=\binom42=6\). Exact
symbolic expansion on the family gives, before inserting the coherence
equation,

\[
\operatorname{tr}(\Gamma_2^2)
=\frac{2}{23^2}
\left(725+14t-4t^2-4\cos\theta\sqrt{Q(t)}\right),
\]

\[
\operatorname{tr}(\Gamma_2^3)
=\frac{6}{23^3}
\left(2932+69t-22t^2-22\cos\theta\sqrt{Q(t)}\right).
\]

Substituting \(\cos\theta\sqrt Q=P/2\) yields

\[
\boxed{\operatorname{tr}(\Gamma_2^2)=\frac{1438}{529}},
\qquad
\boxed{\operatorname{tr}(\Gamma_2^3)
=\frac{6(2899-8t)}{12167}}.
\]

The quadratic moment is blind to motion along the family. The cubic moment has
constant derivative

\[
\frac{d}{dt}\operatorname{tr}(\Gamma_2^3)=-\frac{48}{12167},
\]

and recovers the parameter exactly:

\[
t=\frac{2899-\frac{12167}{6}\operatorname{tr}(\Gamma_2^3)}8.
\]

The spectrum of \(\Gamma_2\), and hence every trace power, is invariant under
the full one-body group \(U(9)\). It follows that two points of this family
with distinct \(t\)-values cannot be related by any one-body unitary. The
known family therefore survives as a genuine one-dimensional locus in the
reduced fiber.

At the two Galois-conjugate real walls, the cubic moments differ by

\[
\operatorname{tr}(\Gamma_2^3)(t_-)
-\operatorname{tr}(\Gamma_2^3)(t_+)
=\frac{1728\sqrt{35}}{1034195}\ne0.
\]

Thus the Galois involution exchanging the two wall states is not a physical
one-body symmetry. Arithmetic conjugacy and degeneracy-stabilizer equivalence
are distinct actions on this fiber.

## 4. Exact first-order coordinates on the full support tangent

The cubic invariant is a global coordinate on the known one-parameter slice.
The support-restricted fixed-spectrum tangent at an interior point has a
second direction. At

\[
t=0,\qquad
\cos\theta=\frac{3\sqrt2}{16},\qquad
\sin\theta=\frac{\sqrt{238}}{16},
\]

the exact spectral projectors of the 1-RDM give a real constraint matrix
\(C\) for

\[
P_\lambda\,d\rho\,P_\lambda=0
\]

at each of the three distinct eigenvalues, together with normalization. In
the 16 real amplitude coordinates on the support,

\[
\operatorname{rank}C=7,\qquad
\dim\ker C=9.
\]

The orbital-phase gauge has exact rank seven and lies in \(\ker C\), so the
support-restricted fixed-spectrum tangent modulo that phase gauge is exactly
two-dimensional. This is the phase-gauge quotient specified by the support
fiber \(F_{S,\lambda}\), not by itself the further quotient under the full
degeneracy stabilizer.

Let

\[
I_3=\operatorname{tr}(\Gamma_2^3),\qquad
I_4=\operatorname{tr}(\Gamma_2^4).
\]

Both differentials annihilate the gauge. Appending \(dI_3\) and \(dI_4\) to
seven independent rows of \(C\) gives rank nine. One exact \(9\times9\)
minor is

\[
\frac{299925504\sqrt{8211}}{1035662780341225}\ne0.
\]

Therefore

\[
(dI_3,dI_4):\ker C/\mathfrak g_{\mathrm{phase}}\longrightarrow\mathbb R^2
\]

is an isomorphism. The cubic and quartic 2-RDM moments form an exact
first-order coordinate pair on the complete two-dimensional
orbital-phase-quotiented support tangent at this point. Because the moments
are invariant under all of \(U(9)\), this also rules out any additional
one-body-orbit direction inside that support tangent.

This is deliberately a tangent statement. It does not prove that every
tangent direction integrates, that the full reduced germ is smooth, or that
\((I_3,I_4)\) is globally injective on the surface. Existing continuation is
numerical evidence for integrability of the second direction. The exact
calculation supplies the next proof target without promoting that evidence.

## 5. Verification

The reproducible files are:

- `scripts/fiber_symmetries.py`, which regenerates all exact minors and the
  symbolic \(v_B\) calculation;
- `results/data/fiber_symmetries.json`, containing 799 integer minor
  certificates and the exact \(v_B\) formulas;
- `scripts/verify_fiber_symmetries_standalone.py`, a standard-library-only
  Bareiss verifier for every minor;
- `tests/test_fiber_symmetries.py`, including regeneration, tamper rejection,
  symbolic moment, exact tangent-rank, and Galois-wall tests.

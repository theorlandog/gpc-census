# Exact scheduled two-body motion on the known v_B carrier

## Result and scope

The known algebraic \(v_B\) circle has an exact, sparse, scheduled
Hamiltonian that preserves its complete eight-determinant carrier. A specified
orbital lift gives a second schedule whose reference trajectory has one
constant complete 1-RDM. Both schedules are at most two-body.

Combined with the exact carrier reconstruction certificate in
`results/data/rdm_observability_census.json`, this gives an exact closed
nonautonomous TD-2RDM trajectory. This is more than kinematic reconstruction:
an ordinary two-body Schrödinger evolution is proved to remain in the chart
where \(\Gamma_3\) is a rational function of \(\Gamma_2\).

Orbital labels are zero-based. For a sorted pair \(p=(i,j)\), use

\[
A_p=a_j a_i,
\qquad
H_K=\sum_{p,q}K_{pq}A_p^\dagger A_q.
\]

The 2-RDM convention is

\[
(\Gamma_2)_{pq}=\langle A_p^\dagger A_q\rangle,
\qquad
\operatorname{Tr}\Gamma_2=\binom42=6.
\]

## Carrier and wall-regular clock

In repository order, the carrier is

\[
\begin{aligned}
D_0&=(0,1,2,4), &D_1&=(0,1,2,8),\\
D_2&=(0,1,3,5), &D_3&=(0,2,3,6),\\
D_4&=(0,3,4,7), &D_5&=(0,3,7,8),\\
D_6&=(1,3,6,8), &D_7&=(2,3,4,8).
\end{aligned}
\]

Introduce a signed coordinate \(s\) on the conic

\[
s^2=119+70t-85t^2.
\]

An explicit periodic clock is

\[
t(\phi)=\frac7{17}+\frac{18\sqrt{35}}{85}\cos\phi,
\qquad
s(\phi)=\frac{18\sqrt{119}}{17}\sin\phi,
\qquad 0\leq\phi\leq2\pi.
\]

Equivalently,

\[
\dot t=-\frac{s}{\sqrt{85}},
\qquad
\dot s=\sqrt{85}\left(t-\frac7{17}\right),
\]

where a dot means \(d/d\phi\). This vector field is tangent to the conic and
is regular at both real walls, where \(s=0\).

In the carrier gauge, put the phase on \(D_0\):

\[
e=e^{i\theta}
=\frac{3+7t-2t^2+i s}
 {2\sqrt{(1+t)(8-t)(2-t)(2+t)}}.
\]

Exact reduction modulo the conic gives \(|e|=1\) and

\[
\dot e=i\omega e,
\qquad
\omega=
\frac{85t^3-105t^2+12t+364}
 {2\sqrt{85}(1+t)(8-t)(2-t)(2+t)}.
\]

All four radicands are strictly positive on the complete physical interval.

## Seven-entry pair kernel

Define the real transfer rates

\[
a=-\frac{s}{2\sqrt{85}\sqrt{(1+t)(8-t)}},
\qquad
b=\frac{s}{2\sqrt{85}\sqrt{(2-t)(2+t)}}.
\]

Set

\[
\mathcal A=-ia\bar e,
\qquad
\mathcal B=ib,
\]

and

\[
\alpha_0=\frac{\mathcal A+\mathcal B}{2},
\qquad
\alpha_1=\frac{\mathcal A-\mathcal B}{2},
\qquad
\alpha_3=-\alpha_1.
\]

The carrier driver is

\[
H_{\rm car}=-\omega n_1n_4+T+T^\dagger,
\]

with density-assisted hopping

\[
T=a_8^\dagger a_4
  \left(\alpha_0n_0+\alpha_1n_1+\alpha_3n_3\right)
=a_8^\dagger a_4
  \left(\alpha_0n_0+\alpha_1(n_1-n_3)\right).
\]

Because \(j<4<8\) for \(j\in\{0,1,3\}\),

\[
a_8^\dagger a_4 n_j=A_{(j,8)}^\dagger A_{(j,4)}.
\]

The symbolic pair-kernel support consists of

\[
K_{14,14}=-\omega,
\qquad
K_{j8,j4}=\alpha_j,
\qquad
K_{j4,j8}=\bar\alpha_j,
\quad j\in\{0,1,3\}.
\]

Thus the Hermitian pair kernel has a seven-entry symbolic support. At either
wall the six hopping coefficients vanish, while the diagonal entry remains
finite.

## Why the carrier is invariant

The hopping \(a_8^\dagger a_4\) has only three relevant actions on or next to
the carrier:

\[
D_0\longrightarrow D_1,
\qquad
D_4\longrightarrow D_5,
\qquad
(1,3,4,6)\longrightarrow D_6.
\]

Their exact signed coefficients are

\[
\alpha_0+\alpha_1=\mathcal A,
\qquad
-(\alpha_0+\alpha_3)=-\mathcal B,
\qquad
-(\alpha_1+\alpha_3)=0.
\]

The last identity cancels the sole leakage channel. Hermiticity cancels the
reverse leakage. The density term is diagonal, so the full carrier is an
invariant subspace of \(H_{\rm car}(\phi)\), not merely tangent at one point.

On the carrier, the target state is

\[
\psi=\frac1{\sqrt{23}}
\left(
\sqrt{1+t}\,e,
\sqrt{8-t},
2,
\sqrt3,
\sqrt{2-t},
\sqrt{2+t},
1,
\sqrt2
\right).
\]

The only five nonzero restricted matrix entries are

\[
H_{00}=-\omega,
\quad H_{01}=iae,
\quad H_{10}=-ia\bar e,
\quad H_{45}=ib,
\quad H_{54}=-ib.
\]

Direct substitution gives the vector identity

\[
i\frac{d\psi}{d\phi}=H_{\rm car}(\phi)\psi(\phi).
\]

There is no unreported projective scalar. At either wall, \(a=b=0\), while
\(\omega\) remains finite, so the schedule closes smoothly.

## Constant complete 1-RDM lift

Let

\[
U=e^{i\beta n_8},
\qquad
u=e^{i\beta}
=\frac{7t-5+i s}{6\sqrt{(2-t)(2+t)}}.
\]

Conic reduction gives \(|u|=1\) and

\[
\dot u=i\nu u,
\qquad
\nu=\frac{28-5t}{\sqrt{85}(4-t^2)}.
\]

The lifted Hamiltonian is

\[
H_{\rm common}=UH_{\rm car}U^\dagger-\nu n_8.
\]

Equivalently, multiply each \(\alpha_j\) in \(T\) by \(u\), retain the
Hermitian conjugate, and add \(-\nu n_8\). The one-body term is allowed in an
at-most-two-body Hamiltonian. On the canonical four-particle sector it also
has the pair-only representation

\[
n_8=\frac13\sum_{j=0}^7 n_jn_8.
\]

The lifted state obeys

\[
i\frac{d(U\psi)}{d\phi}=H_{\rm common}(U\psi),
\]

and its complete 1-RDM is exactly

\[
\rho_B=\frac1{23}
\left[
\operatorname{diag}(20,14,14,14,5,4,4,4,13)
+3\left(|8\rangle\langle4|+|4\rangle\langle8|\right)
\right]
\]

for every \(\phi\). This is a fixed matrix, not only a fixed diagonal or fixed
spectrum.

## Exact TD-2RDM closure

The carrier has a unimodular \(8\) by \(8\) pair-incidence minor and a
connected seven-edge graph of collision-free 2-RDM coherences. Therefore, on
the pure full-support carrier chart,

\[
\Gamma_3=F_{S,3}(\Gamma_2)
\]

with \(F_{S,3}\) rational. Every weight stays strictly positive around the
circle, and both Hamiltonians preserve the carrier exactly. The BBGKY equation
therefore closes as

\[
\frac{d\Gamma_2}{d\phi}
=\mathcal L_{H(\phi)}
\left(\Gamma_2,F_{S,3}(\Gamma_2)\right).
\]

This is exact and nonautonomous. It is an analytically soluble ground-truth
instance for time-dependent 2-RDM reconstruction and propagation. It does not
solve generic TD-2RDM closure.

The current package proves and parameterizes this closure, but does not yet
emit the rational evaluator \(F_{S,3}\) or integrate the reduced equation
without consulting carrier amplitudes. That executable propagator is the next
computational deliverable.

The coefficient-weight tangent is
\((1,-1,0,0,-1,1,0,0)\). Its pair-diagonal image is nonzero, while the lifted
complete 1-RDM is constant. The protocol therefore moves two-body correlation
data that are invisible to every static one-body expectation. Because the
disconnected two-body contribution is fixed by the complete 1-RDM, the
changing part is the two-body cumulant.

There is also an exact basis-free witness:

\[
\operatorname{Tr}\left[(\Gamma_2)^3\right]
=\frac{6(2899-8t)}{12167},
\qquad
\frac{d}{dt}\operatorname{Tr}\left[(\Gamma_2)^3\right]
=-\frac{48}{12167}.
\]

Here \(\Gamma_2[P,Q]=\langle A_P^\dagger A_Q\rangle\) in the unordered-pair
basis, so \(\operatorname{Tr}\Gamma_2=\binom42=6\). A unit-trace convention
would rescale the cubic invariant by \(6^{-3}\).

This spectral invariant changes along the loop. Therefore the trajectory is
not merely a one-body orbital rotation, including a rotation inside the
stabilizer of the fixed 1-RDM.

## Conditional shortcut construction

If a separate at-most-two-body scheduled parent obeys

\[
H_0(\phi)\psi(\phi)=0,
\]

then for any laboratory schedule \(\dot\phi=v(\tau)\),

\[
H_{\rm total}(\tau)
=H_0(\phi(\tau))+v(\tau)H_{\rm car}(\phi(\tau))
\]

tracks the carrier state exactly. The common-1RDM version uses
\(H_{\rm common}\). This conditional identity does not assert that a suitable
parent, gap, or experimental implementation has been supplied here.

## Limits

- The result is a precomputed time-dependent schedule, not an autonomous
  Hamiltonian.
- Neighboring initial states are not claimed to keep the target 1-RDM.
- No spatial locality, noise tolerance, or driver gap is proved.
- The pair-only formula for \(n_8\) is restricted to the \(N=4\) sector.
- The rational TD-2RDM closure is not asserted for mixed states or outside the
  certified carrier chart.
- The closure is certified but an amplitude-free reduced propagator is not yet
  implemented.
- No global statement about the full fixed-1RDM fiber or its dimension is
  made.

The machine-readable certificate is
`results/data/vb_carrier_hamiltonian.json`. It is generated by
`scripts/vb_carrier_hamiltonian.py` and independently checked by
`scripts/verify_vb_carrier_hamiltonian_standalone.py`.

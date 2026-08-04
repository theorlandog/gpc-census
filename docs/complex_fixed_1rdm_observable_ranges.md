# Exact complex-Hermitian observable ranges on the fixed Borland-Dennis fiber

## Result

Fix the nondegenerate Borland-Dennis one-particle density matrix

\[
\gamma=\operatorname{diag}
\left(
\frac{17}{20},\frac34,\frac35,\frac25,\frac14,\frac3{20}
\right).
\]

The pinned selection rule, strict occupation ordering, and incidence certificate
reduce the complete pure-state fiber to

\[
|\Psi(\alpha,\beta)\rangle
=
\sqrt{\frac35}|012\rangle
+\frac{e^{i\alpha}}2|034\rangle
+\sqrt{\frac3{20}}e^{i\beta}|135\rangle.
\]

This note extends the real-symmetric observable-range solver to rational
Gaussian Hermitian carrier matrices

\[
H=
\begin{pmatrix}
h_{00}&h_{01}&h_{02}\\
\overline{h}_{01}&h_{11}&h_{12}\\
\overline{h}_{02}&\overline{h}_{12}&h_{22}
\end{pmatrix},
\qquad
h_{jk}\in\mathbb Q(i),
\quad h_{jj}\in\mathbb Q.
\]

For a checked generic nonzero-flux input, the exact constrained range is

\[
\boxed{
\mathcal I_H(\gamma)
=
\left[E_0+z_{\min},E_0+z_{\max}\right],
}
\]

where \(z_{\min}\) and \(z_{\max}\) are the smallest and largest real roots of
an explicitly generated primitive sextic polynomial in \(\mathbb Z[z]\).
The roots are selected by exact Sturm isolation. No floating-point optimizer,
mesh, or local critical-point search determines either endpoint.

The implementation also emits a rational linear subresultant lift. It proves
that every real sextic root comes from a real point of the phase torus, rather
than from a complex solution introduced by elimination.

## 1. The cycle phase is the new invariant

Write the three weighted edge phasors as

\[
a=\frac{\sqrt{15}}5h_{01},
\qquad
b=\frac35h_{02},
\qquad
c=\frac{\sqrt{15}}{10}h_{12}.
\]

Then

\[
E(\alpha,\beta)
=E_0
+\operatorname{Re}(ae^{i\alpha})
+\operatorname{Re}(be^{i\beta})
+\operatorname{Re}(ce^{i(\beta-\alpha)}),
\]

with

\[
E_0=\frac35h_{00}+\frac14h_{11}+\frac3{20}h_{22}.
\]

A carrier-basis rephasing sends

\[
h_{jk}\longmapsto e^{-i\theta_j}h_{jk}e^{i\theta_k}.
\]

The individual edge phases are gauge-dependent, but the cycle product

\[
\boxed{
\chi=h_{01}h_{12}\overline{h}_{02}
}
\]

is invariant. Its argument is the triangle flux. Real-symmetric matrices have
\(\arg\chi=0\) or \(\pi\); a genuinely complex triangle has
\(\operatorname{Im}\chi\ne0\).

This distinction is not cosmetic. Two matrices with identical diagonal terms
and edge magnitudes can have different constrained ranges because the three
phase preferences cannot all be gauged away around a cycle.

## 2. Eliminate one phase exactly

At fixed \(\alpha\), combine the two terms containing \(\beta\):

\[
\operatorname{Re}(be^{i\beta})
+\operatorname{Re}(ce^{i(\beta-\alpha)})
=
\operatorname{Re}\!\left[(b+ce^{-i\alpha})e^{i\beta}\right].
\]

Therefore the exact upper and lower branches are

\[
E_\sigma(\alpha)
=E_0+\operatorname{Re}(ae^{i\alpha})
+\sigma\left|b+ce^{-i\alpha}\right|,
\qquad \sigma\in\{-1,+1\}.
\]

Introduce

\[
u_\alpha=(\cos\alpha,\sin\alpha),
\]

and real vectors \(p,d\in\mathbb R^2\) such that

\[
\operatorname{Re}(ae^{i\alpha})=p\cdot u_\alpha,
\qquad
|b+ce^{-i\alpha}|^2=R+d\cdot u_\alpha.
\]

Explicitly,

\[
p=(\operatorname{Re}a,-\operatorname{Im}a),
\]

\[
d=
\left(
2\operatorname{Re}(b\overline c),
-2\operatorname{Im}(b\overline c)
\right),
\qquad
R=|b|^2+|c|^2.
\]

Define the four scalar invariants

\[
A=|p|^2,
\qquad
D=|d|^2,
\qquad
C=p\cdot d,
\qquad
R=|b|^2+|c|^2.
\]

For the fixed weights and rational Gaussian matrix entries, all four lie in
\(\mathbb Q\):

\[
\boxed{
\begin{aligned}
A&=\frac35|h_{01}|^2,\\
D&=\frac{27}{125}|h_{02}|^2|h_{12}|^2,\\
C&=\frac9{25}\operatorname{Re}\chi,\\
R&=\frac9{25}|h_{02}|^2+\frac3{20}|h_{12}|^2.
\end{aligned}
}
\]

The oriented area is

\[
\delta=\det(p,d)=\frac9{25}\operatorname{Im}\chi,
\]

and hence

\[
\boxed{AD-C^2=\delta^2.}
\]

Thus nonzero gauge-invariant flux is exactly the condition that \(p\) and
\(d\) form a basis of the phase plane.

## 3. The quartic branch curve

Let

\[
z=E-E_0,
\qquad
y=\sigma\sqrt{R+d\cdot u_\alpha}.
\]

Then

\[
p\cdot u_\alpha=z-y,
\qquad
d\cdot u_\alpha=y^2-R.
\]

For any unit vector \(u\in\mathbb R^2\), its two projections obey the Gram
identity

\[
D(p\cdot u)^2
-2C(p\cdot u)(d\cdot u)
+A(d\cdot u)^2
=AD-C^2.
\]

Substitution gives the quartic plane curve

\[
\boxed{
F(z,y)=
D(z-y)^2
-2C(z-y)(y^2-R)
+A(y^2-R)^2
-(AD-C^2)=0.
}
\]

When \(\delta\ne0\), this equation is equivalent to the original branch
search. It is not merely a necessary elimination relation. Indeed, from any
real \((z,y)\) on the curve, set

\[
U=z-y,
\qquad V=y^2-R.
\]

The unique phase vector is

\[
\boxed{
\begin{aligned}
\cos\alpha&=\frac{d_2U-p_2V}{\delta},\\
\sin\alpha&=\frac{-d_1U+p_1V}{\delta}.
\end{aligned}
}
\]

The Gram equation gives \(\cos^2\alpha+\sin^2\alpha=1\). Once \(\alpha\)
is recovered, \(\beta\) follows from

\[
\boxed{
e^{i\beta}
=\frac{\overline{b+ce^{-i\alpha}}}{y}.
}
\]

This formula includes both signs of \(y\). It gives the maximizing \(\beta\)
when \(y>0\) and the minimizing \(\beta\) when \(y<0\).

## 4. Critical values are a rational sextic

The physical branch curve is compact. At every global maximum or minimum of
\(z\), including singular points, one has

\[
F(z,y)=0,
\qquad
\frac{\partial F}{\partial y}(z,y)=0.
\]

Eliminating \(y\) gives

\[
\operatorname{Res}_y(F,F_y)
=K(A,C,D,R)P_H(z),
\]

where \(K\ne0\) whenever \(A>0\) and \(AD-C^2>0\), and \(P_H\) is a
primitive degree-six polynomial in \(\mathbb Z[z]\).

The fact that \(P_H\) is rational is structural. The original expectation
contains \(\sqrt{15}\), but the critical values depend only on the four
rotational invariants \(A,C,D,R\), all rational for rational Gaussian input.
The complex phase problem therefore moves from \(\mathbb Q(\sqrt{15})\) to a
sextic algebraic extension over \(\mathbb Q\), rather than to uncontrolled
transcendental numerics.

A resultant root can in principle come from complex eliminated variables. The
artifact closes that gap. The subresultant sequence supplies integer
polynomials \(L_A(z)\) and \(L_B(z)\) such that

\[
\boxed{y=-\frac{L_B(z)}{L_A(z)}}
\]

at every root of \(P_H\). The verifier checks exactly that

\[
\gcd(P_H,L_A)=1,
\]

and that, after substitution,

\[
P_H\mid\operatorname{num}F\!\left(z,-\frac{L_B}{L_A}\right),
\]

\[
P_H\mid\operatorname{num}F_y\!\left(z,-\frac{L_B}{L_A}\right).
\]

Therefore every real root of \(P_H\) has a real critical lift. Conversely,
every global endpoint is a real root of \(P_H\). Under the checked squarefree
and coprime conditions,

\[
\boxed{
\min\mathcal I_H(\gamma)=E_0+\min\operatorname{Roots}_{\mathbb R}(P_H),
}
\]

\[
\boxed{
\max\mathcal I_H(\gamma)=E_0+\max\operatorname{Roots}_{\mathbb R}(P_H).
}
\]

The generator isolates all real roots exactly. The standalone verifier rebuilds
the Sylvester determinant with `fractions.Fraction`, verifies the lift by
polynomial division, and independently checks every isolating interval with a
Sturm sequence.

## 5. A continuous unit-edge flux family

Consider

\[
H(\Phi)=
\begin{pmatrix}
0&1&1\\
1&0&e^{i\Phi}\\
1&e^{-i\Phi}&0
\end{pmatrix}.
\]

All three carrier-edge magnitudes are one and the cycle product is
\(\chi=e^{i\Phi}\). Put

\[
u=\cos\Phi.
\]

The complete critical polynomial is

\[
\boxed{
\begin{aligned}
P_u(z)={}&
1{,}000{,}000z^6
+2{,}000{,}000u z^5\\
&+(1{,}110{,}000u^2-2{,}330{,}000)z^4\\
&+(180{,}000u^3-6{,}060{,}000u)z^3\\
&-(5{,}164{,}200u^2+143{,}700)z^2\\
&-(1{,}798{,}200u^3+57{,}600u)z\\
&-218{,}700u^4+6{,}831u^2+3{,}969.
\end{aligned}
}
\]

The spectrum of constrained expectation values is invariant under complex
conjugation, \(\Phi\mapsto-\Phi\), so it depends on the flux only through
\(u=\cos\Phi\). Algebraically,

\[
\boxed{P_{-u}(-z)=P_u(z).}
\]

Consequently,

\[
z_{\max}(u)=-z_{\min}(-u).
\]

At quarter flux, \(u=0\), the interval is symmetric. Away from quarter flux,
the range is generally asymmetric.

### Real-symmetric limits

At zero flux, \(u=1\), the sextic factors over \(\mathbb Q(\sqrt{15})\), and
its extreme real roots are

\[
-1,
\qquad
\frac3{10}(2+\sqrt{15}).
\]

Thus

\[
\mathcal I_{H(0)}(\gamma)
=
\left[-1,\frac3{10}(2+\sqrt{15})\right].
\]

At \(\pi\)-flux, \(u=-1\), its extreme real roots are

\[
-\frac3{10}(2+\sqrt{15}),
\qquad 1,
\]

so

\[
\mathcal I_{H(\pi)}(\gamma)
=
\left[-\frac3{10}(2+\sqrt{15}),1\right].
\]

The sextic family therefore contains both earlier real-symmetric regression
cases as exact boundary points, rather than defining a disconnected new
solver.

## 6. Exact regression 1: quarter flux

For

\[
H_{\pi/2}=
\begin{pmatrix}
0&1&1\\
1&0&i\\
1&-i&0
\end{pmatrix},
\qquad \chi=i,
\]

one has

\[
A=\frac35,
\qquad
D=\frac{27}{125},
\qquad
C=0,
\qquad
R=\frac{51}{100},
\qquad
\delta=\frac9{25}.
\]

The exact critical polynomial is

\[
P_{\pi/2}(z)
=
1{,}000{,}000z^6
-2{,}330{,}000z^4
-143{,}700z^2
+3{,}969.
\]

It has four real roots. The exact range endpoints are the least and greatest,
with rational isolating intervals

\[
-1.54578275998357
<z_{\min}<
-1.54578275998354,
\]

\[
1.54578275998354
<z_{\max}<
1.54578275998357.
\]

Hence

\[
\boxed{
\mathcal I_{H_{\pi/2}}(\gamma)
=
[-\rho,\rho],
}
\]

where \(\rho\) is the unique positive root of \(P_{\pi/2}\) in the stored
positive endpoint interval. The symmetry follows both from the even polynomial
and from \(P_{-u}(-z)=P_u(z)\) at \(u=0\).

The linear lift is

\[
L_A(z)=1900z^2-99,
\qquad
L_B(z)=-900z^3-201z.
\]

Thus every critical root has

\[
y=\frac{900z^3+201z}{1900z^2-99}.
\]

The verifier checks \(\gcd(P_{\pi/2},L_A)=1\) and both exact divisibility
identities.

## 7. Exact regression 2: a rational non-special flux

Choose the Pythagorean phase

\[
e^{i\Phi}=\frac35+\frac45i.
\]

Then

\[
H_{3,4,5}=
\begin{pmatrix}
0&1&1\\
1&0&\frac35+\frac45i\\
1&\frac35-\frac45i&0
\end{pmatrix},
\qquad
\chi=\frac35+\frac45i.
\]

This is neither zero flux, \(\pi\)-flux, nor quarter flux. The invariants are

\[
A=\frac35,
\qquad
D=\frac{27}{125},
\qquad
C=\frac{27}{125},
\qquad
R=\frac{51}{100},
\qquad
\delta=\frac{36}{125}.
\]

The primitive critical polynomial is

\[
\boxed{
\begin{aligned}
P_{3,4,5}(z)={}&
6{,}250{,}000z^6
+7{,}500{,}000z^5
-12{,}065{,}000z^4\\
&-22{,}482{,}000z^3
-12{,}517{,}575z^2\\
&-2{,}643{,}570z
-136{,}971.
\end{aligned}
}
\]

It has four real roots. The exact endpoint isolating intervals are

\[
-1.34617437837967
<z_{\min}<
-1.34617437837964,
\]

\[
1.68507354505325
<z_{\max}<
1.68507354505328.
\]

Therefore

\[
\boxed{
\mathcal I_{H_{3,4,5}}(\gamma)
=
[z_{\min},z_{\max}],
}
\]

with the endpoints defined exactly as the indicated roots of
\(P_{3,4,5}\). The unequal absolute endpoint values certify that this test is
not a disguised sign flip or rescaling of the quarter-flux example.

Its linear lift is

\[
L_A(z)=75{,}000z^3+251{,}000z^2+96{,}750z-3{,}870,
\]

\[
L_B(z)=-112{,}500z^3-157{,}500z^2-85{,}065z-13{,}311.
\]

Again, the verifier proves exact coprimality and divisibility.

## 8. Explicit at-most-two-body realization

Let

\[
D_0=012,
\qquad D_1=034,
\qquad D_2=135,
\]

and use

\[
A_{ij}=a_j a_i,
\qquad i<j.
\]

The collision-free pair channels have fermionic signs

\[
A_{34}^\dagger A_{12}:D_0\mapsto+D_1,
\]

\[
A_{35}^\dagger A_{02}:D_0\mapsto-D_2,
\]

\[
A_{15}^\dagger A_{04}:D_1\mapsto+D_2.
\]

The pair occupations \(n_0n_1\), \(n_0n_3\), and \(n_1n_3\) restrict to the
three diagonal carrier matrix units. Therefore every complex Hermitian carrier
matrix above is realized by the at-most-two-body operator

\[
\boxed{
\begin{aligned}
\widehat H={}&
 h_{00}n_0n_1+h_{11}n_0n_3+h_{22}n_1n_3\\
&+\overline h_{01}A_{34}^\dagger A_{12}
 +h_{01}A_{12}^\dagger A_{34}\\
&-\overline h_{02}A_{35}^\dagger A_{02}
 -h_{02}A_{02}^\dagger A_{35}\\
&+\overline h_{12}A_{15}^\dagger A_{04}
 +h_{12}A_{04}^\dagger A_{15}.
\end{aligned}
}
\]

Each displayed hopping term annihilates the third carrier determinant and maps
only its intended determinant pair. Hence \(\widehat H\) preserves the full
three-dimensional carrier and restricts exactly to \(H\).

This matters for interpretation: the sextic range is not the range of an
abstract matrix pasted onto a variational ansatz. It is the constrained range
of an explicit ordinary at-most-two-body fermionic observable on the complete
fixed-1RDM fiber.

## 9. What the standalone verifier proves

`scripts/verify_complex_fixed_1rdm_observable_ranges_standalone.py` uses only
the Python standard library. For each regression it independently:

1. parses and verifies the rational Gaussian Hermitian matrix;
2. recomputes \(\chi,A,C,D,R,\delta\) with `fractions.Fraction`;
3. rebuilds the quartic \(F(z,y)\);
4. constructs the \(7\times7\) Sylvester matrix for \(F\) and \(F_y\);
5. evaluates its determinant exactly and recovers the stored primitive sextic;
6. verifies squarefreeness;
7. verifies the rational linear lift and its coprimality;
8. constructs a Sturm sequence;
9. proves that every stored rational interval contains exactly one root;
10. proves that the listed intervals exhaust all real roots;
11. selects the least and greatest real roots;
12. rechecks the three fermionic hopping signs and the diagonal incidence minor.

The decimal interval midpoints are display data. They do not participate in
root selection or proof.

## 10. Scope and limitations

The exact generic engine currently accepts rational Gaussian carrier entries.
This is broad enough for exact flux regressions, Pythagorean phases, and many
symbolic model Hamiltonians, but it is not yet an arbitrary algebraic-number
front end.

The nonzero-flux theorem checks three generic conditions:

- \(\delta\ne0\);
- a squarefree sextic;
- a linear-lift denominator coprime to that sextic.

If one fails, the implementation reports the degeneracy rather than silently
using a numerical optimizer. Zero and \(\pi\) flux are already solved by the
real-symmetric branch engine and appear as exact limits of the unit-flux
sextic. Nongeneric repeated-root complex cases require a higher-subresultant
fallback in a later milestone.

The benchmark remains one fixed three-determinant Borland-Dennis fiber. The
method itself is more general: whenever eliminating all but one phase produces
an algebraic branch curve, its critical-value discriminant and a real-lift
certificate can replace floating constrained optimization. Larger carrier
graphs will generally produce higher-dimensional elimination problems.

## 11. Path to electronic-structure Hamiltonians

The next practical step is not merely to feed arbitrary \(3\times3\) matrices
to the solver. It is to project a recognizable second-quantized Hamiltonian
onto the certified carrier and retain the exact map from one- and two-electron
integrals to \(h_{jk}\).

For this carrier, every off-diagonal edge is a double excitation. Under the
sign convention above, the three restricted matrix elements are supplied by
three signed antisymmetrized two-electron integrals. The projection pipeline
should:

1. parse an FCIDUMP-like integral file exactly where possible;
2. build the three carrier diagonal energies by Slater-Condon rules;
3. build the three signed double-excitation couplings;
4. verify that no omitted term leaks the carrier if a carrier-preserving model
   is claimed;
5. compute the cycle product \(\chi\);
6. dispatch real \(\chi\) to the real branch solver and nonreal \(\chi\) to the
   sextic solver;
7. report a certified fixed-1RDM energy interval and the exact projection
   assumptions.

Ordinary nonrelativistic FCIDUMP files in a real orbital basis usually produce
real carrier couplings, so the earlier zero/\(pi\)-flux engine is already the
relevant path for many molecular examples. Genuinely complex flux is expected
for models with magnetic fields, twisted boundary conditions, complex Bloch
orbitals, or spin-orbit structure. A physically honest implementation should
make that distinction explicit rather than presenting complex flux as generic
for standard real molecular Hamiltonians.

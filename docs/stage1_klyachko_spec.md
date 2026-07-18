# Stage 1: Klyachko constraint generator — spec from Altunbulak thesis
## Target: (3,11)-(3,13); MUST first reproduce ranks 6-10 ground truth

## Core theorem 3.2.1 + coefficient machinery [thesis 1955-2065]

(3.12)

v

Theorem 3.2.1 In the above notations all constraints on the occupation numbers
λ of the system Hrν in a state ρν of spectrum µ are given by the inequalities
X
i

ai λv(i) ≤

X

aνk µw(k)

k

for all test spectra a and permutations v, w such that cvw (a) 6= 0.

(3.13)

CHAPTER 3. ν-REPRESENTABILITY PROBLEM

36

Proof : Follows from Proposition 3.1.1 and Theorem 2.2.1. Remember that the
left action of a permutation on “places” is inverse to its right action on indices.
As a result, the permutations v and w, acting on a and f∗ (a) = aν in Theorem
2.2.1, move to the indices of λ and µ in the inequality (3.13). 
Remark 3.2.1 The coefficient cvw (a) depends only on the order in which quantities aT appear in the spectrum aν . The order changes when the test spectrum
a crosses a hyperplane
HT |T 0 :

X
i∈T

ai =

X

aj .

j∈T 0

The hyperplanes cut the set of all test spectra into a finite number of polyhedral
cones called cubicles. For each cubicle one has to check the inequality (3.13)
only for its extremal edges. As a result, the ν-representability amounts to a finite
system of linear inequalities.

3.2.1

Topological nature of the coefficients cvw (a)

Note that the inequalities (3.13) are subject to the topological condition cvw (a) 6=
0. So one has to calculate the coefficients cvw (a) in order to give Theorem 3.2.1
full strength. We borrow from [1] the following calculation of these coefficients.

Canonical generators
To proceed we first need an alternative description of the cohomology of flag variety Fa (Hr ) [4]. Recall that the latter is understood here as the set of Hermitian
operators in Hr of given spectrum a. To avoid technicalities, we assume the spectrum to be simple, i.e., a1 > a2 > · · · > ar . Let Ei be the eigenbundle on Fa (Hr )
whose fiber at X ∈ Fa (Hr ) is the eigenspace of operator X with eigenvalue ai .
Their Chern classes xi = c1 (Ei ) generate the cohomology ring H ∗ (Fa (Hr )) and we
refer to them as the canonical generators. The elementary symmetric functions
σi (x) of the canonical generators are the characteristic classes of the trivial bundle

37

CHAPTER 3. ν-REPRESENTABILITY PROBLEM

Hr and thus vanish. This identifies the cohomology with the ring of coinvariants
H ∗ (Fa (Hr )) = Z[x1 , x2 , . . . , xr ]/(σ1 , σ2 , . . . , σr ).

(3.14)

This approach to the cohomology is more functorial and by that reason leads to
an easy calculation of the morphism (3.11)
ϕ∗a : H ∗ (Faν (Hν )) → H ∗ (Fa (H)).
Recall that the spectrum aν consists of the quantities aT =

P

i∈T ai arranged in

decreasing order, where T runs over all semistandard tableaux of shape ν. We
P
define xT = i∈T xi in a similar way.
Proposition 3.2.1 Let xi and xνk be the canonical generators of H ∗ (Fa (H)) and
H ∗ (Faν (Hν )) respectively. Then
ϕ∗a (xνk ) = xT ,

when

aνk = aT .

(3.15)

In other words, ϕ∗a (xνk ) is obtained from aνk by the substitution ai 7→ xi .
Proof : The eigenbundle Ei is equivariant with respect to the adjoint action
X 7→ uXu∗ of the unitary group U(H). Therefore it is uniquely determined by
the linear representation of the centralizer D = Z(X) in a fixed fiber Ei (X) or
by its character εi : D → S1 = {z ∈ C∗ | |z| = 1}. In the eigenbasis e of
the operator X the centralizer becomes a diagonal torus with typical element
z = diag(z1 , z2 , . . . , zr ) and the character εi : z 7→ zi .
Let now X ν = ϕa (X), Dν = Z(X ν ), and eT be the weight basis of Hν , introduced in the beginning of this section, parameterized by semistandard tableaux
T of shape ν and arranged in the order of eigenvalues aν . Then the characQ
ν
ter of the pull back ϕ−1
a (Ek ) is just the weight
i∈T εi of the k-th vector eT ,
where the tableau T is determined from the equation aνk = aT , cf. (3.7). Thus

## Practical algorithm [thesis 977-1076]
Practical algorithm

Note that, the set of all allowed spectra of the electron density matrix forms a
convex polytope, called Moment Polytope. The above theorem gives an inner
approximation to this polytope, while the Theorem 1.4.1 gives an outer approximation. Combining these two results leads to the following practical algorithm,
which allows to find explicit constraints for the N -representability problem:
1. Find all irreducible components Hλ ⊂ S m (∧N Hr ) for m ≤ M .
e which gives an inner
2. Calculate the convex hull of the corresponding spectra λ
in
⊂ P for the moment polytope P.
approximation PM
in
3. Identify the facets of PM
that are given by the inequalities of Theorem 3.2.1.
out
⊃ P.
They cut out an outer approximation PM
in
out
4. Increase M and continue until PM
= PM
.

1.6

Taking into account spin

Actually, the state space of a single particle with spin splits into the tensor
product H = Hr ⊗ Hs of the orbital component Hr and the spin component Hs .
The total N -fermion space decomposes into spin-orbital components as follows
[38]
∧N (Hr ⊗ Hs ) =

X

t

Hrν ⊗ Hsν ,

(1.11)

|ν|=N
t

where ν t stands for the transpose diagram, and Hrν and Hsν are irreducible representations of unitary groups U (Hr ) and U (Hs ) with Young diagrams ν and ν t ,

14

CHAPTER 1. INTRODUCTION

respectively. In many physical systems, like electrons in an atom or a molecule,
the total spin is a well defined quantity which singles out a specific component of
t

this decomposition. We have to deal with pure states ψ ∈ Hrν ⊗Hsν . From n◦ 1.1.5,
t

t

the reduced states ρνr and ρνs of ψ are isospectral, that is, Specρνr = Specρνs . So
t

we can identify the spectrum Specρνs with Specρνr .
On the other hand, the Schur-Weyl duality
H⊗N =

X

Hν ⊗ V ν ,

(1.12)

|ν|=N

between irreducible representations Hν and V ν of the unitary U(H) and the
symmetric SN groups, respectively, allows to define the ith reduced density matrix
ρi : H for ρν : Hν as the reduced density matrix for ρν ⊗ 1. The operator ρν ⊗ 1
acting on the component Hν ⊗ V ν commutes with SN , and hence the reduced
state ρi is independent of i.
The problem which we address here is the following:
What are the constraints on the spectra of ρν : Hν and its particle
density matrix ρ = N ρi : H?
It is a variation of the N -representability problem. We call it as ν-representability.
As a result, by solving the above problem we may find all constraints on the
t

spectra Specρνr and Specρνs .
In Chapter 3 we give the formal solution of this problem which is the generalization of Theorem 1.4.1, and in Chapter 5 we give another solution which is
the generalization of Theorem 1.5.1. Combining these two approaches gives an
algorithm which is the modified version of the one given in previous section. In
this new algorithm there is a small modification: instead of the symmetric power
S m (∧N H) we have to use the plethysm [Hν ]µ .
As an example let’s consider the constraints on the mixed state ρν and its
reduced matrix ρ of a system of three electrons of the total spin J = 1/2. The
problem is equivalent to ν-representability for ν =

and Spec ρν = (µ1 , µ2 ). A

calculation based on the above algorithm shows that the constraints amount to

## Implementation plan
1. PURE-STATE SPECIALIZATION: mu=(1,0,...,0) => every constraint: sum_i a_i*l_{v(i)} <= a_T (one N-subset T)
2. CUBICLE EDGES: extremal rays of {a1>=...>=ar>=0} cut by hyperplanes a_T=a_T' — small-integer edge vectors
3. SCHUBERT TEST c_vw(a): pullback phi*_a (substitute x^nu_k -> x_{T(k)}), expand sigma_w in Schubert basis
   of Z[x1..xr]/(e1..er), read coefficient of sigma_v. Tooling: lrcalc (has schubmult / flag Schubert ops).
4. Validity: c != 0; AK/Ressayre sharpness: c = 1 for the complete irredundant list.
5. TEST-FIRST GATES: reproduce (3,6): BD equalities + l4<=l5+l6; (3,7): 4; (3,8): 31; (3,9): 52; (3,10): 93.
6. Sandbox tooling confirmed working: lrcalc (pip), lrs (apt). Validation suite: padding/lift/duality/census
   modules from this week — every generator output must pass before use.

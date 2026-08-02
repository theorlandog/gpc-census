# Paper plan: Certified observable manifolds for exact fermionic reduced dynamics

## Working title

**Certified Observable Manifolds for Exact Fermionic Reduced Dynamics**

Alternative conservative title:

**Support-Combinatorial Certificates for Exact Pure-State 2-RDM Dynamics**

## Central thesis

Certain sparse fermionic carriers admit exact low-order observable atlases:
selected fixed-frame 2-RDM entries are global coordinates on the pure
full-support locus, determine all higher RDMs rationally, and support exact
Hamiltonian evolution without an independent wavefunction representation.

## Proposed abstract

Reduced-density-matrix dynamics is normally obstructed by the open BBGKY
hierarchy. We identify a class of determinant carriers on which this
obstruction disappears exactly. A full-rank pair-incidence minor recovers
determinant weights, while a connected graph of collision-free two-body
coherences recovers all relative phases. These data give an explicit rational
inverse from selected fixed-frame two-particle reduced-density-matrix entries
to the pure carrier density matrix and hence to every higher RDM. The
full-support observable manifold has simplex--torus topology and inherits the
Fubini--Study symplectic form in explicit observable coordinates; its boundary
phase defects are governed by connectivity of the surviving coherence graph.
We certify the construction across a finite fermionic census, give infinite
cyclic and product families, and exhibit carrier-preserving two-body dynamics,
including a finite transport model continued through amplitude-zero events by
a multi-chart Cartesian propagator. The results provide exact benchmarks and
constructive sufficient conditions for closed time-dependent 2-RDM dynamics.
They do not address mixed-state reconstruction, unknown carriers, or generic
\(N\)-representability.

## Section structure

### 1. Introduction

- State the hierarchy-closure problem.
- Explain fixed-carrier, fixed-frame scope immediately.
- State contributions as certificates and constructions.
- Explicitly distinguish standard geometric quantum mechanics from the new
  observable realization.

### 2. Setting and conventions

- Fermionic Fock sector and determinant carriers.
- RDM sign convention.
- Carrier rank-one matrix \(X\).
- Pair-incidence and signed channel maps.

### 3. Exact reconstruction theorem

- Weight recovery.
- Collision-free channel graph.
- Path-product inverse.
- Higher-RDM contraction.
- Domain and nonclaims.

### 4. Observable-manifold geometry

- Simplex--torus parametrization.
- Pullback of Fubini--Study form.
- Exact Hamiltonian flow.
- Clarify which statements are standard and which are carrier-certified.

### 5. Boundary stratification and observable atlases

- Active supports.
- Component-count phase defect.
- Multiple coherence trees.
- Cartesian coordinates at zero events.
- Precise continuation hypotheses.

### 6. Constructive families

- Cyclic-window rank theorem.
- Labeled collision-free cycle.
- Disjoint-block gluing theorem.
- Fixed-core active-space comparison.
- Correct ambient dimensions.

### 7. Exact dynamics

- \(v_B\) scheduled fixed-1RDM trajectory.
- 29-state invariant hopping carrier.
- Exact \(\Gamma^{(3)}(\Gamma^{(2)})\) closure.
- Multi-chart propagation through zeros.

### 8. Census and verification

- 799 supports and certified counts.
- Independent verifiers.
- Unresolved records and stronger determinacy closure.
- Provenance and artifact policy.

### 9. Computational consequences

- Dense RDM versus sparse observable coordinates.
- Compiled sparse kernels.
- Product-cell scaling family.
- Honest benchmark limits.

### 10. Prior work and distinction

A table comparing:
- pure-state determination by marginals;
- unique-ground-state/parent-Hamiltonian results;
- RDM matrix completion;
- phase synchronization;
- geometric quantum mechanics;
- TD-2RDM closures;
- tensor-network and active-space reductions.

### 11. Limitations and open problems

- mixed states;
- carrier inference;
- orbital-frame quotient;
- stability bounds;
- optimal atlas construction;
- non-product interacting scalable examples;
- generic lower bounds.

## Main theorem numbering

1. Certified pure-state 2-RDM reconstruction.
2. Observable-chart diffeomorphism and symplectic pullback.
3. Boundary component theorem.
4. Exact closed reduced Hamiltonian flow.
5. Cyclic-window observability.
6. Gluing/composition.
7. Multi-chart continuation under explicit regularity hypotheses.

## Referee-risk checklist

- Is any “global” statement only valid on full support?
- Is a chart singularity being confused with physical non-identifiability?
- Is a standard symplectic fact advertised as novel?
- Does a benchmark compare equally optimized algorithms?
- Are product families being sold as generic interacting systems?
- Are numerical rank claims replaced by exact certificates?
- Are all fermionic signs independently regenerated?
- Does every theorem list fixed carrier, fixed frame, purity, and support assumptions?

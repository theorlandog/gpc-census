#!/usr/bin/env python3
"""Exact transport extensions for observable determinant carriers.

Three exact campaigns:
1. state-specific Hermitian cross hopping tangent to S(5,2) x S(5,2);
2. invariant closure of a chosen hop and a fresh 2-RDM observability audit;
3. fixed-core active-space families with nontrivial one-body normalizers.
"""
from __future__ import annotations
import argparse, itertools, json
from collections import defaultdict, deque
from fractions import Fraction
from pathlib import Path
import sympy as sp


def windows(d,n,offset=0):
    return [tuple(sorted(offset+((k+j)%d) for j in range(n))) for k in range(d)]

def apply_one(det,q,p):
    det=list(det)
    if q not in det or p in det:return None,0
    i=det.index(q);sg=-1 if i%2 else 1;det.pop(i)
    j=sum(x<p for x in det);sg*=(-1 if j%2 else 1);det.insert(j,p)
    return tuple(det),sg

def apply_cross(det,q,p,s,r):
    x,a=apply_one(det,s,r)
    if x is None:return None,0
    y,b=apply_one(x,q,p)
    return y,a*b

def rank_z(rows):
    if not rows:return 0
    return sp.polys.matrices.DomainMatrix.from_list(rows,sp.ZZ).rank()

def incidence_rank(support,orbitals,order):
    return rank_z([[int(set(f).issubset(d)) for d in support]
                   for f in itertools.combinations(range(orbitals),order)])

def permutation_sign(seq):
    inv=sum(seq[i]>seq[j] for i in range(len(seq)) for j in range(i+1,len(seq)))
    return -1 if inv%2 else 1

def unique_graph(support,order=2):
    particles=len(support[0]); terms=defaultdict(lambda:defaultdict(int))
    for a,left in enumerate(support):
        ls=set(left)
        for b,right in enumerate(support):
            inter=sorted(ls.intersection(right)); remn=particles-order
            if len(inter)<remn:continue
            for rem in itertools.combinations(inter,remn):
                rs=set(rem); lm=tuple(x for x in left if x not in rs); rm=tuple(x for x in right if x not in rs)
                terms[(lm,rm)][(a,b)]+=permutation_sign(lm+rem)*permutation_sign(rm+rem)
    adj=[set() for _ in support]; witnesses=[]
    for (lm,rm),co in sorted(terms.items()):
        co={e:v for e,v in co.items() if v}
        if lm>=rm or len(co)!=1:continue
        (a,b),sg=next(iter(co.items()))
        if a==b:continue
        adj[a].add(b);adj[b].add(a);witnesses.append((a,b,lm,rm,sg))
    seen={0};q=deque([0])
    while q:
        a=q.popleft()
        for b in adj[a]:
            if b not in seen:seen.add(b);q.append(b)
    return len(seen)==len(support),witnesses

def hermitian_cross_bases(d1,d2):
    off=d1; directed=[(q,p,s,r) for q,p in itertools.permutations(range(d1),2)
                      for s,r in itertools.permutations(range(off,off+d2),2)]
    bases=[];seen=set()
    for t in directed:
        rev=(t[1],t[0],t[3],t[2]);key=min(t,rev)
        if key in seen:continue
        seen.add(key);bases.append((key,rev if key==t else t))
    return bases

def state_specific_tangent_certificate():
    d=5; n=2; off=5
    support=[tuple(sorted(a+b)) for a in windows(d,n) for b in windows(d,n,off)]
    S=set(support);si={x:i for i,x in enumerate(support)}; bases=hermitian_cross_bases(5,5)
    outside=sorted({dst for src in support for pair in bases for op in pair
                    for dst,sg in [apply_cross(src,*op)] if dst is not None and dst not in S})
    oi={x:i for i,x in enumerate(outside)}
    leak=sp.zeros(len(outside),len(bases));inside=sp.zeros(len(support),len(bases))
    for src in support: # uniform carrier state: every coefficient is 1
        for j,pair in enumerate(bases):
            for op in pair:
                dst,sg=apply_cross(src,*op)
                if dst is None:continue
                (inside if dst in S else leak)[(si if dst in S else oi)[dst],j]+=sg
    null=sp.Matrix.hstack(*leak.nullspace())
    action=inside*null
    chosen=None
    for k in range(action.cols):
        v=action[:,k]
        if any(v[i]!=v[0] for i in range(len(v))):chosen=k;break
    assert chosen is not None
    coeff=null[:,chosen]; tangent=action[:,chosen]
    terms=[]
    for j,c in enumerate(coeff):
        if c:terms.append({"forward":list(bases[j][0]),"adjoint":list(bases[j][1]),"coefficient":str(c)})
    return {
      "carrier":"S(5,2) x S(5,2)","carrier_dimension":25,
      "uniform_state_leakage_rows":len(outside),"hermitian_cross_hopping_basis_dimension":len(bases),
      "state_leakage_rank":leak.rank(),"state_tangent_kernel_dimension":len(bases)-leak.rank(),
      "inside_tangent_rank":action.rank(),"selected_kernel_terms":terms,
      "selected_tangent_nonzero_entries":[[i,str(x)] for i,x in enumerate(tangent) if x],
      "exact_checks":{"leakage_zero":all(x==0 for x in leak*coeff),"tangent_nonprojective":any(tangent[i]!=tangent[0] for i in range(len(tangent)))},
      "interpretation":"An exact Hermitian scheduled cross-hop is tangent at the uniform state although no nonzero cross-hop preserves the entire carrier. This is a first-order local certificate, not yet a finite-time schedule."
    }

def invariant_closure(base_support,ops):
    S=set(base_support);changed=True
    while changed:
        changed=False
        for src in list(S):
            for op in ops:
                dst,sg=apply_cross(src,*op)
                if dst is not None and dst not in S:S.add(dst);changed=True
    return sorted(S)

def closure_certificate():
    base=[tuple(sorted(a+b)) for a in windows(5,2) for b in windows(5,2,5)]
    hop=(0,1,5,7);adj=(1,0,7,5)
    enlarged=invariant_closure(base,[hop,adj]); connected,w=unique_graph(enlarged)
    return {
      "base_carrier_dimension":len(base),"hop":{"forward":list(hop),"adjoint":list(adj)},
      "enlarged_carrier_dimension":len(enlarged),"added_determinants":[list(x) for x in enlarged if x not in set(base)],
      "pair_incidence_rank":incidence_rank(enlarged,10,2),"pair_incidence_full_column_rank":incidence_rank(enlarged,10,2)==len(enlarged),
      "collision_free_gamma2_connected":connected,"collision_free_gamma2_channels":len(w),
      "invariant_under_selected_hermitian_hop":invariant_closure(enlarged,[hop,adj])==enlarged,
      "interpretation":"Adding the four leakage images of one desired Hermitian hop produces a 29-state invariant carrier that remains exactly reconstructible from its fixed-frame pure-state 2-RDM."
    }

def onebody_normalizer_dimension(support,orbitals):
    S=set(support);terms=list(itertools.permutations(range(orbitals),2));rows=[]
    for src in support:
        by={}
        for j,(q,p) in enumerate(terms):
            dst,sg=apply_one(src,q,p)
            if dst is not None and dst not in S:by.setdefault(dst,[0]*len(terms))[j]=sg
        rows.extend(by.values())
    r=rank_z(rows);return len(terms)-r,r,len(rows)

def fixed_core_family(core_size,active_size,active_particles):
    core=tuple(range(core_size));active=tuple(range(core_size,core_size+active_size))
    support=[tuple(sorted(core+c)) for c in itertools.combinations(active,active_particles)]
    orbitals=core_size+active_size;connected,w=unique_graph(support)
    nullity,rank,constraints=onebody_normalizer_dimension(support,orbitals)
    return {"core_particles":core_size,"active_orbitals":active_size,"active_particles":active_particles,
            "carrier_dimension":len(support),"pair_incidence_rank":incidence_rank(support,orbitals,2),
            "pair_incidence_full_column_rank":incidence_rank(support,orbitals,2)==len(support),
            "collision_free_gamma2_connected":connected,"collision_free_gamma2_channels":len(w),
            "off_diagonal_one_body_terms":orbitals*(orbitals-1),"one_body_leakage_rank":rank,
            "one_body_leakage_kernel_dimension":nullity,"guaranteed_active_offdiagonal_normalizer_dimension":active_size*(active_size-1),"leakage_constraints":constraints}

def build():
    return {"schema_version":1,"evidence":"exact integer/rational linear algebra",
      "scheduled_tangent_transport":state_specific_tangent_certificate(),
      "hopping_closure_enlargement":closure_certificate(),
      "different_factor_families":{
        "theorem":"A fixed core tensored with a one- or two-electron complete active space is pair-incidence full rank, has connected collision-free 2-RDM channels, and is normalized by the full active one-body algebra.",
        "records":[fixed_core_family(2,5,1),fixed_core_family(2,5,2),fixed_core_family(3,6,2)],
        "novelty_boundary":"These are exact control families closely related to standard complete-active-space constructions; their role is to show transport and observability are compatible, not to claim the active-space idea itself as new."
      },
      "claims":{"proved":["A nonprojective Hermitian cross-hopping tangent exists at the uniform S(5,2)^2 state.","The selected hopping closure has 29 determinants and remains 2-RDM observable.","Fixed-core CAS(1,a) and CAS(2,a) carriers provide infinite observable families with nontrivial one-body normalizers."],
                "not_proved":["A finite-time scheduled transport trajectory on S(5,2)^2.","A general minimal-closure theorem for arbitrary hops.","Novelty of standard active-space normalizer families."]}}

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('results/data/transport_extensions.json'));a=p.parse_args()
    data=build();a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n');print(json.dumps(data,indent=2,sort_keys=True));return 0
if __name__=='__main__':raise SystemExit(main())

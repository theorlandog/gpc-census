#!/usr/bin/env python3
"""Exact carrier-preserving cross-block pair-hopping audit.

Audits cyclic-window product carriers for two-body terms
(a_p^dag a_q)(a_r^dag a_s) that preserve particle number in each block.
The leakage matrix is integer-valued; its exact rank gives the dimension of
all carrier-preserving cross-block hopping kernels.
"""
from __future__ import annotations
import argparse,itertools,json
from fractions import Fraction
import sympy as sp
from pathlib import Path


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

def rank_q(rows):
    if not rows:return 0
    return sp.polys.matrices.DomainMatrix.from_list(rows, sp.ZZ).rank()

def incidence_rank(support,orbitals,order=2):
    rows=[]
    for f in itertools.combinations(range(orbitals),order):
        fs=set(f);rows.append([int(fs.issubset(d)) for d in support])
    return rank_q(rows)

def audit(d1,n1,d2,n2):
    off=d1; A=windows(d1,n1);B=windows(d2,n2,off)
    support=[tuple(sorted(a+b)) for a in A for b in B]; S=set(support)
    full=[tuple(sorted(a+b)) for a in itertools.combinations(range(d1),n1)
          for b in itertools.combinations(range(off,off+d2),n2)]
    terms=[(q,p,s,r) for q,p in itertools.permutations(range(d1),2)
           for s,r in itertools.permutations(range(off,off+d2),2)]
    leakage=[]
    for src in support:
        by_dst={}
        for j,t in enumerate(terms):
            dst,sg=apply_cross(src,*t)
            if dst is not None and dst not in S:
                by_dst.setdefault(dst,[0]*len(terms))[j]=sg
        leakage.extend(by_dst.values())
    rank=rank_q(leakage); nullity=len(terms)-rank
    return {"factors":[[d1,n1],[d2,n2]],"carrier_dimension":len(support),
            "pair_incidence_rank":incidence_rank(support,d1+d2),
            "pair_incidence_full_column_rank":incidence_rank(support,d1+d2)==len(support),
            "cross_hopping_terms":len(terms),"leakage_constraints":len(leakage),
            "leakage_rank":rank,"carrier_preserving_hopping_dimension":nullity}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=Path('results/data/transport_rigidity.json'))
    args=ap.parse_args()
    cases=[(4,2,4,2),(4,2,5,2),(5,2,5,2),(5,3,5,3),(6,2,6,2),(6,3,6,3),(7,2,7,2),(7,3,7,3)]
    data={"schema_version":1,"scope":"exact integer leakage-rank audit",
          "records":[audit(*c) for c in cases],
          "claims":{"proved_by_certificate":[
              "For S(5,2) x S(5,2), every cross-block pair-hopping kernel preserving the full carrier is zero.",
              "S(4,2) x S(4,2) has a 16-dimensional preserving hopping space but pair incidence rank 11<16."
          ],"nonclaim":"The finite scan is not a general no-go theorem for all cyclic parameters."}}
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Exact lower-rank regression prototype for GPC facet generation.

Phase 0 deliberately separates:
  1. exact facet recovery from a certified vertex oracle;
  2. Bruhat/Gale signatures and inherited-facet labels;
  3. the carrier selected by saturation of a recovered GPC;
  4. exact 2-RDM observability diagnostics on that carrier.

The regression target is the Borland--Dennis system wedge^3 C^6.
No floating-point convex hull routines are used.
"""
from __future__ import annotations
import itertools, json, math
from pathlib import Path
import sympy as sp

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/hybrid_rank_generator_prototype.json"

# Reduced coordinates:
# lambda=(1-z,1-y,1-x,x,y,z), x=lambda4,y=lambda5,z=lambda6.
BD_VERTICES = [
    (sp.Rational(0), sp.Rational(0), sp.Rational(0)),
    (sp.Rational(1,2), sp.Rational(1,2), sp.Rational(0)),
    (sp.Rational(1,2), sp.Rational(1,4), sp.Rational(1,4)),
    (sp.Rational(1,2), sp.Rational(1,2), sp.Rational(1,2)),
]

def primitive_inequality(normal, offset):
    vals=list(normal)+[-offset]
    den=sp.ilcm(*[sp.denom(v) for v in vals])
    ints=[int(v*den) for v in vals]
    g=0
    for v in ints: g=math.gcd(g,abs(v))
    ints=[v//g for v in ints]
    # inequality n.x <= b, so final integer is a.x <= b
    a=tuple(ints[:-1]); b=-ints[-1]
    first=next((v for v in a+(b,) if v),1)
    if first<0:
        a=tuple(-v for v in a); b=-b
        # orientation will be restored by vertex test later
    return a,b

def recover_facets(vertices):
    """Recover irredundant tetrahedral facets exactly."""
    V=[sp.Matrix(v) for v in vertices]
    out=[]
    for face in itertools.combinations(range(len(V)),3):
        p,q,r=[V[i] for i in face]
        u=q-p; v=r-p
        n=sp.Matrix([
            u[1]*v[2]-u[2]*v[1],
            u[2]*v[0]-u[0]*v[2],
            u[0]*v[1]-u[1]*v[0],
        ])
        if n==sp.zeros(3,1): continue
        b=n.dot(p)
        residual=[sp.simplify(n.dot(x)-b) for x in V]
        if all(t<=0 for t in residual):
            pass
        elif all(t>=0 for t in residual):
            n=-n; b=-b; residual=[-t for t in residual]
        else:
            continue
        a,bb=primitive_inequality(tuple(n),b)
        # restore <= orientation after primitive sign normalization
        if not all(sum(a[j]*vertices[i][j] for j in range(3))<=bb
                   for i in range(len(vertices))):
            a=tuple(-x for x in a);bb=-bb
        out.append({"a":list(a),"b":bb,"face":list(face)})
    # canonical unique
    unique={}
    for row in out: unique[(tuple(row["a"]),row["b"])]=row
    return sorted(unique.values(),key=lambda r:(r["a"],r["b"]))

def lift_lambda(v):
    x,y,z=v
    return (1-z,1-y,1-x,x,y,z)

def facet_meaning(a,b):
    a=tuple(a)
    names={
      (0,-1,1): "lambda6 <= lambda5 (Weyl chamber)",
      (-1,1,0): "lambda5 <= lambda4 (Weyl chamber)",
      (1,-1,-1): "lambda4 <= lambda5 + lambda6 (Borland-Dennis GPC)",
      (2,0,0): "2 lambda4 <= 1, i.e. lambda4 <= 1/2 (ordering plus pair equalities)",
    }
    return names.get(a,f"{a} <= {b}")

def bd_selection_carrier():
    # Saturation of D=2-(lambda1+lambda2+lambda4)>=0:
    # kernel of Dhat=2-(n1+n2+n4), zero-based orbitals {0,1,3}.
    return [D for D in itertools.combinations(range(6),3)
            if 2-sum(p in {0,1,3} for p in D)==0]

def incidence_rank(S,k):
    features=list(itertools.combinations(range(6),k))
    M=sp.Matrix([[int(set(P)<=set(D)) for D in S] for P in features])
    return int(M.rank())

def bitstate(D):
    s=0
    for p in D:s|=1<<p
    return s

def ann(state,p):
    if not ((state>>p)&1): return None
    sign=-1 if (state&((1<<p)-1)).bit_count()%2 else 1
    return state^(1<<p),sign

def cre(state,p):
    if (state>>p)&1: return None
    sign=-1 if (state&((1<<p)-1)).bit_count()%2 else 1
    return state|(1<<p),sign

def matrix_element(left,right,P,Q):
    state=bitstate(right);sign=1
    for q in Q:
        out=ann(state,q)
        if out is None:return 0
        state,s=out;sign*=s
    for p in reversed(P):
        out=cre(state,p)
        if out is None:return 0
        state,s=out;sign*=s
    return sign if state==bitstate(left) else 0

def collision_free_graph(S,k=2):
    features=list(itertools.combinations(range(6),k))
    edges=set(); collided=0
    for P in features:
        for Q in features:
            if P==Q:continue
            terms=[]
            for a,b in itertools.combinations(range(len(S)),2):
                if (matrix_element(S[a],S[b],P,Q) or
                    matrix_element(S[b],S[a],P,Q)):
                    terms.append((a,b))
            if len(terms)==1:edges.add(terms[0])
            elif len(terms)>1:collided+=1
    adjacency={i:set() for i in range(len(S))}
    for a,b in edges:
        adjacency[a].add(b);adjacency[b].add(a)
    seen=set()
    if S:
        stack=[0]
        while stack:
            u=stack.pop()
            if u in seen:continue
            seen.add(u);stack.extend(adjacency[u]-seen)
    return {
      "edge_count":len(edges),
      "connected":len(seen)==len(S),
      "collided_oriented_channel_count":collided,
    }

def main(out=None):
    facets=recover_facets(BD_VERTICES)
    for f in facets:f["meaning"]=facet_meaning(f["a"],f["b"])
    S=bd_selection_carrier()
    r1=incidence_rank(S,1);r2=incidence_rank(S,2)
    graph=collision_free_graph(S,2)
    payload={
      "schema_version":1,
      "system":"wedge^3 C^6",
      "prototype_scope":
        "exact V-to-H regression plus GPC-selected-carrier 2-RDM audit",
      "vertices_reduced":[[str(x) for x in v] for v in BD_VERTICES],
      "vertices_lambda":[[str(x) for x in lift_lambda(v)]
                         for v in BD_VERTICES],
      "recovered_facets":facets,
      "nontrivial_gpc_recovered":
        any(f["a"]==[1,-1,-1] and f["b"]==0 for f in facets),
      "pinned_carrier":{
        "constraint":"2-(lambda1+lambda2+lambda4)>=0",
        "support":[list(D) for D in S],
        "size":len(S),
        "singleton_rank":r1,
        "non_toric_phase_complexity":len(S)-r1,
        "pair_incidence_rank":r2,
        "gamma2_graph":graph,
        "pure_full_support_gamma2_certified":
            r2==len(S) and graph["connected"],
      },
      "limitations":[
        "This phase-0 prototype starts from a certified exact vertex oracle.",
        "It does not yet enumerate all Ressayre elements from representation weights.",
        "Bruhat and Schubert modules are next-stage candidate and cross-certificate engines.",
      ],
    }
    text=json.dumps(payload,indent=2)+"\n"
    if out:
        Path(out).parent.mkdir(parents=True,exist_ok=True)
        Path(out).write_text(text)
        print(f"wrote {out}: {len(facets)} facets, carrier of size {len(S)}")
    else:
        print(text)
    return payload

if __name__=="__main__":
    import argparse
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,default=DEFAULT_OUT)
    args=ap.parse_args()
    main(args.out)

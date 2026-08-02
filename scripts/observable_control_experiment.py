import argparse, itertools, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/observable_control_experiment.json"


def dets(d,n): return list(itertools.combinations(range(d),n))
def dist(a,b): return len(set(a)-set(b))

def components(m, edges):
    adj=[set() for _ in range(m)]
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    seen=set(); cs=[]
    for s in range(m):
        if s in seen: continue
        q=[s]; seen.add(s); c=[]
        while q:
            u=q.pop(); c.append(u)
            for v in adj[u]:
                if v not in seen: seen.add(v); q.append(v)
        cs.append(c)
    return cs

def control_q(S,n):
    m=len(S)
    for q in range(1,n+1):
        E={(i,j) for i,j in itertools.combinations(range(m),2) if dist(S[i],S[j])<=q}
        if len(components(m,E))==1: return q
    return n

def unique_edges(S,n,q):
    # A q-RDM channel for pair (a,b) exists iff dist<=q. Collision labels are
    # (removed from a, added to b, common remainder chosen); use exact (P,Q)
    # channel identity from differing modes padded by common subsets.
    channels=defaultdict(set)
    for a,b in itertools.combinations(range(len(S)),2):
        A,B=set(S[a]),set(S[b]); rem=tuple(sorted(A-B)); add=tuple(sorted(B-A)); r=len(rem)
        if r>q: continue
        common=sorted(A&B)
        for pad in itertools.combinations(common,q-r):
            P=tuple(sorted(rem+pad)); Q=tuple(sorted(add+pad))
            channels[(P,Q)].add((a,b))
            channels[(Q,P)].add((a,b))
    return {next(iter(es)) for es in channels.values() if len(es)==1}, sum(len(es)>1 for es in channels.values())

def observe_q(S,n):
    m=len(S)
    for q in range(1,n+1):
        E,_=unique_edges(S,n,q)
        if len(components(m,E))==1: return q
    return n

def random_dynamics(S,q,seed=0,t=1.0):
    rng=np.random.default_rng(seed); m=len(S)
    H=np.diag(rng.normal(size=m)).astype(complex)
    for i,j in itertools.combinations(range(m),2):
        if dist(S[i],S[j])<=q:
            z=rng.normal()+1j*rng.normal(); H[i,j]=z; H[j,i]=z.conjugate()
    # exact diagonalization
    w,V=np.linalg.eigh(H)
    psi=np.ones(m,dtype=complex)/math.sqrt(m)
    out=V@(np.exp(-1j*w*t)*(V.conj().T@psi))
    p=np.abs(out)**2
    entropy=float(-(p[p>1e-15]*np.log(p[p>1e-15])).sum())
    participation=float(1/(p@p))
    return entropy, participation

def run(d=6,n=3,max_m=6,samples_per_size=400,seed=1):
    rng=random.Random(seed); universe=dets(d,n); rows=[]
    for m in range(2,max_m+1):
        combos=list(itertools.combinations(universe,m))
        if len(combos)>samples_per_size: combos=rng.sample(combos,samples_per_size)
        for idx,S in enumerate(combos):
            qc=control_q(S,n); qo=observe_q(S,n)
            coll={q:unique_edges(S,n,q)[1] for q in range(1,n+1)}
            ent,pr=random_dynamics(S,qc,seed=seed*100000+m*1000+idx)
            rows.append({'m':m,'q_control':qc,'q_observe':qo,'gap':qo-qc,'collisions':coll,'entropy':ent,'participation':pr})
    summary={
      'systems':{'d':d,'N':n},'supports':len(rows),
      'q_equal':sum(r['q_control']==r['q_observe'] for r in rows),
      'q_observe_gt_control':sum(r['q_observe']>r['q_control'] for r in rows),
      'q_observe_lt_control':sum(r['q_observe']<r['q_control'] for r in rows),
      'gap_histogram':dict(Counter(str(r['gap']) for r in rows)),
      'pair_histogram':dict(Counter(f"{r['q_control']},{r['q_observe']}" for r in rows)),
      'mean_entropy_by_q_control':{},'mean_participation_by_q_control':{}
    }
    for q in range(1,n+1):
      rr=[r for r in rows if r['q_control']==q]
      if rr:
        summary['mean_entropy_by_q_control'][str(q)]=sum(r['entropy'] for r in rr)/len(rr)
        summary['mean_participation_by_q_control'][str(q)]=sum(r['participation'] for r in rr)/len(rr)
    return {'schema_version':1,'summary':summary,'rows':rows}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=DEFAULT_OUT)
    ap.add_argument('--samples-per-size',type=int,default=400)
    ap.add_argument('--seed',type=int,default=1)
    args=ap.parse_args()
    out=run(samples_per_size=args.samples_per_size,seed=args.seed)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(f"wrote {args.out}: {out['summary']['supports']} carriers")


if __name__=='__main__':
    main()

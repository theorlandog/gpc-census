"""v8s: double-block real(+-)-exclusion sweep with SHARDING + CHECKPOINTS.
Reconstruction of interference8.py, validated against known counts, plus:
  --shard s S : partition search space by first-w inclusion bits (S=2^w shards)
  i8_done.txt : ledger "idx.shard", resumable across reboots/spot kills
Usage:
  python interference8s.py -1                    # list ansatze
  python interference8s.py <idx> [size]          # one ansatz, no shard
  python interference8s.py --worker W K [S]      # stride workers over (idx,shard)
"""
import sys, os
from itertools import combinations
from ortools.sat.python import cp_model

dets=list(combinations(range(9),4)); nd=len(dets); index={T:t for t,T in enumerate(dets)}
DONE_FILE="i8_done.txt"
def _done_set():
    if not os.path.exists(DONE_FILE): return set()
    return {l.strip() for l in open(DONE_FILE) if l.strip()}
def _mark_done(key):
    with open(DONE_FILE,"a") as f: f.write(key+"\n"); f.flush(); os.fsync(f.fileno())

def hop_sign(T,i,j):
    lst=list(T); sj=(-1)**lst.index(j); lst.remove(j)
    pos=sum(1 for k in lst if k<i); si=(-1)**pos
    return si*sj, tuple(sorted(lst+[i]))
def sqfree(x):
    m=1; f=x; d=2
    while d*d<=f:
        while f%(d*d)==0: f//=d*d; m*=d
        d+=1
    return m,f
def splits(e1,e2):
    s=e1+e2; out=[]
    for a in range(1,s//2+1):
        b=s-a; x2=a*b-e1*e2
        if x2>0: out.append((a,b,x2))
    return out
def all_ansatze():
    M=[20,14,14,14,14,4,4,4,4]
    pair_types=[(20,14),(20,4),(14,14),(14,4),(4,4)]
    out=[]
    for i in range(len(pair_types)):
        for j in range(i,len(pair_types)):
            p1,p2=pair_types[i],pair_types[j]
            need=list(p1)+list(p2); pool=list(M); ok=True
            for x in need:
                if x in pool: pool.remove(x)
                else: ok=False; break
            if not ok: continue
            left=sorted(pool,reverse=True)
            s1,s2=splits(*p1),splits(*p2)
            for si_,(a1,b1,x1) in enumerate(s1):
                for sj_,(a2,b2,x2) in enumerate(s2):
                    if i==j and sj_<si_: continue   # same-type mirror dedup
                    d=left+[a1,b1,a2,b2]
                    out.append((p1,p2,a1,b1,x1,a2,b2,x2,d))
    return out
ANS=all_ansatze()

def block_ok(kw,E,mx,fx):
    """exists edge signs with per-sqfree-class sums: fx-class = +-mx, others = 0"""
    if not E: return False
    cls={}
    for a,b,sg in E:
        m,f=sqfree(kw[a]*kw[b]); cls.setdefault(f,[]).append(m)
    for f,ms in cls.items():
        reach={0}
        for m in ms: reach={r+m for r in reach}|{r-m for r in reach}
        if f==fx:
            if mx not in reach and -mx not in reach: return False
        else:
            if 0 not in reach: return False
    if fx not in cls: return False
    return True

def _run(idx, size=None, shard=None, S=1):
    key=f"{idx}.{shard if shard is not None else 'all'}"
    if size is None and key in _done_set():
        print(f"status: SKIP-DONE  ansatz: {idx} shard: {shard}",flush=True); return False
    p1,p2,a1,b1,x1,a2,b2,x2,d=ANS[idx]
    U1,V1,U2,V2=5,6,7,8
    mx1,fx1=sqfree(x1); mx2,fx2=sqfree(x2)
    m=cp_model.CpModel()
    k=[m.NewIntVar(0,23,f"k{t}") for t in range(nd)]
    y=[m.NewBoolVar(f"y{t}") for t in range(nd)]
    for t in range(nd):
        m.Add(k[t]<=23*y[t]); m.Add(k[t]>=y[t])
    for mo in range(9):
        m.Add(sum(k[t] for t in range(nd) if mo in dets[t])==d[mo])
    allowed={frozenset({U1,V1}),frozenset({U2,V2})}
    for p in range(nd):
        Tp=set(dets[p])
        for q in range(p+1,nd):
            Tq=set(dets[q])
            if len(Tp&Tq)==3 and frozenset(Tp^Tq) not in allowed:
                m.AddBoolOr([y[p].Not(),y[q].Not()])
    if size is not None: m.Add(sum(y)==size)
    if shard is not None:                     # bitmask partition on w likely-included dets
        w=S.bit_length()-1
        blockmodes={U1,V1,U2,V2}
        order=sorted(range(nd), key=lambda t: (-len(set(dets[t])&blockmodes), t))
        sd=order[:w]
        for i in range(w):
            m.Add(y[sd[i]]== (shard>>i)&1)
    def edgecat(u,v):
        E=[]
        for t in range(nd):
            T=dets[t]
            if v in T and u not in T:
                s,Tp=hop_sign(T,u,v); tp=index.get(Tp)
                if tp is not None: E.append((tp,t,s))
        return E
    E1,E2=edgecat(U1,V1),edgecat(U2,V2)
    class CB(cp_model.CpSolverSolutionCallback):
        def __init__(s_):
            super().__init__(); s_.n=0; s_.hit=None
        def on_solution_callback(s_):
            s_.n+=1
            if s_.n%25000==0: print(f"  ...{s_.n} solutions",flush=True)
            kw={t:s_.Value(k[t]) for t in range(nd) if s_.Value(y[t])}
            e1=[(a,b,sg) for a,b,sg in E1 if a in kw and b in kw]
            e2=[(a,b,sg) for a,b,sg in E2 if a in kw and b in kw]
            if block_ok(kw,e1,mx1,fx1) and block_ok(kw,e2,mx2,fx2):
                s_.hit=dict(kw); print(f"status: FOUND ansatz {idx} shard {shard} kw={kw}",flush=True); s_.StopSearch()
    solver=cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions=True
    solver.parameters.num_workers=1
    cb=CB()
    tag=f"ansatz {idx}"+(f" shard {shard}/{S}" if shard is not None else "")+(f" size {size}" if size else "")
    print(f"{tag}: block1 {p1}->({a1},{b1})x2={x1}  block2 {p2}->({a2},{b2})x2={x2}",flush=True)
    st=solver.Solve(m,cb)
    if cb.hit:
        return True
    print(f"status: EXHAUSTED-CPSAT  ansatz: {idx}  shard: {shard}  solutions checked: {cb.n}  solver: {solver.StatusName(st)}",flush=True)
    if size is None: _mark_done(key)
    return False

if __name__!="__main__":
    pass
elif len(sys.argv)>1 and sys.argv[1]=="-1":
    print(f"total double-block ansatze: {len(ANS)}")
    for i,(p1,p2,a1,b1,x1,a2,b2,x2,d) in enumerate(ANS):
        print(i,p1,(a1,b1),"x1=",x1," ",p2,(a2,b2),"x2=",x2)
    sys.exit()
elif len(sys.argv)>1 and sys.argv[1]=="--worker":
    W=int(sys.argv[2]); K=int(sys.argv[3]); S=int(sys.argv[4]) if len(sys.argv)>4 else 1
    tasks=[(i,s) for i in range(len(ANS)) for s in (range(S) if S>1 else [None])]
    for n in range(W,len(tasks),K):
        i,s=tasks[n]
        if _run(i,shard=s,S=S): break
    sys.exit()
elif len(sys.argv)>1:
    idx=int(sys.argv[1]); size=int(sys.argv[2]) if len(sys.argv)>2 else None
    _run(idx,size)

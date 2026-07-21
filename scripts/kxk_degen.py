#!/usr/bin/env python3
"""Degenerate 3x3 block solver: the class-count-law keystone.
Ansatz: one 3x3 block on modes (m1,m2,m3) carrying an eigenvalue sub-multiset
E=(e1,e2,e3) with >=2 distinct values (repeats ALLOWED -- the family 2x2
block_ansatze and cliques both miss); integer diagonal (a,b,c) majorized by E;
support with one-hop pairs ONLY on the three block pairs; off-diagonals
s_ij = signed surd sums; phase theta from the exact 3x3 char-poly identity:
  sum s_ij^2 = e2(a,b,c) - e2(E)          (rational, exact)
  cos(theta) = (e3(E) - abc + a s23^2 + b s13^2 + c s12^2) / (2 s12 s13 s23)
Any hit is certified by the full exact char-poly identity (sympy).
Usage: python kxk_degen.py --ints 6,6,1,1,1,1,1,1,1,1,1 --den 7 -N 3 [--budget 200]
"""
import time
import argparse
from itertools import combinations

def sqfree(n):
    m,q,i=1,n,2
    while i*i<=q:
        while q%(i*i)==0: q//=i*i; m*=i
        i+=1
    return m,q
def hop(A,B):
    sA,sB=set(A),set(B); c=sA&sB
    if len(c)!=len(A)-1: return None,None
    (i,)=sA-c; (j,)=sB-c
    return tuple(sorted((i,j))),(-1)**A.index(i)*(-1)**B.index(j)

def surd_add(d_,q,c):
    d_[q]=d_.get(q,0)+c
def surd_sq(a):
    out={}
    for q1,c1 in a.items():
        for q2,c2 in a.items():
            m,q=sqfree(q1*q2)
            surd_add(out,q,c1*c2*m)
    return out
def surd_mul(a,b):
    out={}
    for q1,c1 in a.items():
        for q2,c2 in b.items():
            m,q=sqfree(q1*q2)
            surd_add(out,q,c1*c2*m)
    return {q:c for q,c in out.items() if c}
def rat(a):
    r=a.get(1,0)
    return r, all(q==1 or c==0 for q,c in a.items())

def run(ints, den, N, budget):
    d=len(ints); D=den
    from collections import Counter
    cnt=Counter(ints)
    alldets=list(combinations(range(d),N))
    HOP={}
    for i in range(len(alldets)):
        for j in range(i+1,len(alldets)):
            pr,sg=hop(alldets[i],alldets[j])
            if pr: HOP[(i,j)]=(pr,sg)
    # block eigen-sub-multisets, >=2 distinct, size 3, drawn from spectrum multiset
    Es=set()
    for E in combinations(sorted(ints,reverse=True),3):
        if len(set(E))>=2: Es.add(tuple(sorted(E,reverse=True)))
    configs=[]
    for E in sorted(Es,reverse=True):
        rem=list(ints)
        for e in E: rem.remove(e)
        tot=sum(E)
        # diagonals majorized by E, each in [1,den]
        for a in range(1,min(D,tot-2)+1):
            for b in range(1,min(D,tot-a-1)+1):
                c=tot-a-b
                if c<1 or c>D: continue
                sd=sorted((a,b,c),reverse=True)
                if sd[0]>E[0] or sd[0]+sd[1]>E[0]+E[1]: continue
                if (a,b,c)!=tuple(sd) and (a,b,c)>tuple(sd): continue  # canonical-ish
                # block modes = last 3 modes; others get rem sorted desc
                nv=tuple(list(sorted(rem,reverse=True))+[a,b,c])
                configs.append((E,(a,b,c),nv))
    print(f'{len(configs)} block configs', flush=True)
    bm=(d-3,d-2,d-1)
    bps=[tuple(sorted(p)) for p in combinations(bm,2)]
    t0=time.time(); per=budget/max(1,len(configs)); stats=[]
    for E,(a,b,c),nv in configs:
        deadline=time.time()+per
        conf=[set() for _ in alldets]
        for (i,j),(pr,sg) in HOP.items():
            if pr not in bps:
                conf[i].add(j); conf[j].add(i)
        rows=[[j for j,T in enumerate(alldets) if m in T] for m in range(d)]
        e2E=E[0]*E[1]+E[0]*E[2]+E[1]*E[2]; e3E=E[0]*E[1]*E[2]
        R=(a*b+a*c+b*c)-e2E
        if R<0: continue
        found=[None]; timed=[False]
        def rec(rem,left,chosen,banned):
            if found[0]: return
            if time.time()>deadline: timed[0]=True; return
            if left==0:
                if any(rem): return
                keys=list(chosen)
                terms={bp:[] for bp in bps}
                for x in range(len(keys)):
                    for y in range(x+1,len(keys)):
                        i,j=sorted((keys[x],keys[y]))
                        if (i,j) in HOP:
                            pr,sg=HOP[(i,j)]
                            m0,q0=sqfree(chosen[i]*chosen[j])
                            terms[pr].append((x,y,sg*m0,q0))
                n=len(keys)
                if n>22: return
                core=[x for x in range(n) if any(any(t[0]==x or t[1]==x for t in terms[bp]) for bp in bps)]
                ncore=len(core)
                if ncore>16: return
                for bits in range(1<<max(0,ncore-1)):
                    e=[1]*n
                    for i2,x in enumerate(core[1:]):
                        if bits>>i2&1: e[x]=-1
                    s={}
                    for bp in bps:
                        acc={}
                        for x,y,co,q in terms[bp]: surd_add(acc,q,e[x]*e[y]*co)
                        s[bp]={q:c for q,c in acc.items() if c}
                    tot2={}
                    for bp in bps:
                        for q,cc in surd_sq(s[bp]).items(): surd_add(tot2,q,cc)
                    r_,pure=rat(tot2)
                    if not pure or r_!=R: continue
                    p12,p13,p23=bps
                    lin={}
                    for coef,pp in ((a,p23),(b,p13),(c,p12)):
                        for q,cc in surd_sq(s[pp]).items(): surd_add(lin,q,coef*cc)
                    rl,pl=rat(lin)
                    if not pl: continue
                    prod3=surd_mul(surd_mul(s[p12],s[p13]),s[p23])
                    rhs=e3E-a*b*c+rl
                    num=sum(cc*(q**0.5) for q,cc in prod3.items())
                    if abs(num)<1e-12:
                        if rhs!=0: continue
                        found[0]=(E,(a,b,c),[(alldets[k],chosen[k]) for k in keys],tuple(e),'real-rigid',None)
                        return
                    ct=rhs/(2*num)
                    if abs(ct)<=1+1e-12:
                        found[0]=(E,(a,b,c),[(alldets[k],chosen[k]) for k in keys],tuple(e),'phase',(rhs,dict(prod3)))
                        return
                return
            best=None
            for m in range(d):
                if rem[m]<=0: continue
                avail=[j for j in rows[m] if j not in banned and j not in chosen]
                if not avail: return
                if best is None or (rem[m],len(avail))<(best[1],len(best[2])): best=(m,rem[m],avail)
            m,_,avail=best
            for j in avail:
                T=alldets[j]
                cap=min(left,*(rem[mm] for mm in T))
                for w in range(cap,0,-1):
                    nr=list(rem)
                    for mm in T: nr[mm]-=w
                    ch2=dict(chosen); ch2[j]=w
                    rec(tuple(nr),left-w,ch2,banned|conf[j]|{j})
                    if found[0]: return
        rec(tuple(nv),D,{},set())
        stats.append((E,(a,b,c),'TIMEOUT' if timed[0] else 'EXHAUSTED'))
        if found[0]:
            E_,diag,supp,e,kind,ph=found[0]
            print(f'*** CANDIDATE: block evals {E_}/{den} diag {diag} kind {kind}')
            for (T,w),ee in zip(supp,e): print('   det',T,'w',w,'sign',ee)
            if ph: print('   cos(theta) = ({}) / (2*sqrt-sum {})'.format(ph[0],ph[1]))
            return found[0]
    from collections import Counter as _C
    sc=_C(x[2] for x in stats)
    print(f'no hit: {len(configs)} configs, {dict(sc)}, total {time.time()-t0:.0f}s')
    for x in stats:
        if x[2]=='TIMEOUT': print('  TIMEOUT:', x[0], x[1])
    return None

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--ints',required=True); ap.add_argument('--den',type=int,required=True)
    ap.add_argument('-N',type=int,default=3); ap.add_argument('--budget',type=float,default=200)
    a=ap.parse_args()
    run(tuple(int(x) for x in a.ints.split(',')), a.den, a.N, a.budget)

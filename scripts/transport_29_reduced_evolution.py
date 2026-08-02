"""Exact reduced evolution on the 29-determinant transport carrier.

Evolves a full-support state under the non-diagonal two-body hop the carrier
closes on, and at each sample time rebuilds the carrier matrix, 2-RDM and
3-RDM from certified 2-RDM coordinates alone.
"""
import argparse
import itertools, json, math, time
from collections import defaultdict, deque
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import sys

# Resolve the sibling module from this file's own directory so the script works
# when imported (tests) as well as when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_extensions as te  # noqa: E402

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/transport_29_reduced_evolution.json"


def build_support():
    base=[tuple(sorted(a+b)) for a in te.windows(5,2) for b in te.windows(5,2,5)]
    hop=(0,1,5,7); adj=(1,0,7,5)
    return te.invariant_closure(base,[hop,adj]), hop, adj

def hamiltonian(support, ops):
    idx={d:i for i,d in enumerate(support)}
    H=np.zeros((len(support),len(support)),complex)
    for j,src in enumerate(support):
        for op in ops:
            dst,sg=te.apply_cross(src,*op)
            if dst is not None:
                if dst not in idx: raise RuntimeError('leak')
                H[idx[dst],j]+=sg
    return H

def channels(support,order):
    particles=len(support[0]); terms=defaultdict(lambda:defaultdict(int))
    for a,left in enumerate(support):
        ls=set(left)
        for b,right in enumerate(support):
            inter=sorted(ls.intersection(right)); remn=particles-order
            if len(inter)<remn:continue
            for rem in itertools.combinations(inter,remn):
                rs=set(rem); lm=tuple(x for x in left if x not in rs); rm=tuple(x for x in right if x not in rs)
                terms[(lm,rm)][(a,b)] += te.permutation_sign(lm+rem)*te.permutation_sign(rm+rem)
    return {k:{e:v for e,v in c.items() if v} for k,c in terms.items() if any(c.values())}

def rdm_from_X(support,X,order,orbitals=10):
    feats=list(itertools.combinations(range(orbitals),order)); fi={f:i for i,f in enumerate(feats)}
    G=np.zeros((len(feats),len(feats)),complex)
    for (lm,rm),co in channels(support,order).items():
        G[fi[lm],fi[rm]]=sum(v*X[a,b] for (a,b),v in co.items())
    return G,feats

def pair_incidence(support,orbitals=10):
    feats=list(itertools.combinations(range(orbitals),2))
    B=np.array([[int(set(f).issubset(d)) for d in support] for f in feats],float)
    # choose independent rows via greedy rank
    rows=[]; r=0
    for i in range(len(feats)):
        rr=np.linalg.matrix_rank(B[rows+[i],:],tol=1e-10)
        if rr>r: rows.append(i);r=rr
        if r==len(support):break
    return B,feats,rows

def tree_witnesses(support):
    ch=channels(support,2); adj=defaultdict(list)
    for (lm,rm),co in ch.items():
        if lm>=rm or len(co)!=1:continue
        (a,b),sg=next(iter(co.items()))
        if a==b:continue
        adj[a].append((b,lm,rm,sg));adj[b].append((a,rm,lm,sg))
    tree=[]; seen={0};q=deque([0])
    while q:
        a=q.popleft()
        for b,lm,rm,sg in adj[a]:
            if b in seen:continue
            seen.add(b);q.append(b);tree.append((a,b,lm,rm,sg))
    assert len(seen)==len(support)
    return tree

def reconstruct_X(G2,feats,support,B,rows,tree):
    fi={f:i for i,f in enumerate(feats)}
    diag=np.real(np.diag(G2))
    w=np.linalg.solve(B[rows,:],diag[rows])
    X=np.zeros((len(support),len(support)),complex); np.fill_diagonal(X,w)
    tadj=defaultdict(list)
    for a,b,lm,rm,sg in tree:
        val=G2[fi[lm],fi[rm]]/sg
        X[a,b]=val;X[b,a]=np.conj(val)
        tadj[a].append(b);tadj[b].append(a)
    # root products: recover cbar0*cj along paths using X_uv X_vw / w_v
    seen={0};q=deque([0])
    while q:
        a=q.popleft()
        for b in tadj[a]:
            if b in seen:continue
            seen.add(b);q.append(b)
            if a!=0:
                X[0,b]=X[0,a]*X[a,b]/w[a]
                X[b,0]=np.conj(X[0,b])
    for a in range(len(support)):
        for b in range(len(support)):
            if a==b:continue
            X[a,b]=np.conj(X[0,a])*X[0,b]/w[0] if a!=0 and b!=0 else X[a,b]
    return X

def main():
    support,hop,adj=build_support();H=hamiltonian(support,[hop,adj]);m=len(support)
    # deterministic full-support initial state with nonuniform phases
    k=np.arange(m); c0=(1.0+0.15*np.cos(k))/(np.sqrt(np.sum((1.0+0.15*np.cos(k))**2))) * np.exp(1j*0.17*k)
    B,features,rows=pair_incidence(support); tree=tree_witnesses(support)
    ts=np.linspace(0,4,17); max_x=max_g2=max_g3=0.0; minw=1
    start=time.perf_counter()
    samples=[]
    for t in ts:
        c=expm(-1j*H*t)@c0; X=np.outer(np.conj(c),c)
        G2,f2=rdm_from_X(support,X,2); Xr=reconstruct_X(G2,f2,support,B,rows,tree)
        G2r,_=rdm_from_X(support,Xr,2); G3,_=rdm_from_X(support,X,3); G3r,_=rdm_from_X(support,Xr,3)
        ex=np.max(np.abs(X-Xr));e2=np.max(np.abs(G2-G2r));e3=np.max(np.abs(G3-G3r))
        max_x=max(max_x,ex);max_g2=max(max_g2,e2);max_g3=max(max_g3,e3);minw=min(minw,np.min(np.abs(c)**2))
        samples.append({'t':float(t),'max_X_error':float(ex),'max_Gamma2_error':float(e2),'max_Gamma3_error':float(e3),'min_weight':float(np.min(np.abs(c)**2))})
    elapsed=time.perf_counter()-start
    out={'schema_version':1,'carrier_dimension':m,'ambient_fci_dimension':math.comb(10,4),'pair_basis_dimension':math.comb(10,2),'three_body_basis_dimension':math.comb(10,3),'pair_incidence_rank':int(np.linalg.matrix_rank(B)),'tree_edges':len(tree),'times':len(ts),'max_X_error':max_x,'max_Gamma2_error':max_g2,'max_Gamma3_error':max_g3,'minimum_weight_over_trajectory':minw,'elapsed_seconds':elapsed,'coordinate_counts':{'carrier_wavefunction_complex':m,'full_Gamma2_complex_entries':math.comb(10,2)**2,'certificate_chart_real':m+2*(m-1),'physical_pure_state_real_dimension':2*m-2,'ambient_fci_wavefunction_complex':math.comb(10,4)},'samples':samples}
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=DEFAULT_OUT)
    destination=parser.parse_args().out
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()

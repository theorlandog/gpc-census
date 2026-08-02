#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, statistics, sys, time
from collections import defaultdict, deque
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

# Resolve the sibling modules from this file's own directory so the script
# works when imported (tests) as well as when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multichart_symplectic_propagator as mp  # noqa: E402
import rdm_method_benchmark as base  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results/data/compiled_sparse_benchmark.json"
PREVIOUS_BENCHMARK = ROOT / "results/data/rdm_method_benchmark.json"

class CompiledSparseRHS:
    def __init__(self, atlas: mp.Atlas):
        self.a=atlas; self.m=len(atlas.support); self.e=len(atlas.edge_list)
        self.neigh=[]
        for i in range(self.m):
            self.neigh.append([(int(k),complex(atlas.H[i,k])) for k in np.flatnonzero(np.abs(atlas.H[i])>0)])
        self.calls=0; self.sync_fallbacks=0

    def _root_data(self,y):
        weights,z=self.a.unpack(y); weights=weights.copy(); weights[np.abs(weights)<5e-15]=0
        ci=self.a.choose_chart(weights)
        if ci is None: raise FloatingPointError
        tree=self.a.trees[ci]; root=tree['root']
        ev={edge:value for edge,value in zip(self.a.edge_list,z)}
        children=defaultdict(list)
        for v,p in tree['parent'].items():
            if p>=0: children[p].append(v)
        r=np.zeros(self.m,dtype=complex); r[root]=weights[root]
        q=deque([root])
        while q:
            p=q.popleft()
            for c in children[p]:
                q.append(c); edge=tuple(sorted((p,c))); val=ev[edge]
                oriented=val if edge==(p,c) else np.conj(val)
                if p==root: r[c]=oriented
                elif weights[p]>2e-13: r[c]=r[p]*oriented/weights[p]
                elif weights[c]<=2e-13: r[c]=0
                else: raise FloatingPointError
        if weights[root]<=2e-13: raise FloatingPointError
        return weights,r,root

    @staticmethod
    def _x(weights,r,root,i,j):
        if i==j: return complex(weights[i])
        if weights[i]<=2e-13 or weights[j]<=2e-13: return 0j
        if i==root: return r[j]
        if j==root: return np.conj(r[i])
        return np.conj(r[i])*r[j]/weights[root]

    def __call__(self,_t,y):
        self.calls+=1
        try:
            w,r,root=self._root_data(y)
            x=lambda i,j:self._x(w,r,root,i,j)
            dw=np.empty(self.m)
            for a in range(self.m):
                s1=sum(h*x(k,a) for k,h in self.neigh[a])
                s2=sum(x(a,k)*np.conj(h) for k,h in self.neigh[a])
                dw[a]=float(np.real(1j*(s1-s2)))
            dz=np.empty(self.e,dtype=complex)
            for n,(a,b) in enumerate(self.a.edge_list):
                s1=sum(h*x(k,b) for k,h in self.neigh[a])
                s2=sum(x(a,k)*np.conj(h) for k,h in self.neigh[b])
                dz[n]=1j*(s1-s2)
            return np.concatenate([dw,dz.real,dz.imag])
        except FloatingPointError:
            self.sync_fallbacks+=1
            return self.a.rhs(_t,y)

def run_one(method):
    a=mp.Atlas.build(); c0=base.initial_state(a); X0=np.outer(np.conj(c0),c0); y0=a.pack_from_X(X0)
    Xref=np.outer(np.conj(expm(-1j*a.H*base.T_FINAL)@c0),expm(-1j*a.H*base.T_FINAL)@c0)
    rhs=a.rhs if method=='sparse2rdm' else CompiledSparseRHS(a)
    t=time.perf_counter(); sol=solve_ivp(rhs,(0,base.T_FINAL),y0,rtol=base.RTOL,atol=base.ATOL,method='DOP853',max_step=base.MAX_STEP); elapsed=time.perf_counter()-t
    try: X,_=a.reconstruct_X(sol.y[:,-1])
    except FloatingPointError: X,_=a.reconstruct_by_synchronization(sol.y[:,-1])
    return {'method':method,'elapsed_seconds':elapsed,'nfev':sol.nfev,'steps':len(sol.t)-1,'error':float(np.max(np.abs(X-Xref))), 'fallbacks':getattr(rhs,'sync_fallbacks',None), 'variables':len(y0)}

def orchestrate(out,repeats):
    raw={m:[run_one(m) for _ in range(repeats)] for m in ['sparse2rdm','compiled_sparse2rdm']}
    summary=[]
    for m,rows in raw.items():
        summary.append({'method':m,'median_seconds':statistics.median(x['elapsed_seconds'] for x in rows),'min_seconds':min(x['elapsed_seconds'] for x in rows),'max_error':max(x['error'] for x in rows),'nfev':rows[0]['nfev'],'steps':rows[0]['steps'],'variables':rows[0]['variables'],'fallbacks':rows[0]['fallbacks']})
    by={x['method']:x for x in summary}
    # include previous benchmark medians for context
    prev=json.loads(PREVIOUS_BENCHMARK.read_text())
    old={x['method']:x for x in prev['summary']}
    result={'schema_version':1,'summary':summary,'speedup_compiled_over_original':by['sparse2rdm']['median_seconds']/by['compiled_sparse2rdm']['median_seconds'],'context':{'carrier_wavefunction_seconds':old['carrier']['median_elapsed_seconds'],'full_fci_seconds':old['fci']['median_elapsed_seconds'],'dense2rdm_seconds':old['dense2rdm']['median_elapsed_seconds']},'ratios':{'compiled_vs_carrier':by['compiled_sparse2rdm']['median_seconds']/old['carrier']['median_elapsed_seconds'],'compiled_vs_fci':by['compiled_sparse2rdm']['median_seconds']/old['fci']['median_elapsed_seconds'],'dense_vs_compiled':old['dense2rdm']['median_elapsed_seconds']/by['compiled_sparse2rdm']['median_seconds']},'raw':raw}
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n'); return result

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=DEFAULT_OUT); ap.add_argument('--repeats',type=int,default=7); args=ap.parse_args(); print(json.dumps(orchestrate(args.out,args.repeats),indent=2));
if __name__=='__main__': main()

import sys
import json
import itertools
import numpy as np
from scipy.optimize import minimize

from checkpoint import Checkpointer, extract_opts

def build(r,N):
    dets=list(itertools.combinations(range(r),N)); idx={t:i for i,t in enumerate(dets)}
    D=len(dets); A=np.zeros((r,r,D,D),complex)
    for t in dets:
        for mp in t:
            s1=(-1)**t.index(mp); t2=tuple(x for x in t if x!=mp)
            for m in range(r):
                if m in t2: continue
                tp=tuple(sorted(t2+(m,))); s2=(-1)**tp.index(m)
                A[m,mp,idx[tp],idx[t]]=s1*s2
    return dets,A

def solve(r,N,lam,A=None,dets=None,mask=None,outer=250,tries=5,tol=1e-16):
    if A is None: dets,A=build(r,N)
    D=A.shape[2]; lam=np.sort(np.array(lam,float))[::-1]
    best=(None,1e9)
    for _ in range(tries):
        psi=np.random.randn(D)+1j*np.random.randn(D)
        if mask is not None: psi*=mask
        psi/=np.linalg.norm(psi)
        for it in range(outer):
            rho=np.einsum('i,mnij,j->mn',psi.conj(),A,psi)
            e,U=np.linalg.eigh(rho)
            T=U[:,::-1]@np.diag(lam)@U[:,::-1].conj().T
            def fg(x):
                p=x[:D]+1j*x[D:]
                if mask is not None: p=p*mask
                n=np.linalg.norm(p); p=p/max(n,1e-14)
                rho=np.einsum('i,mnij,j->mn',p.conj(),A,p)
                M=rho-T
                E=float(np.real(np.sum(M*M.conj())))
                B=np.einsum('nm,mnij->ij',M,A)
                g=B@p; g=g-np.vdot(p,g).real*p
                if mask is not None: g=g*mask
                return E,np.concatenate([g.real,g.imag])*2/max(n,1e-14)
            x=np.concatenate([psi.real,psi.imag])
            res=minimize(fg,x,jac=True,method='L-BFGS-B',options={'maxiter':60,'ftol':1e-18,'gtol':1e-14})
            psi=res.x[:D]+1j*res.x[D:]
            if mask is not None: psi=psi*mask
            psi/=np.linalg.norm(psi)
            rho=np.einsum('i,mnij,j->mn',psi.conj(),A,psi)
            e=np.linalg.eigvalsh(rho)[::-1].real
            spec_res=float(np.sum((e-lam)**2))
            if spec_res<tol: return psi,spec_res
        if spec_res<best[1]: best=(psi,spec_res)
    return best

def minimize_support(r,N,lam,psi,dets,A,res_tol=1e-12):
    D=len(psi); mask=(np.abs(psi)>1e-10).astype(float)
    improved=True
    while improved:
        improved=False
        for i in sorted(np.where(mask>0)[0],key=lambda i:abs(psi[i])):
            m2=mask.copy(); m2[i]=0
            if m2.sum()<1: continue
            out=solve(r,N,lam,A,dets,mask=m2,outer=120,tries=3)
            if out[1]<res_tol:
                psi,mask=out[0],m2; improved=True; break
    return psi,mask

def run_one(r,N,den,n,A=None,dets=None):
    lam=[x/den for x in n]
    if A is None: dets,A=build(r,N)
    psi,res=solve(r,N,lam,A,dets)
    if res>1e-12:
        rec={"r":r,"N":N,"den":den,"n":n,"status":"FAIL","residual":res}
        print(json.dumps(rec)); return rec
    psi,mask=minimize_support(r,N,lam,psi,dets,A)
    sup=[i for i in range(len(psi)) if mask[i] and abs(psi[i])>1e-9]
    j0=max(sup,key=lambda i:abs(psi[i])); psi=psi*np.exp(-1j*np.angle(psi[j0]))
    rec={"r":r,"N":N,"den":den,"n":n,"status":"OK","residual":res,
        "support_size":len(sup),
        "support":[[list(dets[i]),round(float(abs(psi[i])),12),round(float(np.angle(psi[i])),12)] for i in sup]}
    print(json.dumps(rec)); return rec

if __name__=="__main__":
    _opts,argv=extract_opts(sys.argv)
    if argv[1]=="one":
        run_one(int(argv[2]),int(argv[3]),int(argv[4]),[int(x) for x in argv[5].split(",")])
    else:
        tasks=json.load(open(argv[2])); cache={}
        # With --out each state is persisted (one JSON record per line) and the
        # run resumes: already-recorded task indices are skipped and the file is
        # checkpointed (fsync + optional S3 mirror) to survive spot interruptions.
        ckpt=Checkpointer.from_opts(None,_opts)
        if ckpt: ckpt.restore()
        out_path=_opts.get("out"); done=set()
        if out_path:
            try:
                for line in open(out_path): done.add(json.loads(line)["index"])
            except (FileNotFoundError,json.JSONDecodeError): pass
        out=open(out_path,"a") if out_path else None
        if ckpt and out: ckpt.attach(out); ckpt.install_signal_handlers()
        for i,(r,N,den,n) in enumerate(tasks):
            if i in done: continue
            if (r,N) not in cache: cache[(r,N)]=build(r,N)
            dets,A=cache[(r,N)]
            rec=run_one(r,N,den,n,A,dets)
            if out:
                rec=dict(rec); rec["index"]=i
                out.write(json.dumps(rec)+"\n"); out.flush(); ckpt.checkpoint()
        if ckpt and out: ckpt.checkpoint(force=True)
        if out: out.close()

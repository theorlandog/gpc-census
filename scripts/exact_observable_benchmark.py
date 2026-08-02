#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, statistics, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm

# Resolve the sibling modules from this file's own directory so the script
# works when imported (tests) as well as when run directly.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multichart_symplectic_propagator as mp  # noqa: E402
import rdm_method_benchmark as base  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "results/data/exact_observable_benchmark.json"
DEFAULT_CONTROL_OUT = ROOT / "results/data/exact_observable_benchmark_with_control.json"
PREVIOUS_BENCHMARK = ROOT / "results/data/rdm_method_benchmark.json"
COMPILED_BENCHMARK = ROOT / "results/data/compiled_sparse_benchmark.json"

class ExactObservablePropagator:
    """Exact observable-coordinate propagator for the four-block transport H.

    It consumes and returns the redundant certified 2-RDM chart. Internally it
    reconstructs the carrier rank-one matrix from observables, applies the
    exact sparse unitary conjugation, and extracts only certified observables.
    No wavefunction amplitudes are read or propagated.
    """
    def __init__(self, atlas: mp.Atlas):
        self.a = atlas
        self.blocks=[]
        used=set()
        for i in range(len(atlas.H)):
            js=np.flatnonzero(np.abs(atlas.H[i])>1e-14)
            if len(js)==1:
                j=int(js[0])
                if i<j:
                    self.blocks.append((i,j,float(np.real(atlas.H[i,j]))))
                    used|={i,j}
        if len(self.blocks)!=4:
            raise RuntimeError(f'expected four 2x2 hopping blocks, got {self.blocks}')

    def propagate(self, y0, t):
        try:
            X,_=self.a.reconstruct_X(y0)
        except FloatingPointError:
            X,_=self.a.reconstruct_by_synchronization(y0)
        c=np.cos(t); s=np.sin(t)
        # X' = exp(iHt) X exp(-iHt). Apply four disjoint 2x2 rotations in-place.
        Y=X.copy()
        for i,j,sign in self.blocks:
            # Left multiplication on rows i,j by exp(+i sign sigma_x t)
            ri=Y[i,:].copy(); rj=Y[j,:].copy()
            Y[i,:]=c*ri + 1j*sign*s*rj
            Y[j,:]=1j*sign*s*ri + c*rj
        for i,j,sign in self.blocks:
            # Right multiplication on columns i,j by exp(-i sign sigma_x t)
            ci=Y[:,i].copy(); cj=Y[:,j].copy()
            Y[:,i]=c*ci - 1j*sign*s*cj
            Y[:,j]=-1j*sign*s*ci + c*cj
        return self.a.pack_from_X(Y)


def build_case():
    a=mp.Atlas.build(); c0=base.initial_state(a); X0=np.outer(np.conj(c0),c0); y0=a.pack_from_X(X0)
    return a,c0,X0,y0


def verify(t=base.T_FINAL):
    a,c0,X0,y0=build_case(); p=ExactObservablePropagator(a); y=p.propagate(y0,t)
    X,_=a.reconstruct_X(y)
    U=expm(-1j*a.H*t); c=U@c0; ref=np.outer(np.conj(c),c)
    return float(np.max(np.abs(X-ref)))


def time_exact(repeats=2000):
    a,c0,X0,y0=build_case(); p=ExactObservablePropagator(a)
    # warmup
    for _ in range(20): p.propagate(y0,base.T_FINAL)
    samples=[]
    for _ in range(repeats):
        t0=time.perf_counter_ns(); y=p.propagate(y0,base.T_FINAL); samples.append((time.perf_counter_ns()-t0)*1e-9)
    X,_=a.reconstruct_X(y)
    U=expm(-1j*a.H*base.T_FINAL); c=U@c0; ref=np.outer(np.conj(c),c)
    return {'method':'exact_observable','median_elapsed_seconds':statistics.median(samples),'minimum_elapsed_seconds':min(samples),'p95_elapsed_seconds':float(np.percentile(samples,95)),'real_variables':len(y0),'max_reference_error':float(np.max(np.abs(X-ref))),'repeats':repeats,'blocks':[list(x) for x in p.blocks]}


def run(out:Path,repeats=2000):
    exact=time_exact(repeats)
    prev=json.loads(PREVIOUS_BENCHMARK.read_text())
    comp=json.loads(COMPILED_BENCHMARK.read_text())
    baseline={r['method']:r for r in prev['summary']}
    compiled={r['method']:r for r in comp['summary']}['compiled_sparse2rdm']
    comparisons={
      'vs_carrier':baseline['carrier']['median_elapsed_seconds']/exact['median_elapsed_seconds'],
      'vs_fci':baseline['fci']['median_elapsed_seconds']/exact['median_elapsed_seconds'],
      'vs_dense2rdm':baseline['dense2rdm']['median_elapsed_seconds']/exact['median_elapsed_seconds'],
      'vs_original_sparse2rdm':baseline['sparse2rdm']['median_elapsed_seconds']/exact['median_elapsed_seconds'],
      'vs_compiled_sparse2rdm':compiled['median_seconds']/exact['median_elapsed_seconds'],
    }
    result={'schema_version':1,'method':exact,'baseline_medians':{k:v['median_elapsed_seconds'] for k,v in baseline.items()},'compiled_sparse_median':compiled['median_seconds'],'speedups':comparisons,'scope':'Specialized exact propagator for the fixed four-block hopping Hamiltonian; not a generic arbitrary-H integrator.'}
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


def specialized_carrier_update(blocks, c, t):
    """The control: the same four analytic rotations on carrier amplitudes.

    This is what an implementer who already has the wavefunction would write
    for this Hamiltonian. It is the fair opponent for the observable
    propagator, and it is deliberately given the same special-case knowledge.
    """
    out=c.copy(); cs=np.cos(t); sn=np.sin(t)
    for i,j,sign in blocks:
        ci=out[i]; cj=out[j]
        out[i]=cs*ci - 1j*sign*sn*cj
        out[j]=-1j*sign*sn*ci + cs*cj
    return out


def time_specialized_carrier(repeats=20000):
    a,c0,_X0,_y0=build_case(); blocks=ExactObservablePropagator(a).blocks
    reference=expm(-1j*a.H*base.T_FINAL)@c0
    error=float(np.max(np.abs(specialized_carrier_update(blocks,c0,base.T_FINAL)-reference)))
    for _ in range(200): specialized_carrier_update(blocks,c0,base.T_FINAL)
    samples=[]
    for _ in range(repeats):
        t0=time.perf_counter_ns(); specialized_carrier_update(blocks,c0,base.T_FINAL); samples.append((time.perf_counter_ns()-t0)*1e-9)
    return {'median_elapsed_seconds':statistics.median(samples),'minimum_elapsed_seconds':min(samples),'max_reference_error':error,'repeats':repeats}


def run_with_control(out:Path,exact:dict,control_repeats=20000):
    """Emit the control comparison as a generated artifact, not by hand."""
    control=time_specialized_carrier(control_repeats)
    result={
      'schema_version':1,
      'exact_observable_median_seconds':exact['method']['median_elapsed_seconds'],
      'exact_observable_error':exact['method']['max_reference_error'],
      'previous_ode_speedups':{
        'carrier':exact['speedups']['vs_carrier'],
        'full_fci':exact['speedups']['vs_fci'],
        'dense_2rdm':exact['speedups']['vs_dense2rdm'],
        'original_sparse_2rdm':exact['speedups']['vs_original_sparse2rdm'],
        'compiled_sparse_2rdm':exact['speedups']['vs_compiled_sparse2rdm'],
      },
      'specialized_exact_carrier_median_seconds':control['median_elapsed_seconds'],
      'specialized_exact_carrier_error':control['max_reference_error'],
      'specialized_carrier_speedup_over_observable':exact['method']['median_elapsed_seconds']/control['median_elapsed_seconds'],
      'work_accounting':'The observable propagator is evaluated once per call; the ODE baselines it is compared against take 314 right-hand-side evaluations across 26 adaptive steps. The comparison is an analytic solution against numerical integration of the same trajectory, not a per-evaluation comparison.',
      'conclusion':'The exact observable propagator beats every prior generic ODE implementation but not an equally specialized exact carrier propagator.',
    }
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=DEFAULT_OUT)
    ap.add_argument('--control-out',type=Path,default=DEFAULT_CONTROL_OUT)
    ap.add_argument('--repeats',type=int,default=2000)
    ap.add_argument('--control-repeats',type=int,default=20000)
    args=ap.parse_args()
    exact=run(args.out,args.repeats)
    control=run_with_control(args.control_out,exact,args.control_repeats)
    print(json.dumps({'exact':exact,'with_control':control},indent=2))

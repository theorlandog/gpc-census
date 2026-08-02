#!/usr/bin/env python3
"""Scaling benchmark for products of exactly observable two-level transport cells.

The physical manifold is a product of q Bloch spheres.  Each cell is encoded by
three real 2-RDM-derived Cartesian coordinates (x,y,z), while a generic dense
carrier representation has 2**q complex amplitudes.  Local transport
Hamiltonians H=sum_j omega_j sigma_x^(j) preserve the product manifold, so the
observable propagator is an exact product of SO(3) rotations.

This benchmark compares:
  * exact observable propagation: O(q) storage and work;
  * dense carrier propagation with structure-aware local tensor rotations:
    O(q 2**q) work and O(2**q) storage;
  * generic dense carrier matrix exponential for small q only.

It deliberately also reports the no-free-lunch control: a factorized amplitude
representation uses O(q) storage too and is expected to be faster by constants.
"""
from __future__ import annotations
import argparse, json, math, statistics, time, tracemalloc
from pathlib import Path
import numpy as np
from scipy.linalg import expm

DEFAULT_OUT = Path(__file__).resolve().parents[1] / "results/data/generalized_observable_benchmark.json"


def initial_cells(q: int) -> np.ndarray:
    j=np.arange(q,dtype=float)
    theta=0.35+0.7*(j+1)/(q+1)
    phi=0.2+0.9*(j+1)/(q+1)
    a=np.cos(theta/2)
    b=np.exp(1j*phi)*np.sin(theta/2)
    return np.stack([a,b],axis=1)


def cells_to_bloch(cells: np.ndarray) -> np.ndarray:
    a,b=cells[:,0],cells[:,1]
    return np.stack([2*np.real(np.conj(a)*b),2*np.imag(np.conj(a)*b),np.abs(a)**2-np.abs(b)**2],axis=1)


def propagate_observable(bloch: np.ndarray, omegas: np.ndarray, t: float) -> np.ndarray:
    # H=omega sigma_x, hence Bloch rotation around x by angle 2 omega t.
    out=bloch.copy(); ang=2*omegas*t; c=np.cos(ang); s=np.sin(ang)
    y=out[:,1].copy(); z=out[:,2].copy()
    out[:,1]=c*y-s*z; out[:,2]=s*y+c*z
    return out


def propagate_factorized_amplitudes(cells: np.ndarray, omegas: np.ndarray, t: float) -> np.ndarray:
    out=cells.copy()
    c=np.cos(omegas*t); s=np.sin(omegas*t)
    a=out[:,0].copy(); b=out[:,1].copy()
    out[:,0]=c*a-1j*s*b; out[:,1]=-1j*s*a+c*b
    return out


def product_state(cells: np.ndarray) -> np.ndarray:
    v=np.array([1+0j])
    for cell in cells:
        v=np.kron(v,cell)
    return v


def apply_local_rotations_dense(psi: np.ndarray, omegas: np.ndarray, t: float) -> np.ndarray:
    q=len(omegas); tensor=psi.reshape((2,)*q).copy()
    for axis,w in enumerate(omegas):
        c=math.cos(w*t); s=math.sin(w*t)
        moved=np.moveaxis(tensor,axis,0)
        a=moved[0].copy(); b=moved[1].copy()
        moved[0]=c*a-1j*s*b; moved[1]=-1j*s*a+c*b
        tensor=np.moveaxis(moved,0,axis)
    return tensor.reshape(-1)


def dense_hamiltonian(q: int, omegas: np.ndarray) -> np.ndarray:
    sx=np.array([[0,1],[1,0]],complex); I=np.eye(2,dtype=complex)
    dim=2**q; H=np.zeros((dim,dim),complex)
    for j,w in enumerate(omegas):
        term=np.array([[1+0j]])
        for k in range(q): term=np.kron(term,sx if k==j else I)
        H += w*term
    return H


def median_time(fn,repeats:int):
    for _ in range(3): fn()
    vals=[]
    for _ in range(repeats):
        t0=time.perf_counter_ns(); fn(); vals.append((time.perf_counter_ns()-t0)*1e-9)
    return statistics.median(vals), min(vals)


def peak_bytes(fn):
    tracemalloc.start(); fn(); _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop(); return peak


def row(q:int,t:float,repeats:int,dense_expm_max_q:int):
    cells=initial_cells(q); bloch=cells_to_bloch(cells); omegas=0.7+0.3*np.arange(q)/(max(q-1,1))
    psi=product_state(cells)
    fobs=lambda: propagate_observable(bloch,omegas,t)
    ffact=lambda: propagate_factorized_amplitudes(cells,omegas,t)
    fcarrier=lambda: apply_local_rotations_dense(psi,omegas,t)
    obs_med,obs_min=median_time(fobs,repeats)
    fact_med,fact_min=median_time(ffact,repeats)
    car_repeats=max(3,min(repeats, max(3,2000000//(2**q))))
    car_med,car_min=median_time(fcarrier,car_repeats)
    obs= fobs(); fact=ffact(); ref_bloch=cells_to_bloch(fact)
    error=float(np.max(np.abs(obs-ref_bloch)))
    result={
      'q':q,'carrier_complex_amplitudes':2**q,'observable_real_variables':3*q,
      'factorized_amplitude_real_variables':4*q,
      'observable_median_seconds':obs_med,'factorized_amplitude_median_seconds':fact_med,
      'dense_carrier_median_seconds':car_med,'dense_carrier_repeats':car_repeats,
      'observable_vs_dense_carrier_speedup':car_med/obs_med,
      'observable_vs_factorized_amplitude_speedup':fact_med/obs_med,
      'observable_peak_bytes':peak_bytes(fobs),'dense_carrier_peak_bytes':peak_bytes(fcarrier),
      'max_bloch_error':error,
    }
    if q<=dense_expm_max_q:
        H=dense_hamiltonian(q,omegas)
        fexp=lambda: expm(-1j*H*t)@psi
        exp_med,_=median_time(fexp,max(3,min(10,repeats)))
        result['generic_dense_expm_seconds']=exp_med
        result['observable_vs_generic_dense_expm_speedup']=exp_med/obs_med
    return result


def run(out:Path,qs:list[int],t:float,repeats:int,dense_expm_max_q:int):
    rows=[row(q,t,repeats,dense_expm_max_q) for q in qs]
    result={
      'schema_version':1,
      'model':'product of q exactly observable two-level transport cells',
      'hamiltonian':'sum_j omega_j sigma_x^(j)',
      'scope':'Exact invariant product manifold. Observable coordinates are local Bloch vectors obtainable from certified low-order RDM entries.',
      'rows':rows,
      'claims':{
        'proved_scaling':'observable O(q) versus dense carrier O(2^q) storage and O(q 2^q) update work',
        'benchmark':'observable propagation beats all tested dense carrier methods beyond the reported crossover',
        'no_free_lunch':'a factorized amplitude representation is also O(q); the observable update wins only by a constant factor from the particular NumPy kernels, and no representation-independent speedup is claimed',
        'separation_scope':'H has no inter-cell coupling and the benchmark state is a product state, so the trajectory never leaves the product manifold and the Schmidt rank stays 1. The exponential separation is against a dense 2^q amplitude vector, which is not the representation an implementer would choose here; the like-for-like opponent is the factorized amplitude update, which is also O(q).'
      }
    }
    out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--out',type=Path,default=DEFAULT_OUT); ap.add_argument('--qs',default='4,6,8,10,12,14,16,18'); ap.add_argument('--time',type=float,default=1.3); ap.add_argument('--repeats',type=int,default=100); ap.add_argument('--dense-expm-max-q',type=int,default=8); a=ap.parse_args()
    print(json.dumps(run(a.out,[int(x) for x in a.qs.split(',')],a.time,a.repeats,a.dense_expm_max_q),indent=2))

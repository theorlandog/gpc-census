#!/usr/bin/env python3
"""Compact independent replay of the rank-7/rank-8 audit certificates."""
from __future__ import annotations
import hashlib
import itertools
import json
from collections import deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "results" / "data" / "full_rank78_compression_audit.json"

def normalized_levels(tau):
    values=sorted(set(tau)); step=0
    import math
    for a,b in zip(values,values[1:]): step=math.gcd(step,b-a)
    step=step or 1
    normalized=[(v-values[0])//step for v in values]
    return normalized, sorted(set(range(normalized[-1]+1))-set(normalized))

def tangent_matrix(tau, n=3):
    d=len(tau); determinants=list(itertools.combinations(range(d),n))
    on=[det for det in determinants if sum(tau[i] for i in det)==0]
    below=[det for det in determinants if sum(tau[i] for i in det)>0]
    roots=[(i,j) for i in range(d) for j in range(i+1,d) if tau[i]<tau[j]]
    row_of={det:r for r,det in enumerate(below)}
    matrix=[[None]*len(roots) for _ in below]
    for variable,det in enumerate(on):
        support=set(det)
        for column,(i,j) in enumerate(roots):
            if i not in support or j in support: continue
            occupied=list(det); annihilation=occupied.index(i); occupied.pop(annihilation)
            creation=sum(value<j for value in occupied)
            sign=-1 if (annihilation+creation)%2 else 1
            target=tuple(sorted(occupied+[j])); row=row_of.get(target)
            if row is not None:
                assert matrix[row][column] is None
                matrix[row][column]=(variable,sign)
    return on,below,roots,matrix

def maximum_matching(matrix):
    adjacency=[[j for j,e in enumerate(row) if e is not None] for row in matrix]
    left=len(adjacency); right=len(matrix[0]) if matrix else 0
    lm=[-1]*left; rm=[-1]*right; dist=[0]*left
    def bfs():
        q=deque(); found=False
        for u in range(left):
            if lm[u]<0: dist[u]=0; q.append(u)
            else: dist[u]=10**9
        while q:
            u=q.popleft()
            for v in adjacency[u]:
                w=rm[v]
                if w<0: found=True
                elif dist[w]==10**9: dist[w]=dist[u]+1; q.append(w)
        return found
    def dfs(u):
        for v in adjacency[u]:
            w=rm[v]
            if w<0 or (dist[w]==dist[u]+1 and dfs(w)):
                lm[u]=v; rm[v]=u; return True
        dist[u]=10**9; return False
    count=0
    while bfs():
        for u in range(left):
            if lm[u]<0 and dfs(u): count+=1
    return count

def sparse_determinant(matrix):
    size=len(matrix); assert all(len(row)==size for row in matrix)
    states={0:{():1}}
    for row_index,row in enumerate(matrix):
        next_states={}
        for used,polynomial in states.items():
            assert used.bit_count()==row_index
            for column,edge in enumerate(row):
                if edge is None or used&(1<<column): continue
                variable,edge_sign=edge
                permutation_sign=-1 if (used>>(column+1)).bit_count()%2 else 1
                target=next_states.setdefault(used|(1<<column),{})
                for monomial,coefficient in polynomial.items():
                    extended=tuple(sorted((*monomial,variable)))
                    updated=target.get(extended,0)+edge_sign*permutation_sign*coefficient
                    if updated: target[extended]=updated
                    else: target.pop(extended,None)
        states={mask:poly for mask,poly in next_states.items() if poly}
        if not states: return {}
    return states.get((1<<size)-1,{})

def main():
    payload=json.loads(ARTIFACT.read_text())
    digest=payload.pop('canonical_sha256_without_this_field')
    assert digest==hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    s7=payload['systems']['(3,7)']; s8=payload['systems']['(3,8)']; summary=payload['summary']
    assert sum(s7['hole_distribution_all'].values())==549
    assert sum(s8['hole_distribution_all'].values())==6607
    assert s7['hall_zero']+s7['hall_feasible']==549
    assert s8['hall_zero']+s8['hall_feasible']==6607
    assert sum(s7['optimal_fusion_width_distribution_all'].values())==549
    assert sum(s8['optimal_fusion_width_distribution_all'].values())==6607
    assert s7['singleton_signature_count']==s7['nonstructural_survivors']==547
    assert s8['singleton_signature_count']==s8['nonstructural_survivors']==6605
    assert summary['total_trace_survivors']==7156
    rows=s8['one_hole_hall_feasible_exact_audit']; assert len(rows)==7
    for record in rows:
        tau=tuple(record['tau']); levels,holes=normalized_levels(tau)
        assert levels==record['normalized_levels']==[0,1,2,4]
        assert holes==record['missing_levels']==[3]
        on,below,roots,matrix=tangent_matrix(tau)
        assert len(on)==record['on_variables']==7
        assert len(below)==len(roots)==record['matrix_size']==12
        assert maximum_matching(matrix)==12
        coefficients=sparse_determinant(matrix)
        assert coefficients=={}
        assert record['determinant_coefficient_count']==0
        assert record['determinant_nonzero'] is False
    print('full rank-7/rank-8 compact certificate verifier: PASS')
    return 0
if __name__=='__main__': raise SystemExit(main())

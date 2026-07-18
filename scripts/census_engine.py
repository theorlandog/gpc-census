"""Rank-9/10 census engine. For each vertex (integer numerator vector + denominator):
  stage 1: integer design at natural denominator  (CP-SAT feasibility)
  stage 2: real-weight design                      (CBC continuous MILP)
  verdict: DESIGN-INT / DESIGN-REAL / INTERFERENCE
Usage:
  python census_engine.py test                 # validation suite vs proven results
  python census_engine.py file vertices.json   # [[r,N,denom,[n...]], ...]
  python census_engine.py one r N denom n1,n2,...
"""
import sys
import json
from itertools import combinations
from ortools.sat.python import cp_model
import pulp

from checkpoint import Checkpointer, extract_opts

def design_int(n, r, N, denom):
    dets=list(combinations(range(r),N)); nd=len(dets)
    m=cp_model.CpModel()
    k=[m.NewIntVar(0,denom,f"k{t}") for t in range(nd)]
    y=[m.NewBoolVar(f"y{t}") for t in range(nd)]
    for t in range(nd):
        m.Add(k[t]<=denom*y[t]); m.Add(k[t]>=y[t])
    for mo in range(r):
        m.Add(sum(k[t] for t in range(nd) if mo in dets[t])==n[mo])
    for a in range(nd):
        Ta=set(dets[a])
        for b in range(a+1,nd):
            if len(Ta&set(dets[b]))==N-1:
                m.AddBoolOr([y[a].Not(),y[b].Not()])
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=120; s.parameters.num_workers=2
    st=s.Solve(m)
    return {cp_model.OPTIMAL:"FEASIBLE",cp_model.FEASIBLE:"FEASIBLE",
            cp_model.INFEASIBLE:"INFEASIBLE"}.get(st,"UNKNOWN")

def design_real(n, r, N):
    dets=list(combinations(range(r),N)); nd=len(dets); D=float(sum(n))/N
    prob=pulp.LpProblem("d",pulp.LpMinimize)
    k=[pulp.LpVariable(f"k{t}",0,D) for t in range(nd)]
    y=[pulp.LpVariable(f"y{t}",cat="Binary") for t in range(nd)]
    for t in range(nd): prob += k[t] <= D*y[t]
    for mo in range(r):
        prob += pulp.lpSum(k[t] for t in range(nd) if mo in dets[t]) == n[mo]
    for a in range(nd):
        Ta=set(dets[a])
        for b in range(a+1,nd):
            if len(Ta&set(dets[b]))==N-1: prob += y[a]+y[b] <= 1
    prob += 0
    st=prob.solve(pulp.PULP_CBC_CMD(msg=0,timeLimit=300))
    return pulp.LpStatus[st]

def classify(r,N,denom,n):
    di=design_int(n,r,N,denom)
    if di=="FEASIBLE": return "DESIGN-INT"
    dr=design_real(n,r,N)
    if dr=="Optimal": return "DESIGN-REAL"
    if dr=="Infeasible" and di=="INFEASIBLE": return "INTERFERENCE"
    return f"UNRESOLVED({di},{dr})"

_opts, argv = extract_opts(sys.argv)

if argv[1]=="test":
    print("=== validation vs proven classifications ===")
    cases=[
        (9,4,21,[16,16,16,6,6,6,6,6,6],"DESIGN-INT"),      # v_A
        (9,4,23,[20,14,14,14,14,4,4,4,4],"INTERFERENCE"),  # v_B
        (8,3,8,[1,1,1,1,1,1,1,0],"DESIGN-REAL"),
        (8,3,17,[3,3,3,3,3,1,1,1],"DESIGN-REAL"),
        (8,3,18,[4,4,4,4,2,1,1,1],"INTERFERENCE"),
        (8,3,26,[5,5,5,5,2,2,2,2],"INTERFERENCE"),
        (8,3,23,[6,5,5,5,2,2,1,1],"INTERFERENCE"),
        (8,3,27,[6,6,4,4,4,1,1,1],"INTERFERENCE"),
        (8,3,30,[7,5,5,5,2,2,2,2],"INTERFERENCE"),
        (8,3,30,[7,7,4,4,4,2,1,1],"INTERFERENCE"),
        (8,3,36,[9,5,5,5,3,3,3,3],"INTERFERENCE"),
        (8,3,45,[9,8,8,8,3,3,3,3],"INTERFERENCE"),
        (8,3,48,[9,9,9,9,4,4,2,2],"INTERFERENCE"),
        (8,3,57,[10,10,10,10,4,4,3,3],"INTERFERENCE"),
        (8,3,54,[12,12,7,7,4,4,4,4],"INTERFERENCE"),
    ]
    allok=True
    for r,N,dn,n,expect in cases:
        got=classify(r,N,dn,n)
        ok = got==expect or (expect=="DESIGN-REAL" and got=="DESIGN-INT")
        allok &= ok
        print(f"  {str(n):34s} -> {got:14s} expect {expect:14s} {'OK' if ok else '*** FAIL ***'}")
    print("ALL VALIDATIONS PASS:" , allok)
elif argv[1]=="file":
    rows=json.load(open(argv[2]))
    # With --out the verdicts are persisted one JSON record per line and the run
    # is resumable: already-recorded indices are skipped and the file is
    # checkpointed (fsync + optional S3 mirror) so a spot kill loses little work.
    ckpt=Checkpointer.from_opts(None,_opts)
    if ckpt: ckpt.restore()
    out_path=_opts.get("out")
    done=set(); tally={}
    if out_path:
        try:
            for line in open(out_path):
                rec=json.loads(line); done.add(rec["index"]); tally[rec["verdict"]]=tally.get(rec["verdict"],0)+1
        except (FileNotFoundError,json.JSONDecodeError): pass
    out=open(out_path,"a") if out_path else None
    if ckpt and out:
        ckpt.attach(out); ckpt.install_signal_handlers()
    for i,(r,N,dn,n) in enumerate(rows):
        if i in done: continue
        v=classify(r,N,dn,n)
        tally[v]=tally.get(v,0)+1
        print(f"{i:4d} {str(n):40s} {v}",flush=True)
        if out:
            out.write(json.dumps({"index":i,"n":n,"verdict":v})+"\n"); out.flush()
            ckpt.checkpoint()
    if ckpt and out: ckpt.checkpoint(force=True)
    if out: out.close()
    print("CENSUS:",tally)
elif argv[1]=="one":
    r=int(argv[2]); N=int(argv[3]); dn=int(argv[4])
    n=[int(x) for x in argv[5].split(",")]
    print(classify(r,N,dn,n))

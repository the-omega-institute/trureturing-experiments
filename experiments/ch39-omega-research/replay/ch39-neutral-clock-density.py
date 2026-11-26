"""Uniform quadratic bounds for the neutral one-loop potential on the entire clock torus."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import arb

pot = runpy.run_path(str(Path(__file__).with_name("ch39-neutral-clock-potential.py")))
edges,sh,q,r = (pot[k] for k in ("edges","sh","q","r"))
colors = {(i,j):d if d in (724,804) else 0 for i,j,d,_ in sh["edges"]}
witnesses = {d:[] for d in (724,804)}
for i,j,k in combinations(range(19),3):
    tri = ((i,j),(j,k),(i,k))
    if all(pair in edges for pair in tri):
        types = sorted(colors[pair] for pair in tri)
        if types in ([0,0,724],[0,0,804]):
            witnesses[types[-1]].append(tri)

def lower(z):
    diag = [z for _ in range(19)]
    weight = {}
    for (i,j),(_,ell,w) in edges.items():
        y = ell*z.sqrt()
        v = y/y.tanh()+w.upper()
        diag[i] += v
        diag[j] += v
        weight[i,j] = y/y.sinh()+w.lower()
    coeffs = []
    for triangles in witnesses.values():
        total = q(0)
        for ij,jk,ik in triangles:
            i,j = ij
            k = jk[1]
            total += 8*weight[ij]*weight[jk]*weight[ik]/(diag[i]*diag[j]*diag[k])
        coeffs.append(total.lower())
    return coeffs

def check():
    mult = 32*arb.pi()**2*r**4
    lo = q(1)
    coeffs = [q(0),q(0)]
    for _ in range(96):
        hi = lo*q(2).sqrt().sqrt()
        coeffs = [v+(hi*hi-lo*lo)*b/2 for v,b in zip(coeffs,lower(hi))]
        lo=hi
    m = pot["short"]["base"]["arc_budget"]()[1].upper()/2
    z0 = 2*m
    wsum = max(sum((1+w.upper()) for pair,(_,_,w) in edges.items() if colors[pair]==d)
               for d in witnesses)
    high = (8*z0*wsum+608*m**3/z0)/mult
    low = min(v.lower() for v in coeffs)/mult
    assert low>0
    print("single-clock neutral triangle counts",[len(v) for v in witnesses.values()])
    print("bounds per sin(p)^2+sin(q)^2",low.lower().str(12),high.upper().str(12))
    return low,high

if __name__ == "__main__":
    check()

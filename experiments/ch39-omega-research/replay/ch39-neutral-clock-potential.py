"""Positive one-loop cost of the two adjoint clocks, certified with closed triangles."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-su3-root-clock-lift.py")))
short,sh,q = m["short"],m["sh"],m["q"]
r = short["r"]
edges,tri = {},[]
weights = short["m"]["base"]["weights"]
for i,j,d,_ in sh["edges"]:
    ids = short["m"]["source_ids"]
    w = weights[ids[i,j]] if (i,j) in ids else q(0)
    f = short["fsbox"][short["band"](d)]
    clock = m["U"] if d==724 else m["V"] if d==804 else m["I"]
    ell = f.upper()*(sh["rho"]*sh["theta"](d)).sqrt()
    edges[i,j] = clock,ell,w

for i,j,k in combinations(range(19),3):
    pair = ((i,j),(j,k),(i,k))
    if not all(p in edges for p in pair):
        continue
    u,v,w = [edges[p][0] for p in pair]
    rot = u*v*w.conjugate().transpose()
    tr = rot[0,0]+rot[1,1]
    loss = 4-abs(tr)**2
    if loss>0:
        tri.append((pair,loss.lower()))

def lower(z):
    coeff = {}
    diag = [z for _ in range(19)]
    for (i,j),(_,ell,w) in edges.items():
        y = ell*z.sqrt()
        c = y/y.tanh()+w.upper()
        diag[i] += c
        diag[j] += c
        coeff[i,j] = y/y.sinh()+w.lower()
    val = q(0)
    for (ij,jk,ik),loss in tri:
        i,j = ij
        k = jk[1]
        val += 2*loss*coeff[ij]*coeff[jk]*coeff[ik]/(diag[i]*diag[j]*diag[k])
    return val.lower()

def check():
    assert tri
    total = q(0)
    # With z in M^2 units all positive closed-walk weights decrease in z.
    ratio = q(2).sqrt().sqrt()
    lo = q(1)
    for _ in range(96):
        hi = lo*ratio
        total += (hi*hi-lo*lo)*lower(hi)/2
        lo = hi
    energy = total/(32*arb.pi()**2*r**4)
    print("positive adjoint triangle witnesses",len(tri),"energy gap per real polarization / Q^4 >",energy.lower().str(12))
    assert energy>0

if __name__=="__main__":
    check()

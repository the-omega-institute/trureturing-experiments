"""Certified leading scalar mode for the quadratic part of the sixth spectral moment."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

folder=Path(__file__).resolve().parent
w=run_path(str(folder/"ch39-shell-gw-quartic-spectral-action.py"))
a,b=w["hops"]()
adj,hs=w["adj"],w["hs"]
count=w["model"]["model"]["count"]
eye=w["m"]["yuk"]["identity"]


def quadratic(left,right,degree):
    s=(degree-2)//2
    powers_a=[eye(18)]
    powers_b=[eye(18)]
    for _ in range(2*s):
        powers_a.append(powers_a[-1]*left)
        powers_b.append(powers_b[-1]*right)
    loc=[(3*i+j,3*i+l) for i in range(6) for j in range(3) for l in range(3)]
    return acb_mat([[degree*sum((-1)**r*powers_a[2*s-r][i,k]*powers_b[r][l,j]
                     for r in range(2*s+1)) for k,l in loc] for i,j in loc])


def rows(mat):
    return [[mat[i,j] for j in range(mat.ncols())] for i in range(mat.nrows())]


def compressed_yukawa(v):
    lib=w["m"]["yuk"]
    l,r=lib["zero_basis"](1),lib["zero_basis"](-1)
    out=acb_mat(3,3)
    for i in range(6):
        mat=acb_mat([[v[9*i+3*j+t,0] for t in range(3)] for j in range(3)])
        L=acb_mat([[l[3*i+j,t] for t in range(3)] for j in range(3)])
        R=acb_mat([[r[3*i+j,t] for t in range(3)] for j in range(3)])
        out+=w["phases"][i]*L.transpose()*mat.transpose()*R
    return out


if __name__=="__main__":
    h=quadratic(a,b,6)
    vals,vecs=np.linalg.eigh(np.array([[complex(z) for z in row] for row in rows(h)]))
    lo,hi = arb("161.055369642"),arb("161.055369644")
    assert count(rows(h),lo)==0 and count(rows(h),hi)==1
    assert count(rows(h),arb("161.15959"))==1
    # Shifted inverse iteration and the gap isolate a unit vector on the true lowest line.
    v=acb_mat([[acb(float(z.real),float(z.imag))] for z in vecs[:,0]])
    e=arb(float(vals[0])); shift=acb(0,arb("1e-7"));system=h-(e+shift)*eye(54)
    for _ in range(4):
        v=system.solve(v)*(-shift)
    v=v/hs(v,v).real.sqrt()
    e=hs(v,h*v).real/hs(v,v).real
    r=h*v-e*v
    # Normalizing the orthogonal projection moves v by at most twice its high component.
    margin=arb(0,(2*hs(r,r).real.sqrt()/(arb("161.15959")-e)).upper())
    v=acb_mat([[v[i,0]+acb(margin,margin)] for i in range(54)])
    singlet=sum(abs(sum(v[9*i+4*j,0] for j in range(3)))**2 for i in range(6))/3
    assert arb(".99894")<singlet<arb(".99896")
    y=compressed_yukawa(v);sq=adj(y)*y
    assert count(rows(sq),arb(".2334")**2)==0 and count(rows(sq),arb(".2352")**2)==3
    fraction=(arb(".016")/(arb("161.15959")-hi)).sqrt()
    gain=arb(".515")
    ratio=(arb(".2334")*(1-fraction**2).sqrt()-gain*fraction)/(arb(".2352")+gain*fraction)
    assert ratio>arb(".027")
    print("sextic onset",e.str(13),"singlet fraction",singlet.str(10),
          "within .016 of onset the light/heavy ratio exceeds",ratio.str(8))

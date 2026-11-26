"""Smallest support in actual zeta heights reproducing the four centered RH samples."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

record=run_path(str(Path(__file__).with_name("ch39-finite-rh-pivot-mass.py")))


def data():
    levels=[arb(1)]
    for j in (1,2,3):
        k=record["bins"][j][0]
        levels.append((arb(k).union(arb(k+1))/record["den"])**2)
    heights=[r.imag for r in acb.zeta_zeros(1,31)]
    zs=[1/x**2 for x in heights]
    indices=(0,1,12,14)
    mat=arb_mat([[1/(1+j*j*zs[i]) for i in indices] for j in range(4)])
    weights=mat.solve(arb_mat([[x] for x in levels])).entries()
    return levels,heights,indices,weights


def main():
    ctx.prec=192
    levels,heights,indices,weights=data()
    z=[1/g**2 for g in heights]
    W=arb_mat([[1,14,49,36],[1,13,36,0],[1,10,9,0],[1,5,4,0]])
    m=W.solve(arb_mat([[x] for x in levels])).entries()
    assert m[3]>z[1]*m[2]
    moments=[z[0]*m[j]-m[j+1] for j in range(3)]
    u=moments[1]/moments[0]
    v=moments[2]/moments[0]-u*u
    assert moments[0]>0 and v>0
    assert z[3]<u+v/u and z[2]>u+v/u
    for x,lo,hi in [(z[1],13,14),(z[2],28,29)]:
        y=u-v/(x-u)
        assert z[lo]>y and y>z[hi]
    assert v>(z[1]-u)*(u-z[13])
    assert min(weights)>arb(".1665") and sum(weights).contains(1)
    print("Four actual heights required, at least up to index 15.")
    print("Indices",[i+1 for i in indices],"weights",[x.str(12) for x in weights])


if __name__=="__main__":
    main()

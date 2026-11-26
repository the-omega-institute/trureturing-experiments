"""Quartic spectral action of the doubled family connection on the six shell sites."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-higgs-low-band.py")))
m["ctx"].prec=288
adj, hs, trace = m["adj"],m["hs"],m["trace"]
s,k,phases,gen = m["transported_spectra"]()
pts,model=m["pts"],m["model"]
h=acb.zeta_zero(43).imag
angles=[(1-arb(d)/26).acos() for d in (3,5,7)]
p,q=(h*(a-angles[0]) for a in angles[1:])


def hops():
    a,b=acb_mat(18,18),acb_mat(18,18)
    for i in range(6):
        for j in range(6):
            if j==i or j==(i^1):
                continue
            link=acb_mat(model["link"](pts[j],pts[i],p,q))
            for row in range(3):
                for col in range(3):
                    b[3*i+row,3*j+col]=link[row,col]
                    a[3*i+row,3*j+col]=-s[i,j]*link[row,col]
    return a,b


def diagonal(blocks):
    return acb_mat([[blocks[i//3][i%3,j%3] if i//3==j//3 else 0
                     for j in range(18)] for i in range(18)])


def full(a,b,x):
    return acb_mat([[a[i,j] if i<18 and j<18 else x[i,j-18] if i<18
                     else x[j,i-18].conjugate() if j<18 else -b[i-18,j-18]
                     for j in range(36)] for i in range(36)])


if __name__=="__main__":
    a,b=hops()
    for z in (a,b):
        zz=z*z
        for i in range(18):
            for j in range(3*(i//3),3*(i//3)+3):
                assert (zz[i,j]-4*int(i==j)).contains(0)
    blocks=[acb_mat([[acb(2*i+j-r,3*r+j-i)/11 for j in range(3)]
                      for r in range(3)]) for i in range(6)]
    x=diagonal(blocks);D=full(a,b,x);D0=full(a,b,acb_mat(18,18))
    mass=hs(x,x).real
    fourth=sum(hs(adj(v)*v,adj(v)*v).real for v in blocks)
    grad=a*x-x*b
    difference=trace(D**4-D0**4).real
    kinetic=hs(grad,grad).real/2
    assert (difference-(16*mass+2*fourth+4*kinetic)).contains(0)
    assert (trace(D**2-D0**2).real-2*mass).contains(0)
    count=model["model"]["count"]
    def rows(a):return [[a[i,j] for j in range(a.ncols())] for i in range(a.nrows())]
    assert count(rows(s),arb(".01377"))==1 and count(rows(k),arb(".01588"))==0
    print("quartic spectral action: 2*site quartic + 16*mass + 4*transport kinetic; singlet-only minimum")

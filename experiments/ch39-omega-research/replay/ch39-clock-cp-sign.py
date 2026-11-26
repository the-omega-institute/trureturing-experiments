"""Nonzero CP-odd witnesses on every edge-selection orbit at the same clock height."""
from pathlib import Path
import runpy
from flint import arb,acb,acb_mat

edge=runpy.run_path(str(Path(__file__).with_name("ch39-clock-edge-source.py")))
src=edge["src"]
shell=edge["s"]
W,J=src["W"],src["J"]
block=src["block"]


def signed(p,q,pair):
    m=shell["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    inv=block(m,J,J).inv()
    b=block(m,W,J)
    k=block(m,W,W)-b*inv*b.conjugate().transpose()
    i,j=pair
    v=-b*block(inv,range(len(J)),(J.index(i),))
    v[W.index(j),0]+=m[i,j].conjugate()
    e,a,b=src["fourier_columns"](k,v)
    t=arb(".001")
    y=acb_mat([[e[i]*(i==j)-t*a[i]*b[j].conjugate() for j in range(3)] for i in range(3)])
    hh=y*y.conjugate().transpose()
    d=[abs(z)**2 for z in e]
    # det[H0,H1] is a simultaneous-unitary invariant; conjugation reverses its sign.
    return 2*(d[0]-d[1])*(d[1]-d[2])*(d[2]-d[0])*(hh[0,1]*hh[1,2]*hh[2,0]).imag


if __name__=="__main__":
    p,q=edge["p"],edge["q"]
    pairs=[(0,7),(0,8),(1,6),(1,8)]
    sites=shell["sites"]
    reflect=lambda i:sites.index((sites[i][1],sites[i][0],sites[i][2]))
    for pair in pairs:
        value=signed(p,q,pair)
        reflected=signed(-p,-q,tuple(map(reflect,pair)))
        assert not value.contains(0) and (value+reflected).contains(0)
        print("clock edge, imaginary commutator determinant",pair,value.str(12))

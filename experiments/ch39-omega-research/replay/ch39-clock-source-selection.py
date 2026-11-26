"""Minimizers for fixed total diagonal source in the original Gaussian determinant."""
from pathlib import Path
import runpy
from flint import arb,acb

source=runpy.run_path(str(Path(__file__).with_name("ch39-single-arc-source.py")))
s=source["shell"]
sites=s["sites"]
rotate=lambda i:sites.index(sites[i][1:]+sites[i][:1])
reflect=lambda i:sites.index((sites[i][1],sites[i][0],sites[i][2]))
even=sorted({0,rotate(0),rotate(rotate(0))})
odd=sorted(reflect(i) for i in even)
assert even==[0,3,4] and odd==[1,2,5]


def covariance(p,q):
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    for l in range(1,20):
        assert source["block"](m,range(l),range(l)).det().real>0
    return m.inv()


if __name__=="__main__":
    p,q=(s[key]*source["ref"] for key in ("a","b"))
    fixed=covariance(p,q)
    winners=[min(group,key=lambda i:float(fixed[i,i].real.mid())) for group in (range(6,9),range(19))]
    assert winners[0] in range(6,9)
    assert min(fixed[i,i].real.lower() for i in range(6))>fixed[6,6].real
    assert all(fixed[14,14].real<fixed[i,i].real for i in range(19) if i!=14)
    radius=arb(0,arb(".00001"))
    puls=covariance(p+radius,q+radius)
    assert all(puls[i,i].real<puls[j,j].real for i in even for j in odd)
    # Matrix coefficients at reflected sites conjugate after negating both clock angles.
    # Transport without a clock conjugates under a coordinate transposition (signed determinant).
    m=s["matrix"](0,{724:acb(0,p).exp(),804:acb(0,q).exp()})
    mirror=s["matrix"](0,{724:acb(0,-p).exp(),804:acb(0,-q).exp()})
    assert all((m[i,j].conjugate()-mirror[reflect(i),reflect(j)]).contains(0)
               for i in range(19) for j in range(19))
    gap=fixed[odd[0],odd[0]].real-fixed[even[0],even[0]].real
    print("selected three clock vertices, diagonal difference, unrestricted winner",even,gap.str(15),14)

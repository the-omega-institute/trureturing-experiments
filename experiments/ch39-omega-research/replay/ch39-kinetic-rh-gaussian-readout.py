"""Positive reference kinetic poles required by finite Euclidean RH response readings."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_poly, arb, arb_mat, arb_poly, ctx

directory=Path(__file__).parent
read=run_path(str(directory/"ch39-finite-nbb-coherent-four-windows.py"))["readings"]


def gauss(count):
    """count nodes reproduce the normalization and 2*count-1 positive offsets."""
    last=2*count-1
    levels=read(last+2)
    # A dyadic scale keeps the low ordinate coordinates moderate.
    unit=256
    x=arb_poly([0,1])
    polys=[1+arb(j*j)/unit*x for j in range(1,last+1)]
    W=prod(polys,start=arb_poly([1]))
    bases=[W]+[prod(polys[:i]+polys[i+1:],start=arb_poly([1])) for i in range(last)]
    mat=arb_mat([[poly[k] for k in range(last+1)] for poly in bases])
    mu=mat.solve(arb_mat([[v] for v in levels[:last+1]])).entries()
    H=arb_mat([[mu[i+j] for i in range(count)] for j in range(count)])
    coeff=H.solve(arb_mat([[-mu[count+i]] for i in range(count)])).entries()
    polynomial=acb_poly(coeff+[1])
    roots=sorted([z.real for z in polynomial.roots(arb(2)**(-ctx.prec//3))],
                 key=lambda z:float(z.mid()),reverse=True)
    assert len(roots)==count and all(g>0 for g in roots)
    assert all(a>b for a,b in zip(roots,roots[1:]))
    # The polynomial is real and the disjoint real-part enclosures cover
    # all roots; a nonreal root would require a conjugate in its own strip.
    mat=arb_mat([[u**j for u in roots] for j in range(count)])
    bare=mat.solve(arb_mat([[v] for v in mu[:count]])).entries()
    weights=[a*W(u) for a,u in zip(bare,roots)]
    assert all(w>0 for w in weights) and sum(weights).contains(1)
    z=[r/unit for r in roots]
    for j in range(last+1):
        assert (levels[j]-sum(w/(1+j*j*v) for w,v in zip(weights,z))).contains(0)
    nextj=last+1
    gap=levels[nextj]-sum(w/(1+nextj**2*v) for w,v in zip(weights,z))
    assert gap>0
    return z,weights,gap


def main():
    ctx.prec,ctx.cap=3072,3
    for count in (2,21):
        z,w,gap=gauss(count)
        average=sum(z)/count
        kinetic=[x/average for x in z]
        print("Reference dimension",count,"minimum spectral weight",min(w).str(12))
        print("Effective frequency extrema",(1/z[0].sqrt()).str(12),(1/z[-1].sqrt()).str(12),
              "largest mean-normalized inertia",max(kinetic).str(12))
        print("Next full RH squared overlap exceeds finite response by",gap.str(15))


if __name__=="__main__":
    main()

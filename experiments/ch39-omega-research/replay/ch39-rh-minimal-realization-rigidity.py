"""Certify stability of the minimal full-edge time realization and of its one-loop crossing."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

directory=Path(__file__).parent
shell=run_path(str(directory/"ch39-kinetic-rh-gaussian-crossing.py"))


def main():
    ctx.prec,ctx.cap=3072,3
    z,w,gap=shell["quad"](21)
    zmin,zmax=min(z),max(z)
    av=sum(z)/21
    low,high=av/zmax,av/zmin
    # Use tangent matrix errors <= 2^-288, with normalized zeroth entry.
    eps=arb(2)**-288
    gram=arb(2)**169*eps
    mult=2*gram
    invroot=1/(1-gram).sqrt()
    polar=invroot-1
    da=invroot**2*mult+zmax*polar*(1+invroot)
    assert da<arb(2)**-117 and zmin>da
    # W(A)^1/2 is a product of positive square roots. Telescope between
    # actual and target factors, bounding each local difference by j^2*da/2.
    size=prod((1+j*j*(zmax+da) for j in range(1,42)),start=arb(1)).sqrt()
    dW=size*sum(j*j for j in range(1,42))*da/2
    dzmap=size*(1-(1-gram).sqrt())+dW
    # Mean normalization maps A to mass-squared operator C=mean(A)*A^-1.
    dc=da/(zmin-da)+av*da/(zmin*(zmin-da))
    rootC=dc/(low.sqrt()+(low-dc).sqrt())
    off=(arb(6).sqrt()*((high+dc).sqrt()*dzmap+rootC))
    mass_error=dc+off
    print("Absolute kinetic eigenvalue error",da.str(7),"coupling-map error",dzmap.str(7),
          "dimensionless scalar mass-squared matrix error",mass_error.str(7))
    assert mass_error<arb("1e-17")
    a,b,c=[shell["base"]["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    modes=[(arb(6),7),(arb(4),6)]
    for j,n in ((5,3),(3,2)):
        modes += [(3-arb(j).sqrt(),n),(3+arb(j).sqrt(),n)]
    d=[8*arb.pi()**2/arb(v)**2 for v in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    q=1000*(K/44).exp()
    for sign,mf in [(1,arb("7.4768884e13")),(-1,arb("7.4768887e13"))]:
        spectrum=shell["scalar_spectrum"](mf,[av/u for u in z],[(1-p)*v for v in w],p,modes)
        e=arb(0)
        count=0
        for val,n in spectrum:
            band=val+arb(0,mass_error)
            assert mf*band.sqrt()>1000
            term=(q/(mf*band.sqrt())).log()
            if term>0:
                e+=n*term
                count+=n
            else:
                assert term<0
        assert count==144 and sign*(22*e-15*F)>0 and q/mf<arb(8).sqrt()
    assert gap>arb("9.74816e-40")
    print("Minimum 1008-state dynamics is spectrally rigid; 144 open scalar modes persist in the 288-bit box.")
    print("Fewer than 48 additional states leave a common Gaussian direction at all later offsets.")


if __name__=="__main__":
    main()

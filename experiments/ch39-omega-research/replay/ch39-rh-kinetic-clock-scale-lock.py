"""Test a common rescaling of all true kinetic poles with arbitrary positive weights."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

directory=Path(__file__).parent
cross=run_path(str(directory/"ch39-kinetic-rh-gaussian-crossing.py"))


def main():
    ctx.prec,ctx.cap=1536,3
    gauss,weights,_=cross["quad"](21)
    true=[1/acb.zeta_zero(k).imag**2 for k in range(1,112)]
    windows=[]
    for index,count,m,rel in [(0,7,8,"0.000000087"),(1,6,9,"0.000099"),(2,6,9,"0.0026")]:
        center=true[index]
        degree=2*(count-1+m)+2
        nodes=[(1+round(40*i/degree))**2 for i in range(degree+1)]
        assert len(set(nodes))==degree+1
        def polynomial(x):
            q=x**m*prod(x-v for i,v in enumerate(true[:count]) if i!=index)
            return q*q*(arb(rel)**2-(x/center-1)**2)
        def W(x):
            return prod(1+s*x for s in nodes)
        target=sum(w*polynomial(n)/W(n) for n,w in zip(gauss,weights))
        norm=arb(0)
        for s in nodes:
            p=-arb(1)/s
            norm+=abs(polynomial(p)/prod(1+v*p for v in nodes if v!=s))
        assert target>arb(2)**-96*norm
        low,high=center*(1-arb(rel)),center*(1+arb(rel))
        windows.append((low,high))
    # Fitting data within tolerance forces at least one support point into
    # EACH window. A uniform dial then has to satisfy three finite lists.
    max_scale=arb(64)
    choices=[]
    for low,high in windows:
        assert max_scale*true[-1]<low
        choices.append([(i+1,low/v,high/v) for i,v in enumerate(true) if low/v < max_scale])
    possible=[]
    from itertools import product
    for a,b,c in product(*choices):
        lo=a[1].max(b[1]).max(c[1])
        hi=a[2].min(b[2]).min(c[2]).min(max_scale)
        if lo>hi:continue
        assert lo<hi
        possible.append((a[0],b[0],c[0]))
    print("Allowed index triples below scale 64",possible)
    assert possible==[(1,2,3)]
    print("On 0 < scale <= 64, all 96-bit fits on scaled true heights require |scale-1| < 8.7e-8.")


if __name__=="__main__":
    main()

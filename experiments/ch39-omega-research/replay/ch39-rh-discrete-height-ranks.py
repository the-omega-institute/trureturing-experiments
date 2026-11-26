"""Positive three-axis weight at the first true heights, for 42 calibrated readings."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

cross=run_path(str(Path(__file__).with_name("ch39-kinetic-rh-gaussian-crossing.py")))


def main():
    ctx.prec,ctx.cap=1536,3
    gauss,weights,_=cross["quad"](21)
    true=[1/acb.zeta_zero(i).imag**2 for i in range(1,42)]
    cases=[(2,2),(3,4),(4,5),(5,8),(6,9),(7,11),(10,10),(9,13),
           (12,14),(13,15),(16,13),(13,21),(20,15)]
    for i,(count,power) in enumerate(cases):
        degree=count-1+power
        js=sorted({1+round(40*j/min(40,degree)) for j in range(min(40,degree)+1)})
        nodes=[j*j for j in js]
        def polynomial(x):
            return (-1)**i*x**power*prod(x-z for k,z in enumerate(true[:count]) if k!=i)
        def denominator(x):
            return prod(1+s*x for s in nodes)
        assert (count-1-i)%2==1 and degree<=len(nodes) and max(js)<=41
        # The target below is evaluated against the Gaussian rule, which has
        # exactly the same readings at all the poles used by this rational filter.
        target=sum(w*polynomial(x)/denominator(x) for x,w in zip(gauss,weights))
        norm=arb(1)/prod(nodes) if len(nodes)==degree else arb(0)
        for s in nodes:
            v=-arb(1)/s
            norm+=abs(polynomial(v)/prod(1+b*v for b in nodes if b!=s))
        value=polynomial(true[i])/denominator(true[i])
        assert value>0 and target>0
        precision=96 if i<8 else 128 if i<11 else 192
        bound=(target-norm*arb(2)**(-precision))/value
        assert bound>0
        print(i+1,"absolute bits",precision,"guaranteed minimum population",bound.str(9))


if __name__=="__main__":
    main()

"""An exact finite partition sum for a nearest-step bias among the ternary zero fields."""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import arb,ctx
from sympy import symbols, expand

directory=Path(__file__).parent
stability=run_path(str(directory/"ch39-rh-heap-finite-stability.py"))
arc=run_path(str(directory/"ch39-rh-internal-arc-calibration.py"))


def main():
    coefficients=[Counter(),Counter()]
    for field in stability["representatives"]():
        triad=bool(stability["classify"](field)[0])
        step=sum(1-stability["dot"](field[p],field[q]) for p,q,r in stability["trails"])
        coefficients[triad][step]+=192 if triad else 8
    u=symbols("u")
    polynomials=[]
    for table in coefficients:
        assert all(e%32==0 for e in table)
        polynomials.append(sum(n*u**(e//32) for e,n in table.items()))
    c,t=polynomials
    assert expand(t-1152*u**3*(1+u)**3)==0
    assert expand(c-8*((1+u**3)**3+3*u**3*(1+u)**3))==0
    ctx.prec,ctx.cap=1536,3
    data=arc["fixed"]["inputs"]()
    a=arc["argument_for"](arb(33)/34,data)
    d=arc["argument_for"](arb(14)/17,data)/(3*a)
    share=(33*d-89)/64
    assert arb(1)/37<share<1
    h=(144*share/(1-share)-3).root(3)
    u=2/(1+h+((1+h)**2-4).sqrt())
    rate=u.root(32)
    print("Required collinear fraction",share.str(14),"nearest-step weight",rate.str(14))
    for numerator in (18,21,25,27,28,31,33):
        x=arb(numerator)/34
        factor=x**32+1/x**32-1
        theta=(3+factor**3)/(147+factor**3)
        predicted=(89+64*theta)/33
        assert abs(predicted-d)>arb("0.19")
    print("No single nontrivial shell half-angle provides this microstate weight.")
    neutral=arc["response"](a*9/11,data)
    print("Unit base and unit stiffness give neutral response",neutral.str(14))


if __name__=="__main__":
    main()

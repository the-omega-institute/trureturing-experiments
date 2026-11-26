"""Momentum-dependent Dirac threshold moments and the amplified zero height."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import cancel, symbols


def identities():
    weights = [(F(1,5),F(3,10),F(4,5),F(1,10),F(3,5)),
               (0,F(1,2),0,F(3,2),0),(F(1,2),0,F(1,2),1,0)]
    assert [5*u-3*v-2*w for u,v,w in zip(*weights)]==[0,0,3,-6,3]
    m,n,s = symbols("m n s",nonzero=True)
    # Let |C|=m, |E|=n, t^2=mn and d^2=(m-n)^2/2.
    tau, dsq = m*n,(m-n)**2/2
    frac = (1+s/m**2)*(1+s/n**2)/(1+s/tau)**2
    assert cancel(frac-1-2*dsq*s/(tau+s)**2)==0
    b1,b2 = [3*(m**(-2*j)+n**(-2*j)-2*tau**(-j)) for j in (1,2)]
    assert cancel(b1-6*dsq/tau**2)==0
    assert cancel(b2-4*b1/tau-b1**2/3)==0
    assert cancel(2*b1**2/(3*b2-b1**2)-dsq/tau)==0


def verify():
    identities()
    ctx.prec = 160
    match = run_path(str(Path(__file__).with_name("ch39-quartic-soft-adjoint-match.py")))
    a = match["a"]
    rows = []
    for n,kind,r,v,state in match["matched_points"]():
        if kind!="center":
            continue
        gamma = acb.zeta_zero(n).imag
        d = r*(a("13/10")/state[0]).sqrt()
        ell2 = ((gamma/624)**2+a("1/4")).sqrt()
        middle = d*acb(((ell2+a("1/2"))/2).sqrt(),((ell2-a("1/2"))/2).sqrt())
        c,t,e = [abs(u)*v for u in (middle+d,middle,middle-d)]
        b1,b2 = [3*(c**(-2*j)+e**(-2*j)-2*t**(-2*j)) for j in (1,2)]
        curvature = 3*b2/(b1*b1)
        recovered = 312*((curvature-1)**2-1).sqrt()
        assert b1>0 and b2>0 and (recovered-gamma).contains(0)
        predicted = 1+(1+(gamma/312)**2).sqrt()
        assert (curvature-predicted).contains(0)
        rows.append(predicted)
        print(n,"normalized moment ratio",predicted.str(13))
    gap = min((b-a).lower() for a,b in zip(rows,rows[1:]))
    assert gap>arb(".00148") and gap<arb(".00149")
    print("Smallest normalized moment gap",gap.str(10))


if __name__ == "__main__":
    verify()

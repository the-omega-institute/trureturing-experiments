"""Invert three heavy fermion moments after canceling every bridge mass block."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, acb, ctx
from sympy import cancel, symbols


def identities():
    weights = [(F(1,5),0,F(1,2)),(F(3,10),F(1,2),0),
               (F(4,5),0,F(1,2)),(F(1,10),F(3,2),1),(F(3,5),0,0),
               (F(4,5),0,F(1,2)),(F(7,10),F(3,2),1)]
    assert [5*x+3*y-8*z for x,y,z in weights]==[-3,3,0,-3,3,0,0]
    p,q = symbols("p q", positive=True)
    x,y = (1+p)/(1-p),(1+q)/(1-q)
    m = [(1-p**j)*(1-q**j) for j in (1,2,3)]
    H,J = m[1]/m[0]**2,m[2]/m[0]**3
    assert cancel(H-x*y)==0
    assert cancel(16*J-9*H**2-1-3*(x*x+y*y))==0
    assert cancel((x*x-1)*(1-p)**2-4*p)==0


def recover(moments):
    t,u,v = moments
    H,J = u/t**2,v/t**3
    total = (16*J-9*H**2-1)/3
    big = (total+(total*total-4*H**2).sqrt())/2
    return H,big,(big*(big-2)).sqrt()/2


def verify():
    identities()
    ctx.prec = 180
    match = run_path(str(Path(__file__).with_name("ch39-quartic-soft-adjoint-match.py")))
    a = match["a"]
    readings = []
    for n,kind,r,v,state in match["matched_points"]():
        if kind!="center":
            continue
        gamma = acb.zeta_zero(n).imag
        d,h = [v*r*(a(t)/state[0]).sqrt() for t in ("13/10","7/2")]
        ell2 = ((gamma/624)**2+a("1/4")).sqrt()
        D = d*acb(((ell2+a("1/2"))/2).sqrt(),((ell2-a("1/2"))/2).sqrt())
        C,E = D+d,D-d
        P,W = -h*D/d,-h*C/d
        moments = [abs(W)**(-2*j)-abs(P)**(-2*j)+abs(E)**(-2*j)-abs(D)**(-2*j)
                   for j in (1,2,3)]
        H,big,eff = recover(moments)
        assert (eff-gamma/624).contains(0)
        assert (big-1-(1+4*(gamma/624)**2).sqrt()).contains(0)
        readings.append(H)
        print(n,"second moment divided by first moment squared",H.str(12))
    gaps = [y-x for x,y in zip(readings,readings[1:])]
    assert all(g>0 for g in gaps)
    print("Minimum gap",min(g.lower() for g in gaps).str(10))


if __name__ == "__main__":
    verify()

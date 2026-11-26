"""Bounds uniform over edge magnitudes, and exact winding sums for a retuned unweighted graph."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
RATIO = run_path(str(ROOT/"ch39-quartic-adjoint-ratio.py"))
a = RATIO["a"]


def equal_edge_return(t,u):
    omega = acb(-1,arb(3).sqrt())/2
    val = acb(0)
    for j in range(3):
        x = 1/(t*omega**j)
        denom = x**4-5*x*x-4*x*u+4*(1-u*u)
        val += omega**(-j)*x*(4*x**3-10*x-4*u)/denom
    return val.real+6*t**4/(1-t**6)


def equal_edge_weights(t, m):
    omega = acb(-1,arb(3).sqrt())/2
    left = acb(0)
    right = acb(0)
    for j in range(3):
        x = 1/(t*omega**j)
        for sign in (1,-1):
            u = (-x+sign*(x*x-2))/2
            s = (u**2-1).sqrt()
            outer, = [root for root in (u+s,u-s) if abs(root)>1]
            tau = 1/outer
            residue = omega**(-j)*x*(4*x**3-10*x-4*u)/(-4*x-8*u)
            leading = -2*residue*tau**(m+1)/(1-tau**2)
            left += leading*(1+tau+tau**2)
            right += leading/(1-tau)
    assert left.imag.contains(0) and right.imag.contains(0)
    return 2*left.real,2*right.real


def verify():
    ctx.prec = 260
    for n,height,eps in RATIO["high_zero_readings"]():
        cosine = 1-1/(2*(height**2+a("1/4")))
        # A Perron vector of the entrywise magnitudes on a component tests B.
        q = eps.sqrt().sqrt()/cosine
        norm = 19*q**265/(1-q**3)
        symmetric = (eps/3).sqrt().sqrt()/cosine
        repeated = 12*symmetric**265/(1-symmetric**3)
        assert q<a(".083") and norm<a("1e-285") and repeated<a("1e-315")
        lo,hi = F("0.01"),F("0.04")
        assert equal_edge_return(a(lo),cosine)<eps<equal_edge_return(a(hi),cosine)
        for _ in range(125):
            mid = (lo+hi)/2
            value = equal_edge_return(a(mid),cosine)-eps
            assert value>0 or value<0
            if value>0:
                hi=mid
            else:
                lo=mid
        t = a(lo).union(a(hi))
        m = 131 if n<=11486 else 132
        low,high = equal_edge_weights(t,m)
        assert low>0 and high-low>0 and (high-low)/low<a("1e-8")
        print(n,"arbitrary, repeated bounds",norm.str(6),repeated.str(6),
              "equal hopping",t.str(10),"compatible weight",low.union(high).str(9))


if __name__=="__main__":
    verify()

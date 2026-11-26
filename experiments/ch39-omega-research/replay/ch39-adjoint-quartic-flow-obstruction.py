"""Rational fourth moments and the one-loop obstruction in the real SU(5) adjoint."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import Rational, diag, expand, eye, ones, symbols

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
a = MATCH["a"]


def traces():
    z = list(symbols("d0:4"))
    z.append(-sum(z))
    t2, t4 = (sum(d**k for d in z) for k in (2,4))
    vectors = sum((x-y)**4 for x in z for y in z)
    exterior = sum((x+y)**4 for i,x in enumerate(z) for y in z[i+1:])
    assert expand(vectors-10*t4-6*t2**2) == 0
    assert expand(exterior+3*t4-3*t2**2) == 0
    # Hessian of Tr Phi^4 in the real root planes and traceless diagonal subspace.
    roots = sum(4*(x*x+x*y+y*y)**2 for i,x in enumerate(z) for j,y in enumerate(z) if i!=j)
    projector = eye(5)-ones(5)/5
    diagonal = 6*projector*diag(*(x*x for x in z))*projector
    assert expand(roots+(diagonal**2).trace()-Rational(128,5)*t4-Rational(336,25)*t2**2) == 0
    # The only SU(5)-charged Weyl multiplets have center classes -1 and 2.
    assert all((s+t)%5 for s in (-1,2) for t in (-1,2))


def radial_bound():
    n = 24
    # Axes and normalized sign vectors give a positive finite fourth-moment rule.
    axis_weight, sign_weight = F(2,n+2),F(n,n+2)
    assert axis_weight/n+sign_weight/n**2 == F(3,n*(n+2))
    assert sign_weight/n**2 == F(1,n*(n+2))
    moment = F(MATCH["SCATTER"]["group_traces"](),n*(n+2))
    assert moment == F(47,520)
    gauge = 6*(10*moment+F(6,4))
    assert gauge == F(375,26)
    assert 2*(n+8)-F(60**2,4)/gauge == F(8,5)
    assert gauge-F(60**2,4*64) == F(75,208)
    # Integrate a lower barrier as long as the scalar flow exists below the target.
    low = a(0)
    target, z0, step = a(".38"), a("1.16"), a("1/64")
    for j in range(627):
        ahead = low+a(".0025")
        assert 128*ahead < 60*z0 and 60*ahead < 2*a(gauge)*z0
        slope = (64*ahead**2-60*z0*ahead+a(gauge)*z0**2)/(16*arb.pi()**2)
        assert 0 < slope*step < a(".0025")
        low = (low+step*slope).lower()
        if low > target:
            break
    assert low > target
    reached = (j+1)*step
    assert 10*arb.pi()**2*(1/target-26/(8*arb.pi())) < 159
    for _,_,_,v,state in MATCH["matched_points"]():
        assert 1/state[0] > z0 and (a(MATCH["PLANCK"]["TARGET"])/v).log() > reached
    print("Double-trace source: 9 g5^4 + 30 lambda_XPhi^2 + 168 e^2/25.")
    print("Certified radial lower bound at the Planck scale",low.str(8),"by t",reached.str(8))


if __name__ == "__main__":
    ctx.prec = 150
    traces()
    radial_bound()

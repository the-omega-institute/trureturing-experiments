"""Bracket a three-threshold Newton match using heat-trace interlacing for the neutral graph."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

shell = runpy.run_path(str(Path(__file__).with_name("ch39-curvature-spectrum-rigidity.py")))
q, r, A = (shell[k] for k in ("q", "r", "weyl"))
pi = arb.pi()
edge_count = len(shell["shell"]["shell"]["edges"])
vertices = sum(shell["ranks"])
assert edge_count == 135 and vertices == 19
# The subspace with all dynamical vertex values zero has codimension V.
# Each Dirichlet edge heat sum is between its Weyl integral minus 1/2
# (Jacobi inversion) and the Weyl integral (decreasing Gaussian).
L = 4*pi.sqrt()/3 * (1-4*q(2).sqrt()+3*q(3).sqrt())
B = 3*q(3).log()-4*q(2).log()
d = q(1)/12
scale = (1000*q("26.8220").exp()).union(1000*q("27.0455").exp())
mp = q("1.220890e19")
baseline = 2*q(".05").union(q("1.99779"))/(3*pi) + q(".000048797").union(q(".00030629"))/(4*pi)
target = pi/d*((mp/scale)**2-baseline)

def moment(nu):
    leading = A*L*nu*nu.sqrt()
    return leading + (-q(edge_count)/2).union(q(vertices))*B*nu

def bracket():
    w = (target.lower()/(A*L)).root(3)**2
    v = (target.upper()/(A*L)).root(3)**2
    # Initial enclosure and monotone bisection on the two scalar envelopes.
    lo = q(0)
    hi = w
    for _ in range(60):
        mid = (lo+hi)/2
        if moment(mid).upper()<target.lower():
            lo=mid
        else:
            hi=mid
    left = lo.lower()
    lo = v
    hi = 2*v
    assert moment(hi).lower()>target.upper()
    for _ in range(60):
        mid = (lo+hi)/2
        if moment(mid).lower()>target.upper():
            hi=mid
        else:
            lo=mid
    right = hi.upper()
    assert moment(left).upper()<target.lower()
    assert moment(right).lower()>target.upper()
    return left.union(right)

def check():
    assert 2*(1+2+1)+12 == 20 and (1+2+1)+12 == 16
    assert 12*(F(1,6)-F(1,6)) + F(12,6) - F(1,2) == F(3,2)
    assert 2*(F(1,6)-F(1,6))+F(1,6) == F(1,6)
    assert 2*(F(1,6)-F(1,3))+F(1,6) == -F(1,6)
    # UV scalar square weights: twelve conformal massless-family towers,
    # four shifted conformal scalar towers, and four shifted xi=1/3 towers.
    raw = F(2,36)+F(4,72)+2*(F(5,900)+F(2,72))
    slope = A*q(raw)/(16*pi**2)
    nu = bracket()
    masses = scale*nu.sqrt()
    heavy = q(3).sqrt()*masses
    center_error = (q(edge_count)*B/(2*A*L*nu.lower().sqrt())).upper()
    for label,val in (("threshold nu",nu),("light offset GeV",masses),("heavy offset GeV",heavy),
                      ("heavy offset / Planck",heavy/mp),("R^2 UV slope",slope)):
        print(label,val.lower().str(14),val.upper().str(14))
    print("relative error of leading Weyl moment",center_error.str(14))

if __name__ == "__main__":
    check()

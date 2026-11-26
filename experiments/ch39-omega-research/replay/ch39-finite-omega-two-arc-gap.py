"""Sharp three-reading discrepancy for finite subgroups of the unit quaternions."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from sympy import Rational, sqrt

heap = run_path(str(Path(__file__).with_name("ch39-finite-omega-heap-window.py")))
mul,bar = heap["times"],heap["bar"]


def main():
    a = (0,1,3,4)
    neighbors = [(0,1,4,3),(1,0,3,4)]
    u,v = [tuple(Q(k,26) for k in mul(a,bar(p))) for p in neighbors]
    comm = mul(mul(mul(u,v),bar(u)),bar(v))
    defect = (1-comm[0])/2
    assert u[0]==v[0]==Q(25,26) and defect==Q(25,4394)
    lower_scalar = u[0]-2*defect
    assert lower_scalar>0 and 4*(1-lower_scalar**2)<1
    # A 24th root and its inverse have the same first two readings and commute.
    # Their common error on each single step is below the commutator defect.
    qreal = (sqrt(6)+sqrt(2))/4
    assert abs((1+qreal)/2-Rational(51,52))<Rational(defect)
    print("Two shortest generators and the commutator:",u,v,comm)
    print("Sharp max error on these three fidelities, over finite unit-quaternion groups:",defect)


if __name__=="__main__":
    main()

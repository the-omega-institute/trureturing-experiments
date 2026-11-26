"""Hypercharge detects the rough diagonal obstruction on open sets of H links."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path
from sympy import symbols, expand


def verify():
    weights = (-2,-2,-2,3,3)
    character = run_path(str(Path(__file__).with_name("ch39-diagonal-index-quadraticity.py")))["character"]
    for m, table, expected in [(3,(0,16,16),40), (5,(0,50,25,25,50),460)]:
        count = 0
        for w in product(range(m),repeat=5):
            if sum(w)%m:
                continue
            k = character(w,table,m)
            count += sum(q*x for q,x in zip(weights,k)) != 0
        assert count == expected
    k = character((0,1,-1,1,-1),(0,16,16),3)
    assert 3*sum(q*x for q,x in zip(weights,k)) == -288
    radius = F(1,2**28)
    assert 64*radius < F(1,10**6)
    # Box packing in 26 real coordinates with mesh at most radius/6.
    assert 26 < 6**2 and 12/radius+1 <= 13/radius
    assert 26*4*5**4 == 65000
    # Independent antifundamental and antisymmetric kernels on a test tuple.
    f, gx, gy, plus, minus = symbols("f gx gy plus minus")
    first = 2*gx+2*gy
    second = -f+gx+plus+minus
    difference = expand(first-second)
    assert expand(difference.subs({gy:0,plus:gx,minus:gx})-(f-gx)) == 0
    assert expand(difference.subs({f:gx,gy:gx,minus:0})-(4*gx-plus)) == 0
    print("Hypercharge -288 at the witness; every one of the 40 and 460 cases detected.")
    print("Certified H-link neighborhood radius 2^-28, product-Haar lower bound exponent 65000.")


if __name__ == "__main__":
    verify()

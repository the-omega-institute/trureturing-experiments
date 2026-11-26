"""Gauge costs for minimal and symmetric charged scalar completions above Q."""
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq

model = run_path(str(Path(__file__).with_name("ch39-shell-gw-proton-abelian-lift.py")))
q = arb


def symmetric_frame():
    # Map X_il to sum_(ij,kl) P_sym[ij,kl] X_il conjugate(X_kj).
    for i in range(9):
        for j in range(9):
            for k in range(9):
                for l in range(9):
                    sym = fmpq(int(i == k and j == l)+int(i == l and j == k), 2)
                    identity = fmpq(int((i,l) == (k,j)), 2)
                    trace = fmpq(int(i == l and k == j), 2)
                    assert sym == identity+trace
    assert 9*(9+1)//2 == 45


def running():
    Q = q("4.352e13").union(q("4.473e13"))
    inv = q("38.82").union(q("38.98")) / (4*arb.pi())
    planck = q("2.435e18")
    t = (planck/Q).log() / (8*arb.pi()**2)
    ky = 5*inv/3
    for n in (17, 25, 30, 45):
        y = ky - q(47+n)*t/3
        assert y > 0
        minimum = 267*t + (16*t)**2/y
        level = minimum/inv
        limit = 1/minimum.sqrt()
        low, high = map(float, (level.lower(), level.upper()))
        print("charged species", n, "critical kX", low, high,
              "max gX in", float(limit.lower()), float(limit.upper()))
        assert (12 < level < 13) if n <= 30 else (13 < level < 14)
        xx = 100 - 267*t
        cross = -16*t
        det = xx*y-cross**2
        charges = [(q(a.numerator)/a.denominator, q(b)) for _,a,b,_,_ in model["fields"]]
        charges += [(q(1)/2, q(0)), (q(0), q(9)), (q(1), q(0))]
        size = max(((xx*a*a-2*cross*a*b+y*b*b)/det).upper() for a,b in charges)
        print("largest abelian charge norm squared at gX=.1 below", size.str(8))
    for kx in (None, 13, 14):
        xx = (q(100) if kx is None else kx*inv)-267*t
        assert xx > 0
        allowed = 3*(ky/t-256*t/xx)-47
        print("gX=.1" if kx is None else "kX="+str(kx),
              "critical number of charged species", float(allowed.lower()), float(allowed.upper()))
    # Landau constraint from hypercharge alone applies for every gX.
    absolute = 3*ky/t-47
    assert absolute < 66
    print("hypercharge-only ceiling", float(absolute.lower()), float(absolute.upper()))


if __name__ == "__main__":
    symmetric_frame()
    running()

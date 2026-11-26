"""Exact scale dependence and certified finite-resolution Planck matches."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import Poly, cancel, diff

ROOT = Path(__file__).parent
RATIO = run_path(str(ROOT/"ch39-quartic-adjoint-ratio.py"))
a, Q = RATIO["a"], RATIO["Q"]
TARGET = 12209*10**15
HALF_UNIT = 5*10**14


def mass_exponents(products):
    powers = {p: sum(v*c for v,c in zip((-220, 648, -428), indices))/504
              for p, indices in products[1].items()}
    powers[Q+2] += F(1, 2)
    assert powers == {Q: F(-109,14), Q+1: F(-107,42),
                      Q+2: F(369,14), 2*Q+1: F(-1,42)}
    cubic = 1344*Q**3+491*Q**2-1400*Q-654
    derivative = sum(c*diff(p,Q)/p for p,c in powers.items())
    assert cancel(derivative-cubic/(42*Q*(Q+1)*(Q+2)*(2*Q+1))) == 0
    assert Poly(cubic,Q).count_roots(0, 4) == 1
    return tuple(powers[p] for p in (Q,Q+1,Q+2,2*Q+1))


def mass_ratio(r, exponents):
    q = r*r
    normalized = (q, (q+1)/2, (q+2)/3, (2*q+1)/3)
    return sum(a(e)*x.log() for e,x in zip(exponents, normalized)).exp()


def enclose_root(mass, target, lo, hi, increasing):
    def side(r):
        residual = mass(a(r))-target
        assert residual > 0 or residual < 0
        return residual < 0 if increasing else residual > 0
    assert side(lo) and not side(hi)
    for _ in range(50):
        mid = (lo+hi)/2
        if side(mid):
            lo = mid
        else:
            hi = mid
    return lo, hi


def least_denominator(left, right):
    # left and right enclose the two endpoints. Reject ambiguous comparisons.
    bound = int(1/(right[0]-left[1]))+1
    for n in range(1, bound+1):
        first = int(n*left[0])+1
        last = int(n*right[1])
        if first <= last:
            assert first == last and left[1] < F(first,n) < right[0]
            return F(first,n)
    raise AssertionError("No rational in the tested denominator range")


def verify():
    ctx.prec = 180
    products = RATIO["product_polynomials"]()
    exponents = mass_exponents(products)
    spectra = RATIO["thresholds"]
    unit = spectra(products, arb(1))
    match = RATIO["MODEL"]["matched"]
    rows = []
    for n, height, eps in RATIO["high_zero_readings"]():
        scale = (2*arb.pi()*(height*height+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        low, v, _ = match(unit, eps, arb(1))
        anchor = v/low*scale
        def mass(r):
            return anchor*mass_ratio(r, exponents)
        for ref in (a("11/25"), a(2)):
            low, v, _ = match(spectra(products,ref), eps, arb(1), ref)
            assert (v/low*scale-mass(ref)).contains(0)
        low_branch = (F(".45"), F(".46"), False)
        high_branch = (F("2.04"), F("2.06"), True)
        lower_root = enclose_root(mass, TARGET, *low_branch)
        upper_root = enclose_root(mass, TARGET, *high_branch)
        left = enclose_root(mass, TARGET+HALF_UNIT, *low_branch)
        right = enclose_root(mass, TARGET-HALF_UNIT, *low_branch)
        left_high = enclose_root(mass, TARGET-HALF_UNIT, *high_branch)
        right_high = enclose_root(mass, TARGET+HALF_UNIT, *high_branch)
        def fit_interval(left, right, k):
            r = a(left[0]).union(a(right[1]))
            return match(spectra(products,r), eps, k, r)
        low, v, log_x = fit_interval(left, right, arb(1))
        assert low*low/arb.pi() < a("1/4") and log_x.exp()**2 < 4*arb.pi()
        assert v*eps*log_x.exp() > 1000
        _, _, minimum = fit_interval(left_high, right_high, 1/arb(2).sqrt())
        assert minimum.exp()**2 > 4*arb.pi()
        simplest = least_denominator(left, right)
        square_left = tuple(r*r for r in left)
        square_right = tuple(r*r for r in right)
        simplest_q = least_denominator(square_left, square_right)
        mass_witness = mass(a(simplest))
        assert abs(mass(a(simplest_q).sqrt()) - TARGET) < HALF_UNIT
        assert abs(mass_witness - TARGET) < HALF_UNIT
        print(n, "roots", float(sum(lower_root)/2), float(sum(upper_root)/2),
              "lowest denominator for r,q", simplest, simplest_q, "M/GeV", mass_witness.str(12))
        rows.append((left, right, simplest, simplest_q))
    assert all(early[0][1] < late[0][0] and early[1][1] < late[1][0]
               for early,late in zip(rows, rows[1:]))
    for size in (2,3,4,5,6):
        shared = []
        for start in range(len(rows)-size+1):
            left = max(row[0] for row in rows[start:start+size])
            right = min(row[1] for row in rows[start:start+size])
            assert left[1] < right[0] or right[1] < left[0]
            if left[1] < right[0]:
                shared.append((11485+start,11485+start+size-1))
        assert shared == ([(11485,11486),(11487,11488),(11489,11490)] if size==2 else [])
        print("Simultaneous rounded labels",size,shared)


if __name__ == "__main__":
    verify()

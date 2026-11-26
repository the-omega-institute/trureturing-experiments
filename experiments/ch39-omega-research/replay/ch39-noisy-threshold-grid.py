"""Finite timelike scan with bounded amplitude error in two isolated threshold bands."""
from fractions import Fraction as F
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

MATCH = run_path(str(Path(__file__).with_name("ch39-joint-moment-separation.py")))
a = MATCH["a"]


def mass_ranges(center, point=None, scale=None):
    bounds = [MATCH["WINDOW"], F(120, 10**6)]
    low = [-b for b in bounds]
    if point is not None:
        low, bounds = point, point
    q = center.root(*low).union(center.root(*bounds))
    if scale is None:
        scale = a(".9999").union(a("1.0001"))
    sigma, lam = [(1+a(l).union(a(b)))/k for l, b, k in zip(low, bounds, (256,12))]
    # Positive integer powers from the center matching identity.
    def coupling_squared(q, sigma, lam):
        factors = (q, (1+sigma*q)/2, (q+2)/3, (2*q+1)/3, sigma, 12*lam)
        return center.low**2*prod(f**p for f,p in zip(factors,(339,109,-1164,2,-111,-2))).root(42)
    g2 = coupling_squared(q, sigma, lam)
    k = a("8/5")
    z, w = g2*(1+3/k**2), g2*(3+k**2)
    ell = ((center.gamma/624)**2+a("1/4")).sqrt()
    real = ((ell+a("1/2"))/2).sqrt()
    d2, h2 = q*z*a("13/10"), q*z*a("7/2")
    fermions = [(d2*(ell+1-2*real),1), (d2*ell,-1), (h2*ell,-1), (h2*(ell+1+2*real),1)]
    sum_v, prod_v = (3+q)*z+w, (2+q)*z*w
    def roots(s, p):
        high = (s+(s*s-4*p).sqrt())/2
        return p/high, high
    vectors = list(zip(roots(sum_v,prod_v), (6,6)))+[(3*z+w,-6)]
    # Signed scalar indices already divided by three.
    scalars = [(4*lam,-16), ((q+2)*lam,-22), ((4*q+2)*lam,6)]
    for s, p, wt in ((q+7, 3*q+6, 6), (q+6, 3*q, -4), (4*q+6, 12*q, 14)):
        scalars += [(lam*m,wt) for m in roots(s,p)]
    scalars += [(lam*sigma*q,2),(lam*(1+sigma*q),-8),(arb(1),-2)]
    return [[(4*scale*m, w) for m, w in group] for group in (fermions, vectors, scalars)]


def models():
    ratio = MATCH["SOURCE"]["RATIO"]
    products = ratio["product_polynomials"]()
    powers = MATCH["SOURCE"]["PLANCK"]["mass_exponents"](products)
    for n, gamma, epsilon in ratio["high_zero_readings"]():
        yield n, MATCH["Center"](gamma, epsilon, products, powers, None)


BANDS = [(F("2.0569"), F("2.0602")), (F("10.782"), F("10.798"))]
NOISE = F("0.000050")


def cut_terms(groups):
    for group, (linear, cubic) in zip(groups, [(F(1,4),F(-1,12)),(F(-1,2),F(1,16)),(0,F(1,48))]):
        for t, w in group:
            yield t, a(w*linear), a(w*cubic)


def band_certificate(groups, which):
    target = groups[0][which][0]
    left, right = BANDS[0 if which == 0 else 1]
    assert a(left) < target < a(right)
    first = (right-left)/8
    # The first scan goes two steps past its right endpoint. Later brackets lie
    # at most two preceding steps from the true threshold on each side.
    domain = a(left-3*first).union(a(right+3*first))
    lipschitz = arb(0)
    for i, (t, linear, cubic) in enumerate(cut_terms(groups)):
        if i == which or t > domain:
            continue
        assert t < domain
        beta = (1-t/domain).sqrt()
        lipschitz += abs((linear/beta+3*cubic*beta)*t/(2*domain**2))
    bound = F("0.30") if which == 0 else F("0.125")
    assert lipschitz < a(bound)
    for j in range(4):
        h = a(first/F(4**j))
        beta = (h/(target+h)).sqrt()
        signal = beta/4-beta**3/12
        assert signal/2 > 2*a(bound)*h+4*a(NOISE)
    assert a(first/64)/target < a("3.42e-6")
    return bound


def locate(read, band, bound):
    left, right = band
    for _ in range(4):
        step = (right-left)/8
        previous = read(left)
        for i in range(10):
            x = left+i*step
            value = read(x+step)
            diff = value-previous
            cutoff = a(bound*step+2*NOISE)
            assert diff > cutoff or diff <= cutoff
            if diff > cutoff:
                left, right = x-step, x+step
                break
            previous = value
        else:
            raise AssertionError("No threshold found")
    return (left+right)/2


def demonstration(groups, slopes, parity):
    terms = list(cut_terms(groups))
    cache = {}
    def read(s):
        if s not in cache:
            signal = arb(0)
            for t, linear, cubic in terms:
                if t < a(s):
                    beta = (1-t/a(s)).sqrt()
                    signal += beta*(linear+cubic*beta**2)
                else:
                    assert t > a(s)
            sign = 1 if (len(cache)+parity)%2 == 0 else -1
            cache[s] = signal+sign*a(NOISE)
        return cache[s]
    estimates = [locate(read, band, slope) for band, slope in zip(BANDS, slopes)]
    for threshold, estimate in zip((groups[0][0][0], groups[0][3][0]), estimates):
        assert abs(a(estimate)/threshold-1) < a("3.42e-6")
    assert len(cache) <= 88
    return a(estimates[0]/estimates[1]), len(cache)


def verify():
    ctx.prec = 168
    candidates, trials = {}, []
    for n, center in models():
        groups = mass_ranges(center)
        slopes = [band_certificate(groups, which) for which in (0,3)]
        ell = ((center.gamma/624)**2+a("1/4")).sqrt()
        u = (2*ell+1).sqrt()
        ratio = a("13/35")*(ell+1-u)/(ell+1+u)
        tol = a("3.42e-6")
        candidates[n] = (ratio*(1-tol)/(1+tol)).union(ratio*(1+tol)/(1-tol))
        for i, delta in enumerate((-MATCH["WINDOW"], MATCH["WINDOW"])):
            for j, eta in enumerate((F(-120,10**6), F(120,10**6))):
                for sign in (-1, 1):
                    scale = 1+a(sign)/10000
                    spectrum = mass_ranges(center, (delta, eta), scale)
                    estimate, calls = demonstration(spectrum, slopes, i+j+(sign+1)//2)
                    trials.append((n, estimate, calls))
    for n, estimate, calls in trials:
        assert [m for m, interval in candidates.items() if interval.overlaps(estimate)] == [n]
    assert all(candidates[n] < candidates[n+1] for n in range(11485,11490))
    print("Boundary trials", len(trials), "largest number of distinct queries",
          max(calls for _,_,calls in trials), "certified maximum", 88)


if __name__ == "__main__":
    verify()

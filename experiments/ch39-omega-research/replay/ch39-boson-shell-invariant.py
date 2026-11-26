"""A rational bosonic readout of the matched shell and Planck normalization."""
from copy import copy
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

GRID = run_path(str(Path(__file__).with_name("ch39-noisy-threshold-grid.py")))
a = GRID["a"]
TARGET = GRID["MATCH"]["SOURCE"]["PLANCK"]["TARGET"]


def invariant(q, t, gauge):
    # t=sigma*q and gauge=g_low^2/(12*lambda).
    return (gauge**2*q**890*((1+t)/2)**219*((2*q+1)/3)**2
            /(((q+2)/3)**2264*t**221))


def reconstruct(two, small, weak, heavy):
    k = a("8/5")
    unit = two/2
    x = small/unit
    q = x*(6-x)/(3-x)
    t = weak/unit
    gauge = heavy/(12*(3+k*k)*(1+3/k**2)*unit)
    return invariant(q, t, gauge)


def derivative_bounds(x, t):
    q = x*(6-x)/(3-x)
    # Derivatives of log I with respect to log x and log t at fixed gauge.
    dq = 1-x/(6-x)+x/(3-x)
    dq_factor = 890-2264*q/(q+2)+4*q/(2*q+1)
    return dq*dq_factor, 219*t/(1+t)-221


def verify():
    ctx.prec = 220
    from sympy import Rational as R, cancel, symbols
    qsymbol, ssymbol = symbols("q sigma", positive=True)
    t = ssymbol*qsymbol
    # Eliminate the variable common potential coefficient from the two 42nd powers.
    B = qsymbol**339*((1+t)/2)**109*((qsymbol+2)/3)**(-1164)*((2*qsymbol+1)/3)**2*ssymbol**(-111)
    mass = qsymbol**(-327)*((1+t)/2)**(-107)*((qsymbol+2)/3)**1107*((2*qsymbol+1)/3)**(-1)*ssymbol**108
    s = invariant(qsymbol,t,1)
    assert cancel(s**21*B*mass**44) == 1
    assert R(88,42)*21 == 44
    amplitudes, centers = [], []
    windows = [(a(".665"),a(".668")), (a(".00430"),a(".00436")),
               (a(".0000337"),a(".0000341")), (a("26.3"),a("26.6"))]
    for n, center in GRID["models"]():
        centers.append((n,center))
        _, vectors, scalars = GRID["mass_ranges"](center)
        two = 4*a(".9999").union(a("1.0001"))*2*(1+a("-.000120").union(a(".000120")))/12
        pool = [two]+[v[0] for v in scalars[:-1]]+[v[0] for v in vectors]
        small, weak, heavy = scalars[5][0], scalars[-3][0], vectors[-1][0]
        for (lo,hi), node in zip(windows,(two,small,weak,heavy)):
            assert lo < node < hi
            assert sum(not (t < lo or hi < t) for t in pool) == 1
        expected = center.low**4*(center.anchor/TARGET)**88
        amplitudes.append(expected)
        for delta, eta in [(F(0),F(0)),(GRID["MATCH"]["WINDOW"],F(-120,10**6))]:
            cal = a("1.0001")
            _, vg, sg = GRID["mass_ranges"](center,(delta,eta),cal)
            reference = 4*cal*2*(1+a(eta))/12
            value = reconstruct(reference,sg[5][0],sg[-3][0],vg[-1][0])
            assert value.overlaps(expected) and value.rad() < abs(expected)*a("1e-9")
        # Give each of the four locations an independent relative error.
        e = a("1.0e-6")
        widen = ((1-e)/(1+e)).union((1+e)/(1-e))
        x_box = small/(two/2)*widen
        t_box = weak/(two/2)*widen
        ax, at = derivative_bounds(x_box,t_box)
        # Ordered coordinates: small, weak, heavy and the shared reference two.
        log_norm = abs(ax)+abs(at)+2+abs(-ax-at-2)
        assert log_norm < 1740
    gaps = [v/u-1 for u,v in zip(amplitudes,amplitudes[1:])]
    assert all(g > a(".00366") for g in gaps)
    error = 1740*(-(1-e).log())
    assert all(u*error.exp() < v*(-error).exp() for u,v in zip(amplitudes,amplitudes[1:]))
    print("Adjacent relative invariant intervals", [g.str(12) for g in gaps])
    half = GRID["MATCH"]["SOURCE"]["PLANCK"]["HALF_UNIT"]
    pairs = []
    for i, (n, left) in enumerate(centers):
        for m, right in centers[i+1:]:
            gain = (right.low/left.low).root(22)
            target_ratio = right.anchor/left.anchor*gain
            if target_ratio > a(F(TARGET+half,TARGET-half)):
                continue
            ti, tj = TARGET/target_ratio.sqrt(), TARGET*target_ratio.sqrt()
            assert TARGET-half < ti < tj < TARGET+half
            potential_ratio = gain**42
            calibration_ratio = 1/potential_ratio
            assert abs(potential_ratio-1) < a("1e-12")
            assert abs(calibration_ratio-1) < a("1e-12")
            ci = copy(left)
            ci.target42 = (ti/ci.anchor)**42
            q = ci.root(0, refine=True)
            assert a(".16")**2 < q < a(".17")**2
            sigma = a("1/256")
            factors = (q**-327*((1+sigma*q)/2)**-107*((q+2)/3)**1107
                       *((2*q+1)/3)**-1*sigma**108)
            assert (right.anchor**42*potential_ratio*factors).overlaps(tj**42)
            assert (right.low**2*potential_ratio**(-a("1/21"))*calibration_ratio).overlaps(left.low**2)
            pairs.append((n,m))
    assert pairs == [(11485,11486),(11487,11488),(11489,11490)]
    print("Rounded-scale pairs with identical second-current boson data",pairs)


if __name__ == "__main__":
    verify()

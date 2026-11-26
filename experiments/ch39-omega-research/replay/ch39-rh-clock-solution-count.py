"""Thresholds at which the number of finite one-loop matchings changes."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-rh-fixed-clock-crossings.py")))
scales = fixed["scales"]
H = fixed["HIGGS2"]


def at_unit(data, unit):
    poles, rho, p, modes, k, f = data
    return [c*(arb(unit)/fixed["CLOCK"])**2 for c in poles], rho, p, modes, k, f


def at_first_pole(data):
    poles, _, _, _, k, _ = data
    a, b = arb(8).log()/2, arb(24).log()/2
    return (k+12*a+4*b-44*(poles[0].sqrt()/1000).log())/16


def light_roots(w, data):
    poles, rho, p, modes, k, _ = data
    _, mf, _, _ = scales(w, k)
    roots = []
    for lam, n in modes:
        load = lam*mf**2
        def secular(x):
            return x-H-load*(p+x*sum(r/(x-c) for r, c in zip(rho, poles)))
        lo, hi = arb(H), H+p*load
        assert hi < poles[0] and secular(lo) < 0 and secular(hi) > 0
        for _ in range(80):
            mid = (lo+hi)/2
            value = secular(mid)
            if value < 0:
                lo = mid
            elif value > 0:
                hi = mid
            else:
                break
        roots.append((lo.union(hi), n))
    return roots


def minimum_residual(data):
    *_, k, f = data
    star = at_first_pole(data)
    q, _, u, v = scales(star, k)
    e = sum(n*(q/x.sqrt()).log() for x, n in light_roots(star, data))
    return 22*e-402*u-46*v-15*f


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    low = at_unit(data, 200_000_000_000)
    high = at_unit(data, 2_750_000_000_000)
    cmin, cmax = low[0][0], high[0][0]
    _, rho, p, modes, k, f = data
    a, b = arb(8).log()/2, arb(24).log()/2
    smin, smax = at_first_pole(high), at_first_pole(low)
    end = arb("9.1")
    assert a < b < smin < smax < end
    assert p-H*sum(r/(c-H) for r, c in zip(rho, low[0])) > 0
    assert H < cmin and 0 < p and 6*p < 1-arb(1)/64
    # Uniform bound at and after the turning point, including every physical
    # w above 7 for clock units above the entire finite parameter rectangle.
    ratio_m = (-2*smin).exp()
    ratio_x = H/cmin+6*p*ratio_m
    moment = sum(r*cmin/c for r, c in zip(rho, low[0]))
    assert ratio_x < 1
    elasticity = (1-arb(1)/64)/(1+6*ratio_m*moment/(1-ratio_x)**2)
    slope = 22*23*((1+arb(4)/11)*elasticity-arb(4)/11)-448
    assert slope > 0
    print("Uniform late residual slope exceeds", slope.lower().str(10))
    # Common finite upper witness and mass-floor test, even at the smallest
    # unit. For unbounded units the static upper mass bound is still valid.
    q, mf, u, v = scales(end, k)
    assert min(x for x, n in light_roots(end, low)) > 1000**2
    assert mf*arb(8).sqrt() > 1000 and q > 1000
    floor = sum(n*(q/(H+lam*p*mf**2).sqrt()).log() for lam, n in modes)
    assert 22*floor-402*u-46*v-15*f > 0
    high_at_a = fixed["endpoint"](a, high)[-1]
    high_at_seven = fixed["endpoint"](arb(7), high)[-1]
    assert high_at_a < 0 and high_at_seven < 0 and smin < 7
    q7, mf7, _, _ = scales(arb(7), k)
    assert q7**2 < cmax and mf7**2/cmax < ratio_m
    # The earlier absolute-clock proof supplies the two plateau thresholds.
    # A lower bound on the boundary t=1/M_F^2=0 is enough at the small unit.
    q0 = scales(a, k)[0]
    pure = 25*sum(fixed["positive"]((q0/c.sqrt()).log()) for c in low[0])
    assert 22*pure-15*f > 0
    for unit, sign in [(275_782_218_980, 1), (275_782_218_981, -1)]:
        g = minimum_residual(at_unit(data, unit))
        assert sign*g > 0
        print("First reference threshold, unit =", unit, "residual", g.str(14))
    print("Critical units ordered at 2.7578221898e11, 1.9612014e12, 2.4033935e12 GeV.")
    print("Successive counts including critical points: 0, 1, 2, 2, 3, 2, 1.")


if __name__ == "__main__":
    main()

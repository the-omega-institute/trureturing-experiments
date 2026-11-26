"""Calibrate a dimensionless shell reading to the first fermion threshold."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-rh-fixed-clock-crossings.py")))
positive = fixed["positive"]


def response(argument, data):
    c, rho, p, *_ = data
    z = [fixed["CLOCK"]**2/v for v in c]
    return sum(r/((1-p)*(1+argument*v)) for v, r in zip(z, rho))


def argument_for(value, data):
    # The part above 1,000,000 is irrelevant for the two book arcs.
    lo, hi = arb(0), arb(1_000_000)
    assert response(lo, data) > value > response(hi, data)
    for _ in range(100):
        mid = (lo+hi)/2
        excess = response(mid, data)-value
        if excess > 0:
            lo = mid
        elif excess < 0:
            hi = mid
        else:
            break
    return lo.union(hi)


def at_ratio(w, y, data):
    poles, rho, p, modes, k, f = data
    q, mf, u, v = fixed["scales"](w, k)
    c = [t/(fixed["CLOCK"]*y)**2 for t in poles]
    spectrum = fixed["shell"]["scalar_spectrum"](mf, c, rho, p, modes)
    masses = [(mf*x.sqrt(), n) for x, n in spectrum]
    total = sum(n*positive((q/m).log()) for m, n in masses)
    return 22*total-402*u-46*v-15*f, q, mf, masses


def universal_checks(data):
    _, _, p, modes, k, f = data
    w = arb("8.20")
    q, mf, u, v = fixed["scales"](w, k)
    limit = 402*u+46*v+15*f
    assert 22*sum(n*(q/(125**2+lam*p*mf**2).sqrt()).log() for lam, n in modes) > limit
    assert q > 3000 and 6/(3-arb(5).sqrt()) < 9
    assert 22*23*(q/3000).log() > limit and 6*p < 1-arb(1)/64
    assert 22*23*(1-arb(15)/(11*64)) > 448
    assert f > 0


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    universal_checks(data)
    records = []
    for numerator, denominator, left_w, right_w in ((33,34,"2.51009","2.51012"), (14,17,"1.45912","1.45916")):
        target = arb(numerator)/denominator
        t = argument_for(target, data)
        y = (t/8).sqrt()
        records.append(t)
        print("Arc", numerator, "/", denominator, "argument", t.str(14), "M_F/mu", y.str(14),
              "second fermion response", response(3*t, data).str(14), flush=True)
        lo, hi = arb(left_w), arb(right_w)
        left = at_ratio(lo, y, data)
        right = at_ratio(hi, y, data)
        assert left[0] < 0 < right[0] and right[2]*arb(8).sqrt() > 1000
        active = 0
        for (heavy, n), (light, nn) in zip(left[3],right[3]):
            assert light > 1000 and n == nn
            if heavy < right[1]:
                active += n
            else:
                assert left[1] < light
        assert active == 71
        mf = left[2].union(right[2])
        print("w", lo.union(hi).str(9), "Q", left[1].union(right[1]).str(9),
              "M_F", mf.str(9), "mu", (mf/y).str(9), "open", active)
    print("Squared threshold ratio needed to recover both book arcs:", (records[1]/records[0]).str(14))


if __name__ == "__main__":
    main()

"""Test the same zero-label ground state against the 128-letter gauge and Yukawa source."""
from itertools import islice
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
clock = run_path(str(root/"ch39-static-pulse-clock.py"))
gauge = run_path(str(root/"ch39-schur-winner-gauge-yukawa.py"))
pulses = clock["pulses"]
model = gauge["shortest"]


def verify():
    ctx.prec = 192
    offset = model["q"](model["model"]["fixed"]["OFFSET"])
    window = {}
    for n, zero in enumerate(acb.zeta_zeros(10815, 1278), 10815):
        series = pulses["pulse_records"](pulses["generator"](zero.imag))
        score = next(islice(series, 6, 7))[0]
        window[n] = score, 1/(arb(1)/4+zero.imag**2), 1/(offset+zero.imag**2)
    p0,w0,h0 = window.pop(11431)
    total = w0+sum(row[1] for row in window.values())
    others = [(p0-p,w/w0,h/h0-1) for p,w,h in window.values()]
    gap = min(d.lower() for d,w,h in others)

    def geometry(t):
        energy = clock["free_ground"](p0-t)
        a, da = clock["endpoint"](energy, 1)
        b, db = clock["endpoint"](energy, 3)
        return a-b,(da-db)/(db+p0*(da-db))

    def shift(t):
        _, q = geometry(t)
        norm, delta = arb(1),arb(0)
        for d,w,h in others:
            ratio = w*(1-d*q)*(t/(d-t))**2
            norm += ratio
            delta += ratio*h
        return delta/norm

    def secular(t, inverse):
        span, _ = geometry(t)
        return 4*span*(inverse-arb(1)/2)*total/w0*t-1+sum(w*t/(d-t) for d,w,h in others)

    def isolate(inverse):
        lo, hi = arb(0), (gap*arb(".99")).mid()
        assert secular(hi,inverse) > 0
        for _ in range(30):
            mid = ((lo+hi)/2).mid()
            sign = secular(mid, inverse)
            if sign < 0:
                lo = mid
            else:
                assert sign > 0
                hi = mid
        return lo.union(hi)

    def check_range(lo, hi, condition):
        todo = [(lo,hi)]
        while todo:
            a,b = todo.pop()
            if condition(shift(a.union(b))):
                continue
            assert b-a > arb("1e-14")
            mid = ((a+b)/2).mid()
            todo.extend(((a,mid),(mid,b)))

    end = isolate(0).upper()
    check_range(arb(0),isolate(1024).upper(),lambda h: abs(h) < arb("1.5e-6"))
    check_range(isolate(32).lower(),end,lambda h: h < arb("-.000085"))

    time = arb("26.81676714125598")
    state = gauge["gauge_state"](gauge["word"],time,(18,6,3))
    target = model["model"]["target"](11431,time)
    ratios = model["ratios"]
    values = [1/(r*u) for r,u in zip(ratios,state)]
    uncertainty = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    residual = [abs(v-target)+u for v,u in zip(values,uncertainty)]
    assert all(x < arb("1.005e-7") for x in residual)
    # Check actual speed bounds for arbitrary quark mixing throughout the scale window.
    flow = model["model"]
    coeff, loop = flow["running"], flow["loop"]
    initial = gauge["rotor_inputs"](gauge["word"],(18,6,3))
    def rates(step, squared_gauges, traces):
        speed = [(-2/(r*loop))*(beta+sum(b*x-d*y for b,x,d,y in
                 zip(row,squared_gauges,yukawa,traces))/loop)
                 for r,beta,row,yukawa in zip(ratios,coeff["BETA"],coeff["B"],coeff["D"])]
        assert all(abs(v) < arb(".13") for v in speed)
        assert speed[2]-speed[0] > arb(".13")
    flow["evolve"]([1/x for x in initial[:3]], [v.log() for v in initial[3:]],
                   arb(33),8192,diagonal_quarks=False,observe=rates)

    tolerance = arb(2)**-21
    delta_time = -(-arb("1.5e-6")).log1p()/2
    assert all(e+arb(".13")*delta_time < tolerance for e in residual)
    time_window = (2*tolerance+residual[0]+residual[2])/arb(".13")
    assert time_window < arb(".000009")
    strong_gap = target*(1-(1-arb(".000085"))*(2*time_window).exp())-residual[0]-arb(".13")*time_window
    assert strong_gap > arb(".00021") > tolerance
    print("weak reference coupling <= 1/1024: all three residuals <",
          max(e+arb(".13")*delta_time for e in residual).upper().str(8))
    print("strong reference coupling >= 1/32: any possible common scale must be within",
          time_window.upper().str(8), "of the anchor; first gauge mismatch >",strong_gap.lower().str(8))
    print("unit coupling unscaled target relative shift",shift(isolate(1)).str(10))


if __name__ == "__main__":
    verify()

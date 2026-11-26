"""Interval calculation of the joint shell-phase, gauge and bare gravity readouts."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import fraction


def verify():
    ctx.prec = 180
    source = Path(__file__).parent
    graph = run_path(str(source/"ch39-weighted-arc-zero-phase.py"))
    a, inputs, fit = run_path(str(source/"ch39-shell-zero-mass-hierarchy.py"))["parameters"]()
    num, den = fraction(graph["response"](3))
    target = a("1.2209e19")
    geom = 3*arb(5).sqrt()/2
    pi = arb.pi()

    def prediction(gamma, data):
        u = 1-1/(2*(gamma**2+a("1/4")))
        evaluate = graph["evaluated"]
        epsilon = evaluate(num, u)/evaluate(den, u)
        low, v, log_x, g5, g10 = fit(epsilon, a(1), data)
        mp = v/low*((2*pi/geom)*(gamma**2+a("1/4")+30+a("2/15"))).sqrt()
        return mp, g5, log_x.exp(), v*log_x.exp()*epsilon

    def for_label(j):
        z = acb.zeta_zero(j)
        assert z.real == a("1/2")
        return z.imag

    lo, hi = a("1.22085e19"), a("1.22095e19")
    m2, m3, m4 = [prediction(for_label(j), inputs)[0] for j in (7002, 7003, 7004)]
    assert m2 < lo < m3 < target < hi < m4 and target-m3 < m4-target
    assert -a(".00003933") < m3/target-1 < -a(".00003931")
    assert -a(".00006651") < m2/target-1 < -a(".00006648")
    assert a(".00007468") < m4/target-1 < a(".00007471")
    left, right = a("7267.30833"), a("7267.30835")
    assert prediction(left, inputs)[0] < target < prediction(right, inputs)[0]
    for j in (7003, *range(11485, 11491)):
        gamma = for_label(j)
        mp, *_ = prediction(gamma, inputs)
        delta = (target/mp).log()*7/(216*pi)
        shifted = [inputs[0], inputs[1]+delta, inputs[2]]
        ratio = (inputs[1]/shifted[1]).sqrt()
        tuned, g, x, light = prediction(gamma, shifted)
        assert (tuned-target).contains(0) and g*g < pi and x*x < 4*pi and light > 1000
        if j == 7003:
            assert abs(ratio-1) < a("7e-9")
        else:
            assert a("1.870e19") < mp < a("1.873e19")
            assert a(".000070") < ratio-1 < a(".000073")
        print(j, "bare Planck / GeV", mp.str(12), "relative shift in g2", (ratio-1).str(9))


if __name__ == "__main__":
    verify()

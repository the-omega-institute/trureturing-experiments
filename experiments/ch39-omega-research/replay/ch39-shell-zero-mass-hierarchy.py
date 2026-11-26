"""Certified zero selection and shortest-arc count in the one-loop matching model."""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def parameters():
    previous = run_path(str(Path(__file__).with_name("ch39-locking-coupling-match.py")))
    a = previous["a"]
    _, idx = run_path(str(Path(__file__).with_name("ch39-locking-scalar-spectrum.py")))["spectrum"]()
    scalar = [-sum(a(v[i])*arb(abs(f[-1])).log() for f, v in idx.items() if f[-1])/12
              for i in range(3)]
    h = (0, 2, 3)
    pi = arb.pi()
    initial = [4*pi/a(g)**2 for g in (".467453", ".63811", "1.0583")]

    def fit(epsilon, k, inputs=None):
        inv = initial if inputs is None else inputs
        d = [2*pi*ai-si+a(35-hi)/6 for ai, si, hi in zip(inv, scalar, h)]
        t0 = (5*d[0]+3*d[1]-8*d[2])/67
        r0 = a("15/16")*(d[0]-d[1]-a("109/15")*t0)
        z = a("134/7")*(r0/3-epsilon.log())
        low = z.exp()/arb(2).sqrt()
        t = t0-a("63/67")*z
        g10 = low*(3+k*k).sqrt()
        g5 = g10/k
        heavy = low*(3+k*k)/k
        vector = a("7/2")*(5*heavy.log()+3*z+25*g10.log())-a("33/6")
        log_x = (8*pi**2/low**2 - a("19/6")*t - 6*epsilon.log()
                 - 2*pi*inv[1]+scalar[1]+vector)/6
        return low, 1000*t.exp(), log_x, g5, g10

    return a, initial, fit


def shell_event():
    shell = [p for p in product(range(30), repeat=3) if sum(x*x for x in p) == 867]
    arcs = [(sum((x-y)**2 for x, y in zip(p, q)), p, q) for p, q in combinations(shell, 2)]
    smallest = min(d for d, p, q in arcs)
    assert len(shell) == 19 and smallest == 12
    words = {(a, b, a, b) for d, p, q in arcs if d == smallest for a, b in ((p, q), (q, p))}
    assert len(words) == 6
    return arb(len(words))/len(shell)**4


def verify():
    ctx.prec = 160
    a, inv, fit = parameters()
    pi = arb.pi()
    def strength(j):
        zero = acb.zeta_zero(j)
        assert zero.real == a("1/2")
        return 1/(zero.imag**2+a("1/4"))
    # For fixed epsilon the minimum of log X over k occurs at k^2=1/2.
    low51, _, min_log_x, _, _ = fit(strength(51), 1/arb(2).sqrt())
    assert low51 < a(".442") and min_log_x > 26
    slope_bound = -16*pi**2/a(".442")**2 + a("19/6")*63/67 + a("42/134") + a("231/2")
    assert slope_bound < 0
    assert fit(strength(55), a(1))[0] > pi.sqrt()
    for j in (52, 53, 54):
        epsilon = strength(j)
        low, v, log_x, g5, g10 = fit(epsilon, a(1))
        x = log_x.exp()
        assert g5.contains(g10) or (g5-g10).contains(0)
        alpha = low**2/pi
        assert (alpha < a("1/4")) if j == 52 else (alpha > a("1/4"))
        assert alpha < 1 and x*x < 4*pi and v*x*epsilon > 1000
        print(j, "g_low", low.str(10), "v/GeV", v.str(10), "x", x.str(9),
              "alpha_min", alpha.str(9))
    for j in (51, 53):
        e = strength(j)
        low = fit(e, a(1))[0]
        delta = -21/(333*pi)*(a(".64")/low).log()
        g2_ratio = (inv[1]/(inv[1]+delta)).sqrt()
        shifted = [inv[0], inv[1]+delta, inv[2]]
        result, v, log_x, g5, g10 = fit(e, a(1), shifted)
        x = log_x.exp()
        assert (result-a(".64")).contains(0) and abs(g2_ratio-1) < a(".00023")
        assert x**2 < 4*pi and v*x*e > 1000
    e = shell_event()
    assert a("4.3e-5") < e < a("4.9e-5") < a("1/6859")
    assert 5 < a("4.3e-5")*19**4 < 6 < a("4.9e-5")*19**4 < 7
    low, v, log_x, g5, g10 = fit(e, a(1))
    assert a(".63") < low < a(".64") and log_x.exp() < 1 and v*log_x.exp()*e > 1000
    print("Shortest-arc word", e.str(12), "v/GeV", v.str(12))


if __name__ == "__main__":
    verify()

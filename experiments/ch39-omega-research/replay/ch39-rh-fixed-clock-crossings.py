"""Count all physical one-loop matchings at a fixed arithmetic clock unit."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

shell = run_path(str(Path(__file__).with_name("ch39-kinetic-rh-gaussian-crossing.py")))
CLOCK = 2_100_000_000_000
HIGGS2 = 125**2


def positive(value):
    return value.max(arb(0))


def inputs():
    nodes, weights, _ = shell["quad"](21)
    a, b, c = [shell["base"]["gram"](i, j) for i, j in ((1, 0), (1, 2), (0, 2))]
    p = (1 + 2*a*b*c - a*a - b*b - c*c) / (1-c*c)
    poles = [CLOCK**2 / z for z in nodes]
    rho = [(1-p)*v for v in weights]
    modes = [(arb(6), 7), (arb(4), 6)]
    for j, n in ((5, 3), (3, 2)):
        modes.extend(((3-arb(j).sqrt(), n), (3+arb(j).sqrt(), n)))
    d1, d2, d3 = [8*arb.pi()**2/arb(g)**2 for g in (".467453", ".63811", "1.0583")]
    k = 5*d1-3*d2-2*d3
    f = arb(109)/15*(d2-d3)-arb(23)/6*(d1-d2)
    return poles, rho, p, modes, k, f


def scales(w, k):
    u = positive(w-arb(8).log()/2)
    v = positive(w-arb(24).log()/2)
    q = 1000*((k-12*u-4*v)/44).exp()
    return q, q*(-w).exp(), u, v


def endpoint(w, data):
    poles, rho, p, modes, k, f = data
    q, mf, u, v = scales(w, k)
    spectrum = shell["scalar_spectrum"](mf, [c/mf**2 for c in poles], rho, p, modes)
    masses = [(mf*x.sqrt(), n) for x, n in spectrum]
    total = sum(n*positive((q/m).log()) for m, n in masses)
    return q, mf, masses, 22*total-402*u-46*v-15*f


def monotonicity(data):
    poles, rho, p, _, k, _ = data
    c1 = poles[0]
    a, b = arb(8).log()/2, arb(24).log()/2
    # The unique scale w_star at which Q crosses the lowest pure pole.
    star = (k+12*a+4*b-44*(c1.sqrt()/1000).log())/16
    _, mf_star, _, _ = scales(star, k)
    assert 0 < a and a < b and b < star
    assert min(rho) > 0 and 0 < p and p < 1 and HIGGS2 < c1
    assert p-HIGGS2*sum(r/(c-HIGGS2) for c, r in zip(poles, rho)) > 0
    # Before star, N >= 48 and the total logarithmic mass elasticity <= 23.
    for s, penalty in [(arb(3)/11, 402), (arb(4)/11, 448)]:
        assert 22*((1+s)*23-s*48)-penalty < 0
    # After star, the only possible active root in each block is the lightest.
    # Its static trial bound and the 1 TeV floor yield this uniform elasticity.
    ratio = (HIGGS2+6*p*mf_star**2)/c1
    assert ratio < 1 and 6*p < 1-arb(1)/64
    elastic = (1-arb(1)/64)/(1+6*mf_star**2/c1/(1-ratio)**2)
    slope = 22*23*((1+arb(4)/11)*elastic-arb(4)/11)-448
    assert slope > 0
    print("w at first pure threshold", star.str(14), "late slope exceeds", slope.lower().str(10))
    return a, star


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = inputs()
    a, star = monotonicity(data)
    windows = [(".31499", ".31502", 1, 119), ("1.64484", "1.64488", -1, 71),
               ("8.19470", "8.19480", 1, 23)]
    for low, high, sign, expected in windows:
        left, right = arb(low), arb(high)
        q0, m0, spectrum0, r0 = endpoint(left, data)
        q1, m1, spectrum1, r1 = endpoint(right, data)
        print("w", low, high, "residuals", r0.str(10), r1.str(10), flush=True)
        assert sign*r0 < 0 and sign*r1 > 0
        assert (right < a, a < left and right < star, star < left).count(True) == 1
        opened = 0
        for (heavy, n), (light, multiplicity) in zip(spectrum0, spectrum1):
            assert n == multiplicity and light > 1000
            if heavy < q1:
                opened += n
            else:
                assert light > q0
        assert opened == expected and m1*arb(8).sqrt() > 1000
        print("Q in", q1.str(13), q0.str(13), "M_F in", m1.str(13), m0.str(13),
              "open", opened, "lightest lower endpoint", min(m for m, n in spectrum1).str(10))
    assert arb(windows[0][1]) < a < arb(windows[1][0]) < arb(windows[1][1]) < star < arb(windows[2][0])


if __name__ == "__main__":
    main()

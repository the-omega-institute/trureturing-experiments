"""Interval one-loop matching with masses from the locking vacuum."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def a(x):
    value = F(x)
    return arb(value.numerator)/value.denominator


def dot(x, y):
    return sum(u*v for u, v in zip(x, y))


def verify():
    ctx.prec = 160
    masses = run_path(str(Path(__file__).with_name("ch39-locking-scalar-spectrum.py")))
    _, indices = masses["spectrum"]()
    scalar = [-sum(a(idx[i])*arb(abs(f[-1])).log() for f, idx in indices.items()
                   if f[-1])/12 for i in range(3)]
    h = (0, 2, 3)
    p = (a("4/5"), a(0), a("1/2"))
    beta = (a("41/10"), -a("19/6"), -a(7))
    w = (5, 3, -8)
    pi = arb.pi()
    inv = [4*pi/a(g)**2 for g in (".467453", ".63811", "1.0583")]

    def reduced(g, finite):
        vector = [-a("7/2")*(i*(1/(2*g)).log()
                              + (5-i)*(1/(arb(2).sqrt()*g*g)).log()
                              + 25*(1/g).log()) - a(finite*(35-i))/6 for i in h]
        d = [2*pi*alpha - s - v for alpha, s, v in zip(inv, scalar, vector)]
        t = dot(w, d)/67
        r = a("15/16")*(d[0]-d[1] - (beta[0]-beta[1])*t)
        residuals = [di-bi*t-2*pi*16*pi/g**2 - (a("4/3")*py-2)*r
                     for di, bi, py in zip(d, beta, p)]
        return residuals, t, (r/3).exp()

    cases = [(0, "1.2722518353", "1.2722518355", "7.0953e14", "7.0955e14",
              "4.5908e-5", "4.5910e-5"),
             (1, "1.2467580789", "1.2467580791", "7.5628e14", "7.5630e14",
              "4.6071e-5", "4.6073e-5")]
    for finite, lo, hi, vlo, vhi, ylo, yhi in cases:
        assert reduced(a(lo), finite)[0][1] < 0 < reduced(a(hi), finite)[0][1]
        box = a(lo).union(a(hi))
        residuals, t, epsilon = reduced(box, finite)
        v = 1000*t.exp()
        assert all(e.contains(0) for e in residuals)
        assert a(vlo) < v < a(vhi) and a(ylo) < epsilon < a(yhi)
        print("Finite vector term", finite, "v [GeV]", v.str(12), "epsilon", epsilon.str(12))
        if finite:
            # epsilon(g_low)*(sqrt(2)*g_low)^(7/134) is constant.
            constant = epsilon*((box/arb(2).sqrt()).log()*a("7/134")).exp()
            low = a(".3").union(a("1.5"))
            general = constant*(-(arb(2).sqrt()*low).log()*a("7/134")).exp()
            assert a("4.3e-5") < general < a("4.9e-5")
            print("Allowed hierarchy for 0.3 <= g_low <= 1.5:", general.str(6))


if __name__ == "__main__":
    verify()

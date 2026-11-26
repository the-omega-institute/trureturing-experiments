"""SU(5) vectorlike adjoint masses and finite two-loop meeting intervals."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import arb, ctx

run = run_path(str(Path(__file__).with_name("ch39-two-loop-mirror-crossing.py")))
q = run["q"]


def doubled(slope, matrix):
    return [2*x for x in slope], [[2*x for x in row] for row in matrix]


def masses():
    light, heavy = q(3), q(11)
    x5 = (heavy-light)/5
    m5 = light+3*x5
    assert (m5+2*x5-heavy).contains(0) and (m5-3*x5-light).contains(0)
    a2 = (heavy**2-light**2)/25
    assert (light**2+25*a2-heavy**2).contains(0)
    # The in-plane displacement is imaginary in the 10. Real magnitudes then agree.


def verify():
    ctx.prec = 160
    masses()
    light = doubled(run["db"],run["dB"])
    heavy = doubled(list(map(q,(F(16,15),0,F(2,3)))),
                    run["sm"]["matrix"](["88/75 0 8/3","0 0 0","1/3 0 38/3"]))
    # The real-adjoint electron threshold can open slightly before the other singlets.
    lo, hi = F(57,4), F(29,2)
    width = F(1,65536)
    assert -(1-2*(q(2)-q(lo)).exp()).log() < q(width)
    for a,b in ((lo,lo),(hi,hi),(lo,hi)):
        t,diff,alpha = run["threshold_crossing"]([(F(2),F(2),*light),(a-width,b,*heavy)])
        if a == b == lo:
            assert diff > q(".008")
        elif a == b == hi:
            assert diff < q("-.037")
        else:
            assert q("5.27e14") < 1000*t.exp() < q("5.81e14")
            assert q("33.13") < alpha < q("33.27")
        print(a,b, "difference", diff.lower().str(10),diff.upper().str(10),
              "scale", (1000*t.exp()).lower().str(10),(1000*t.exp()).upper().str(10),
              "alpha",alpha.lower().str(10),alpha.upper().str(10))


if __name__ == "__main__":
    verify()

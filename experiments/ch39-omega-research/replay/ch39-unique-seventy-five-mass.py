"""单份 75 的全局匹配唯一性、无裸质量的输运障碍及中性子图谱隙。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

crossing = run_path(str(Path(__file__).with_name("ch39-unit-seventy-five-crossing.py")))
soft, smooth = crossing["soft"], crossing["smooth"]


def verify():
    ctx.prec = 256
    masses = crossing["shell_spectra"]()
    levels = {q: row if q else [arb(0), *row] for q, row in masses.items()}
    assert all(len(row) == 19 for row in levels.values())
    assert min(v.lower() for q in (5, 10) for v in levels[q]) > arb("1.9")
    geo = soft["geo"]
    radius2 = (geo["theta"](804)/geo["theta"](12))**2
    mz, em, weak = arb("91.1876"), arb("127.930"), arb(".23129")
    log_min = (arb("1e13")/mz).log()
    log_max = (arb("1e17")/mz).log()
    ref2 = radius2*(-2*log_min).exp()
    pi = arb.pi()
    c_base = 3*(1-weak)*em/5-3*weak*em+2/arb(".1180")
    c_low, c_high = (c_base+L/(5*pi) for L in (log_min, log_max))
    indices = crossing["index_table"]()
    coeff_c = {q: -y+3*w-2*c for q,(y,w,c) in indices.items()}
    coeff_d = {q: -y+c for q,(y,w,c) in indices.items()}
    assert coeff_c == {0:24, 5:-12, 10:-12}
    assert coeff_d == {0:12, 5:-3, 10:-9}

    def high(q, x):
        return sum(smooth["sheet"](v+x, radius2, False) for v in levels[q])

    def reference_bound(q, x):
        return sum(ref2/(8*(v+x)) for v in levels[q])

    # F(m; r, r exp(-L)) decreases with m, increases with L. At m=0
    # the difference extends to L. Use the positive endpoint for the lower bound.
    left, right = arb(15)/1000, arb(30)/1000
    small = (24*(high(0,left)-reference_bound(0,left))
             -12*(high(5,0)+high(10,0)))
    assert c_low + 2*small/(3*pi) > 0
    far = arb(128)
    assert c_high + 2*24*19*radius2/(3*pi*8*far) < 0

    # On a middle interval the positive neutral term is bounded at its
    # left end, while the two negative charged terms are bounded at its right.
    a = right
    boxes = 0
    while a < far:
        b = min(a*arb(17)/16, far)
        upper = 24*high(0,a)-12*sum(high(q,b)-reference_bound(q,b) for q in (5,10))
        assert c_high + 2*upper/(3*pi) < 0
        a = b
        boxes += 1

    # -dF/dx on the neutral zero mode: restrict the angular integral to y<=3/4.
    cutoff = arb(3)/4
    weight = (3*cutoff-cutoff**3)/2
    neutral_loss = weight*(1/(right+ref2/4)-1/(left+radius2*(1-cutoff**2)/4))/2
    charged_loss = 19/(2*arb("1.9"))
    assert neutral_loss > charged_loss
    dx_upper = {name: 2*coeff[0]*(charged_loss-neutral_loss)/(3*pi)
                for name, coeff in (("C", coeff_c), ("D", coeff_d))}
    assert dx_upper["C"] < -50 and dx_upper["D"] < -25

    # The reference derivatives in L are positive and at most ref2/(4*m).
    # Each charged or neutral tower has at most nineteen modes; bound both
    # residuals with absolute coefficients on the whole remaining rectangle.
    dref = 19*ref2/(4*left)
    c_error = 2*sum(abs(w) for w in coeff_c.values())*dref/(3*pi)
    d_error = (2*sum(abs(w) for w in coeff_d.values())/3+7*(3+3)/4)*dref/pi
    assert 1/(5*pi)-c_error > 0
    assert -arb(111)/(20*pi)+d_error < 0
    print("Mass outside (0.015,0.030) excluded, finite middle intervals", boxes)
    print("Upper mass slopes C,D", *(v.upper().str(8) for v in dx_upper.values()))

    # At x=0, the Dirac kernel dominates the vector kernel. Their cumulative
    # weight difference is 5*y*(1-y*y)/14. Compare beta coefficients in W.
    assert 7*2*3 <= 8*indices[0][1] and 7*3*3 <= 8*indices[5][1]
    slope = (arb(19)/12-arb(2)*indices[0][1]/3)/pi
    assert slope < 0
    weak_upper = weak*em + slope*log_min
    assert weak_upper < -arb("43.91")
    print("Zero bare mass, arbitrary shared transport: weak inverse below", weak_upper.upper().str(11))

    h = acb.zeta_zero(11431).imag
    leaf_weight = 1-arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    gap19 = 2-2*(pi/19).cos()
    assert leaf_weight*gap19 > arb(".02727739")
    assert arb(".02251149")/gap19 < arb(".825281")
    assert leaf_weight*(2-2*(pi/20).cos()) > arb(".02251149")
    assert (2-2*(pi/21).cos()) < arb(".02251147")
    print("Positive ordinary spectrum on shell subgraphs exceeds", (leaf_weight*gap19).lower().str(11))


if __name__ == "__main__":
    verify()

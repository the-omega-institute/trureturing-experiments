"""实际零点的标准算术权中，有多少能满足十九点壳的五荷三轻模条件。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_series, ctx


root = Path(__file__).parent
window = run_path(str(root/"ch39-common-step-zero-window.py"))
shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))


def time_data():
    s = arb_series([arb(1)/2,1],prec=3)
    xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    total = xi[2]/xi[0]
    s = arb_series([arb(3)/2,1],prec=2)
    xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    return total,xi[1]/(2*xi[0]*total)


def three_mode_test():
    levels = [(15-arb(189).sqrt())/2, (15+arb(189).sqrt())/2, (15+arb(261).sqrt())/2]
    eps = arb("3e-8")
    assert shell["soft"]["point"]["study"]["offset"] >= 0
    delta = arb(23)/25*arb(2)**-20/14**2
    assert 4*delta + 2*delta/3 < eps
    assert len(shell["geo"]["edges"]) == 135
    count_at = {q: [[window["supplier"](q,t+sign*eps) for t in levels] for sign in (-1,1)]
                for q in (6,3,4,2,1)}
    def accepted(height):
        assert height>14
        phase=window["geometry"]["phases"](height)
        for low,high in count_at.values():
            below=[f(phase) for f in low]
            above=[f(phase) for f in high]
            assert below==above
            if 1+below[0]+below[2]-below[1]!=3:
                return False
        return True
    return accepted


def verify(count=8192):
    ctx.prec,ctx.cap = 112,3
    three_modes=three_mode_test()
    total, first_time = time_data()
    zeros = acb.zeta_zeros(1,count)
    assert len(zeros)==count and zeros[0].imag > 14
    seen, accepted = arb(0),arb(0)
    hits = []
    for n,zero in enumerate(zeros,1):
        h = zero.imag
        score = 1/h**2
        seen += score
        if three_modes(h):
            hits.append(n)
            accepted += score
    assert hits[0] == 200
    if count == 8192:
        assert len(hits) == 67
    mass_low = accepted/total
    mass_hi = (accepted+total-seen)/total
    hmin = zeros[hits[0]-1].imag
    minimum_first = 1/(1+1/hmin**2)
    hill = 1/(1+1/zeros[0].imag**2)
    outside = (minimum_first-first_time)/(minimum_first-hill)
    scale = hmin*((1-first_time)/first_time).sqrt()
    if count == 8192:
        assert arb(".0009364") < mass_low < mass_hi < arb(".007719")
    assert minimum_first > arb(".99999363") and outside > arb(".32092") and scale > arb("15.88")
    print(count,len(hits),"首个候选的序号、高度、全谱权份额的下上界",
          hits[0],hmin.str(12),mass_low.str(11),mass_hi.str(11))
    print("未缩放首项下界、RH 首项、所需外部份额下界、逆高度所需共同因子下界",
          minimum_first.str(13),first_time.str(13),outside.str(10),scale.str(10))
    return [zeros[n-1].imag for n in hits]


if __name__ == "__main__":
    verify()

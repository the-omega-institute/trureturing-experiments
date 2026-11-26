"""由荷二曲率的薄本征值带给三窗零点的相对上密度界。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

window = run_path(str(Path(__file__).with_name("ch39-common-step-zero-window.py")))


def cover():
    low, high = arb("18.94085"), arb("18.94416")
    pi = arb.pi()
    limit = 14

    @cache
    def counts(depth):
        error = 4*pi/2**depth
        return [window["supplier"](2, t) for t in (low-error, high+error, low+error, high-error)]

    @cache
    def phase(depth,index):
        return acb(0,pi*(2*index+1)/2**(depth+1)).exp()

    pending = [(3,i,j) for i in range(8) for j in range(8)]
    inside, uncertain = 0, 0
    while pending:
        depth,i,j = pending.pop()
        clocks = {724:phase(depth,i), 804:phase(depth,j)}
        low_out, high_out, low_in, high_in = counts(depth)
        if low_out(clocks) >= 15 or high_out(clocks) < 15:
            continue
        if depth > 12 and low_in(clocks) < 15 and high_in(clocks) >= 15:
            inside += 4**(limit-depth)
        elif depth == limit:
            uncertain += 1
        else:
            pending.extend((depth+1,2*i+a,2*j+b) for a in range(2) for b in range(2))
    assert (inside,uncertain) == (90230,155790)
    return arb(inside+uncertain)/4**limit


def restrict_corners():
    # 先将偶数荷的环面按 π 折叠，再把每一只钟限制在离 0 至多 π/16 的区间。
    pi = arb.pi()

    @cache
    def center(depth, index):
        return acb(0, pi*(2*index+1)/2**(depth+1)).exp()

    @cache
    def band_counts(depth):
        error = 4*pi/2**depth
        return (window["supplier"](2, arb("18.94085")-error),
                window["supplier"](2, arb("18.94416")+error))

    pending = [(3, i, j) for i in range(8) for j in range(8)]
    accepted_area = arb(0)
    excluded_area = arb(0)
    while pending:
        depth, i, j = pending.pop()
        scale = 2**depth
        if all(16*min(s+1, scale-s) <= scale for s in (i, j)):
            accepted_area += arb(1)/scale**2
            continue
        clocks = {724: center(depth, i), 804: center(depth, j)}
        low, high = band_counts(depth)
        if low(clocks) >= 15 or high(clocks) < 15:
            excluded_area += arb(1)/scale**2
            continue
        assert depth < 11
        pending.extend((depth+1, 2*i+a, 2*j+b) for a in range(2) for b in range(2))
    assert accepted_area == arb(1)/64 and accepted_area + excluded_area == 1

    # 两钟在 π 附近时，荷一至少有四个轻模，按谱分界两侧的严格惯性计数。
    depth = 6
    error = 2*pi/2**depth + arb("3e-8")
    lower = (15-arb(189).sqrt())/2
    upper = (15+arb(261).sqrt())/2
    cuts = [window["supplier"](1, cut) for cut in (lower-error, upper-error, 15-lower+error)]
    for i in range(60, 68):
        for j in range(60, 68):
            clocks = {724: center(depth, i), 804: center(depth, j)}
            low_count, upper_count, middle_count = [count(clocks) for count in cuts]
            assert 1 + low_count + upper_count - middle_count >= 4


def eigen_band():
    unit = 32*arb("2.435e18")/724**2
    k = (9*arb(3).log()-12*arb(2).log())/2
    ratio = arb("246.22")**2/(unit*arb("1.22089e19"))
    fraction = 1-(ratio/9)**2/(arb(".03")*k)
    mass_lo = ratio*(fraction*k/(arb.pi()*arb("1.002e-122")))**(arb(1)/4)
    mass_hi = ratio*(k/(arb.pi()*arb(".998e-122")))**(arb(1)/4)
    delta = arb(23)/25*arb(2)**-20/14**2
    center_lo = 15-delta/3
    lower = (center_lo+(center_lo**2+4*mass_lo).sqrt())/2
    upper = (15+(225+4*mass_hi).sqrt())/2+4*delta
    return lower, upper


def verify():
    ctx.prec = 144
    lower, upper = eigen_band()
    assert lower > arb("18.94085") and upper < arb("18.944146")
    cap = cover()
    assert cap < arb(".00091650")
    restrict_corners()
    cap *= arb(3)/4
    assert cap < arb(".000687372")
    print("三代与曲率的必要本征值带、全部合格正零点的相对上密度界",
          lower.str(12),upper.str(12),cap.str(13))


if __name__ == "__main__":
    verify()

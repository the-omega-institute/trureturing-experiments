"""原壳最小割分界下的三窗角距、六面味坐标及零点上密度。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
sign = run_path(str(root / "ch39-three-window-phase-sign.py"))
inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))


def cell_count(charge, level):
    positive, coefficients = inertia["schur_coefficients"](charge, level)
    radius = charge*arb.pi()/16
    second = [abs(a)+abs(b)+4*abs(c) for _, a, b, c in coefficients]

    @cache
    def phase(n, i):
        return acb(0, radius*((2*arb(i)+1)/n-1)).exp()

    def counts(n, i, j):
        u, v = phase(n, i), phase(n, j)
        w = u.conjugate()*v
        half = radius/n
        low = high = 16-positive
        for (constant, a, b, c), second_order in zip(coefficients, second, strict=True):
            au, bv, cw = a*u, b*v, c*w
            value = constant+2*(au+bv+cw).real
            # 对两个相角分别求导。二阶余项用 |Δx|、|Δy|≤half 统一控制。
            error = 2*half*(abs((au-cw).imag)+abs((bv+cw).imag)) + half**2*second_order
            low += int(value < -error)
            high += int(not (value > error))
        return low, high

    return counts


def cover():
    # 仅用 T<10 的必要条件。软辐边对单位边低端 λ 的误差容许 3×10^-8。
    curvature = sign["density"]["eigen_band"]()
    lower, upper = arb("18.94088212"), arb("18.94414547")
    assert lower < curvature[0] < curvature[1] < upper
    entry = (15-arb(185).sqrt())/2
    epsilon = arb("3e-8")
    delta = arb(23)/25*arb(2)**-20/14**2
    c = 15-delta/3
    assert (c-(c*c-40).sqrt())/2+4*delta < entry+epsilon
    middle = ((15+arb(185).sqrt())/2-epsilon, (15+arb(265).sqrt())/2+epsilon)
    assert middle[0] < (c+(c*c-40).sqrt())/2 and (15+arb(265).sqrt())/2+4*delta < middle[1]
    for level, expected in zip(middle, (7, 9), strict=True):
        sign["certify_no_crossing"](sign["schur_coefficients"](6, level), 6, -1)
        assert sign["window"]["supplier"](6, level)({724: 1, 804: 1})+1 == expected

    six = cell_count(6, entry+epsilon)
    low, high = (cell_count(2, x) for x in (lower, upper))
    aperture = arb(".256255")
    resolution = 2048
    pending = [(1, 0, 0)]
    units = 0
    while pending:
        n, i, j = pending.pop()
        six_low, six_high = six(n, i, j)
        if six_high < 1:
            continue
        a, b = low(n, i, j)
        if a >= 16:
            continue
        e, f = high(n, i, j)
        if f < 16:
            continue
        separated = arb.pi()*(abs(i-j)+1)/(8*n) < aperture
        if separated and (n >= resolution or (b < 16 <= e and six_low >= 1)):
            # 在 resolution 层以后继续细分只是为了角距界；面积用精确的二进有理数。
            units += arb(1)/n**2
            continue
        assert n < 2**24, (n, i, j)
        pending.extend((2*n, 2*i+x, 2*j+y) for x in (0, 1) for y in (0, 1))
    # 两钟零相位方框的环面归一面积为 1/256。
    bound = units/256
    flavor = sign["physical"]["joint"]["flavor"]
    a, b = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    bottom = a*aperture.cos()-abs(b)*aperture.sin()
    print("角距和味坐标的单侧界、面积外界", aperture.str(12), bottom.str(16), bound.str(17))
    assert bottom > arb(".96619523") and bound < arb(".00012050834")


def open_witness():
    physical = sign["physical"]
    model, initial = physical["model_inputs"]()
    ctx.prec = 224
    joint = physical["joint"]
    chiral = joint["chiral"]
    supply = chiral["source"]["spectral_supplier"](model)
    angles = {724: arb(".1624062"), 804: arb("-.0938475")}
    radius = arb(2)**-46
    phases = {d: acb(0, t).exp() for d, t in angles.items()}
    rows = supply(arb(10)**13, clocks=phases)
    bands = {}
    for q in (0, *model["charges"]):
        values = chiral["ordered"]([v for v, charge in rows if q == charge])
        bands[q] = [v+arb(0, 320*q*radius) for v in values]
    # spectral_supplier 对全体高度 >10000 的三条软辐边使用统一误差。
    assert all(bands[q][2] < arb("9.99999999") and bands[q][3] > 10 for q in model["charges"])
    gauge = joint["gauge"]
    _, _, alpha = gauge["crossing_enclosure"](model, initial, bands, False, 724,
                                             "29.04", ".1184", recenter=True)
    vacuum = joint["vacuum"]
    heavy = 32*arb("2.435e18")/724**2*bands[2][15]
    squared = (arb("246.22")**2/(heavy*arb("1.22089e19")))**2
    paired = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](bands), min_gap="1e-18")
    fields = paired+vacuum["low_sector"](4, squared)
    assert vacuum["newton"](fields, arb(".03")) < 1 < vacuum["newton"](fields, arb(".033"))
    coeff = (9*arb(3).log()-12*arb(2).log())/2
    curve = squared**2*coeff/(arb.pi()*arb("1e-122"))*(1-arb(0).union(squared/(arb(".03")*coeff)))
    assert arb(".998") < curve < arb("1.002") and arb(".117") < alpha < arb(".119")
    box = {d: acb(0, t+arb(0, radius)).exp() for d, t in angles.items()}
    coordinate = joint["flavor"]["singlet_at"](box)[1]
    assert angles[724]-angles[804]-2*radius > arb(".25625369")
    assert coordinate < arb(".96619557")
    for n, lo, hi in ((8236, ".21780052", ".21780053"), (148341, ".23147506", ".23147507")):
        actual = chiral["source"]["geometry"]["phases"](acb.zeta_zero(n).imag)
        difference = abs((actual[724]/actual[804]).arg())
        assert arb(lo) < difference < arb(hi)
    print("开方窗：强耦合、曲率比和六面味坐标", alpha.str(10), curve.str(13), coordinate.str(14))


if __name__ == "__main__":
    ctx.prec = 160
    cover()
    open_witness()

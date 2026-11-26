"""不指定零点谱权时，三窗与精确平均跃迁概率所容许的最大 CP 亏额。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
sign = run_path(str(root / "ch39-three-window-phase-sign.py"))
cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
window = sign["window"]
aperture = ".215915"


def counts_on_cell(q, level):
    definite = sign["determinant_sign_on_cell"](sign["schur_coefficients"](q, level), q, -1)
    inertia = window["supplier"](q, level)
    radius = arb.pi()/16

    def count(n, i, j):
        if not definite(n, i, j):
            return None
        phases = {d: acb(0, radius*((2*arb(index)+1)/n-1)).exp()
                  for d, index in ((724, i), (804, j))}
        return 1 + inertia(phases)

    return count


def angular_bound():
    # 荷六在靠近 15 的放宽轻模谱段内恰有两根，三模要求最低本征值落入另一段。
    tolerance = arb("3e-8")
    edges = ((15+arb(189).sqrt())/2-tolerance, (15+arb(261).sqrt())/2+tolerance)
    for level, expected in zip(edges, (7, 9)):
        sign["certify_no_crossing"](sign["schur_coefficients"](6, level), 6, -1)
        assert 1 + window["supplier"](6, level)({724: 1, 804: 1}) == expected
    cutoff = (15-arb(189).sqrt())/2 + tolerance
    lower, upper = sign["density"]["eigen_band"]()
    assert arb("18.94085") < lower < upper < arb("18.944146")
    six = counts_on_cell(6, cutoff)
    two_low = counts_on_cell(2, arb("18.94085"))
    two_high = counts_on_cell(2, arb("18.944146"))
    pending = [(1, 0, 0)]
    while pending:
        n, i, j = pending.pop()
        if arb.pi()*(abs(i-j)+1)/(8*n) < arb(aperture):
            continue
        if six(n, i, j) == 0:
            continue
        low = two_low(n, i, j)
        if low is not None and low >= 16:
            continue
        high = two_high(n, i, j)
        if high is not None and high < 16:
            continue
        assert n < 524288
        pending.extend((2*n, 2*i+a, 2*j+b) for a in (0, 1) for b in (0, 1))


def open_windows(centers=((".15674163", "-.05916834"), (".072", ".072"))):
    physical = sign["physical"]
    model, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    supply = chiral["source"]["spectral_supplier"](model)
    ctx.prec = 224
    radius = arb(2)**-42
    coordinates = []
    for a, b in centers:
        angles = {724: arb(a), 804: arb(b)}
        phases = {d: acb(0, t).exp() for d, t in angles.items()}
        rows = supply(arb(10)**13, clocks=phases)
        masses = {}
        for q in (0, *model["charges"]):
            levels = chiral["ordered"]([v for v, charge in rows if charge == q])
            masses[q] = [v + arb(0, 320*q*radius) for v in levels]
        _, alpha, curvature = physical["layer_observables"](model, initial, masses)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curvature < arb("1.002")
        box = {d: acb(0, t + arb(0, radius)).exp() for d, t in angles.items()}
        coordinates.append(physical["joint"]["flavor"]["singlet_at"](box)[1])
        print("相位窗、强耦合与曲率比", a, b, alpha.str(9), curvature.str(12))
    return coordinates


def verify():
    ctx.prec = 160
    angular_bound()
    witnesses = open_windows()
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    a, b = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    lo = (a*arb(aperture).cos()-abs(b)).lower()
    hi = (a*a+b*b).sqrt().upper()
    assert lo > arb(".97561882")
    data = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())

    def point(coordinate):
        up = cp["two"]["probabilities"](coordinate, data)[0][0]
        current = cp["cp_value"](coordinate, polys)
        return up, current

    target, reference = point(flavor["singlet"](11661)[1])

    def chord_loss(left, right):
        up_left, j_left = point(left)
        up_right, j_right = point(right)
        weight = (target-up_right)/(up_left-up_right)
        assert 0 < weight < 1
        return 1 - (weight*j_left + (1-weight)*j_right)/reference

    pending = [(lo, hi)]
    while pending:
        left, right = pending.pop()
        c = left.union(right)
        up = cp["two"]["probabilities"](c, data)[0]
        current = cp["cp_value"](c, polys, jet=True)
        if up[1] < 0 and current[2]*up[1]-current[1]*up[2] > 0:
            continue
        assert right-left > arb("1e-9")
        middle = (left+right)/2
        pending.extend(((left, middle), (middle, right)))
    upper = chord_loss(lo, hi)
    actual = chord_loss(*witnesses)
    assert actual > arb(".1656788") and upper < arb(".1656886")
    print("味坐标下界、零点两窗可达亏额与全谱亏额上界", lo.str(12), actual.str(14), upper.str(14))


if __name__ == "__main__":
    verify()

"""六项时间、味驻点与逐高度带两项跃迁均值固定时，只在任意远尾重配 CP。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, arb_poly, ctx


def limiting_pair(windows):
    centers = [windows["inner_centers"][i] for i in (0, 1, 3, 5)]
    # 区间证书覆盖四个开窗。下列极限取各窗的确切中心。
    windows["previous"]["open_windows"](centers)
    jets = windows["observables"]()
    cp = windows["cp"]
    data = cp["two"]["spectral_data"]()
    columns = []
    probabilities = []
    for a, b in centers:
        phases = {724: acb(0, arb(a)).exp(), 804: acb(0, arb(b)).exp()}
        c = cp["flavor"]["singlet_at"](phases)[1]
        columns.append([1, *[f[0] for f in jets(c)]])
        probabilities.append([series[0] for series in cp["two"]["probabilities"](c, data)])
    matrix = arb_mat(columns)
    assert arb(".4894581") < matrix.det() < arb(".4894583")
    first = arb_mat([row[:3] for row in columns[:3]]).transpose()
    signed = first.solve(arb_mat([[-v] for v in columns[3][:3]])).entries()+[arb(1)]
    assert all(v > 0 or v < 0 for v in signed)
    mass = sum(w for w in signed if w > 0)
    halves = [[(sign*w/mass if sign*w > 0 else arb(0)) for w in signed] for sign in (1, -1)]
    assert all(sum(part).contains(1) for part in halves)
    means = [[sum(w*pair[j] for w, pair in zip(part, probabilities)) for j in (0, 1)]
             for part in halves]
    assert all(a.overlaps(b) for a, b in zip(*means, strict=True))
    losses = [sum(w*row[3] for w, row in zip(part, columns)) for part in halves]
    assert losses[0]-losses[1] > arb(".05015")
    return means[0], losses


def global_margin(rows, flavor_masses, invariants, ref):
    # t=c-c₀。归一行列式为 1+t*a(t)，线性项由精确驻点约束消去。
    # 对 x>1000 的行保留 log 的负二次余项，其余各行使用 log R≤R−1。
    t = arb(0, 1)-ref
    gap = arb(0)
    for row, mass in zip(rows, flavor_masses, strict=True):
        x = row[1]
        determinant = 1+sum(e*x**i for i, e in enumerate(invariants, 1))
        ratio = (determinant*(1/determinant(ref)))(arb_poly([ref, 1]))
        gap -= mass*arb_poly(ratio.coeffs()[2:])(t)
        if x > 1000:
            divided = arb_poly(ratio.coeffs()[1:])(t)
            assert divided > 0
            maximum = max(arb(1), (1+t*divided).upper())
            gap += mass*divided.lower()**2/(2*maximum**2)
    return gap


def verify(root):
    windows = run_path(str(root / "ch39-three-window-two-probability-cp.py"))
    pinned = run_path(str(root / "ch39-pinned-cp-zero-support.py"))
    (up, charm), tail_losses = limiting_pair(windows)
    joint = pinned["joint"]
    seeds = pinned["previous"]["indices"][:10]
    rows, times, invariants, ref, probabilities, current = joint["columns_at"](2_903_600_000, seeds)
    ctx.prec = 640
    columns, target = pinned["problem"](rows, times, probabilities)
    basis = arb_mat(list(zip(*columns, strict=True)))
    tail = arb(".1684")
    at_infinity = [1]*6+[0, up-probabilities[6], charm-probabilities[7]-arb("2e-12"), 1]
    rhs = arb_mat([[v-tail*z] for v, z in zip(target, at_infinity, strict=True)])
    weights = basis.solve(rhs).entries()
    assert all(v > arb("4.66e-11") for v in weights)
    assert (sum(weights)+tail).contains(1)
    flavor_masses = [w*r[2] for w, r in zip(weights, rows, strict=True)]
    total = sum(flavor_masses)+tail
    assert arb("1.00208167") < total < arb("1.00208169")
    # 无穷处 x=0，行列式恒为 1，对二次余项的贡献为零。
    gap = global_margin(rows, flavor_masses, invariants, ref)
    assert gap > arb("4.0e-13")
    seed_loss = sum(w*(1-r[4]/current) for w, r in zip(flavor_masses, rows, strict=True))
    high, low = [(seed_loss+tail*l)/total for l in tail_losses]
    assert high > arb(".07179808") and low < arb(".06336985")
    assert high-low > arb(".00842824")
    print("两支的极限 CP 相对亏额、差额与未归一作用的二次余量",
          *[v.str(14) for v in (low, high, high-low, gap)])


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)

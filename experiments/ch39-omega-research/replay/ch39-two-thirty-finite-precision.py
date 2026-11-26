"""第 2、30 枚真实零点共用非负有效截止的最优相对误差对数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb


def verify(root):
    fields = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    expdiff = run_path(str(root / "ch39-two-zero-cutoff-support.py"))["exponential_difference"]
    _, _, (first, other) = fields(root, (2, 30))
    for g, alpha, mass in (first[0], other[0]):
        assert g == -62 and (alpha-arb(1)/6).contains(0) and mass.is_zero()
    pi = arb.pi()
    k, w = 1/(12*pi), 31/(4*pi)
    eps = arb("1e-122")
    log10 = arb(10).log()
    eta = arb("1e-15")
    X = arb("1e-14")

    # 两套谱中最轻的重行在第 30 枚。以它的指数为公共因子。
    pairs = []
    for (g, alpha, z), (h, beta, mass) in zip(first[1:], other[1:], strict=True):
        assert h == g and (alpha-beta).contains(0)
        pairs.append((-g, z, mass))
    lead = min(range(len(pairs)), key=lambda i: float(pairs[i][2]))
    spin, _, ref = pairs[lead]
    remaining = [(abs(s), m) for s, z, h in pairs for m in (h, z)]
    del remaining[2*lead]
    gap = min((m-ref).lower() for _, m in remaining)
    assert gap > 0

    # 可行上界：限定为两点截止测度，剖面是红外单位高度和一段窄高台。
    r = 1-eta
    B, C = (w/k)*(1+r), w*r
    t = 2*eps/(B+(B*B-4*C*eps).sqrt())
    s = r*t
    height = (1/k-s)/(t-s)
    assert 0 < s < t < 1 and height > 1
    abs_alpha = max(sum(abs(a) for _, a, _ in rows[1:]).upper() for rows in (first, other))
    abs_spin = max(sum(abs(g) for g, _, _ in rows[1:]) for rows in (first, other))
    thermal_error = max((6*abs_alpha).upper(), (arb(abs_spin)/62).upper())
    log_floor = (ref.lower()/t-thermal_error.log())/log10

    # 第 2 枚的重场自旋前缀非负。正自旋延后到质量上端，负自旋提前到下端。
    events = sorted(((z.upper() if g < 0 else z.lower(), -g) for g, _, z in first[1:]))
    balance = 0
    for _, charge in events:
        balance += charge
        assert balance >= 0
    # Abel 求和给 W'_2 >= 2wx；高于 X 的牛顿负部可由实测真空读数控制。
    tau = eta/100
    negative_k = sum(abs(a) for _, a, _ in first if a < 0)/(2*pi)
    high_k_loss = negative_k*(1+tau)*eps/(2*w*X)
    assert high_k_loss < arb("1e-105")
    # 对 0<x<=X，两个导数同无质量导数的相对偏差均小于 tiny。
    tiny = arb("1e-200")
    assert thermal_error*(-ref/X).exp() < tiny
    a = (1-eta)*k*eps/(2*w)
    assert 0 < a < X
    tail_fraction = 1-tau - (2*w/k)*(a/eps)*(1+tiny)/(1-tiny)*(1+tau+high_k_loss)
    assert tail_fraction > arb("9e-16")

    # 证明 e^(ref/x) * (4π/x) * (W'30-W'2) 在整段均大于 margin。
    margin = arb("1e-8")
    assert spin-sum(weight for weight, _ in remaining)*(-gap/X).exp() > margin
    left = X
    leaves = 0
    while left < 1:
        right = min(2*left, arb(1))
        todo = [(left, right)]
        while todo:
            lo, hi = todo.pop()
            interval = lo.union(hi)
            diff = sum(spin*expdiff(interval, z, h, ref) for spin, z, h in pairs)
            if diff > margin:
                leaves += 1
            else:
                assert hi-lo > arb("1e-22"), (lo, hi)
                midpoint = (lo+hi)/2
                todo.extend(((lo, midpoint), (midpoint, hi)))
        left = right

    spin_total = sum(abs(g) for g, _, _ in first)
    prefactor = tail_fraction*margin/(2*spin_total)
    log_ceiling = (ref.upper()/a-prefactor.log())/log10
    assert -tau.log()/log10 < log_ceiling
    assert arb("5.88583332912177e112") < log_floor
    assert log_ceiling < arb("5.88583332912180e112")
    print("第 30 枚的最轻重质量平方比与两谱间隙", ref.str(20), gap.str(16))
    print("正剖面的高台端点、峰值与重场误差常数", t.str(17), height.str(15), thermal_error.str(12))
    print("红外以上的真空份额下界、差值覆盖段数、牛顿负部上界",
          tail_fraction.str(12), leaves, high_k_loss.str(12))
    print("log10(1/E) 的下界与上界", log_floor.str(21), log_ceiling.str(21))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)

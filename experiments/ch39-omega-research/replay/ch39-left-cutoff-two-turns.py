"""左邻荷六第一枢轴同时匹配三个读数所需的最少能序转向。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    gravity = run_path(str(root / "ch39-three-index-free-atoms.py"))
    thermal = run_path(str(root / "ch39-common-response-required-drawdown.py"))["field_bounds"]
    high_weak = run_path(str(root / "ch39-full-neutral-free-cutoff.py"))["weak_column"]
    left = 2_306_888_100_084
    moments, heat, fields = gravity["broken_fields"](root, (left,))
    _, _, pivots, _ = next(v for v in native["shell_weights"](root, left) if v[0] == 6)
    ctx.prec = 256
    c = [z.real for z in pivots[0]]
    assert sum(c).contains(1)
    low = native["weak_series"]()(arb(1))
    mp, mz = arb("1.22089e19"), arb("91.1876")
    high = high_weak(arb(1), (mz / mp)**2, (arb("1.176397e14") / mp)**2)

    def weak_bound(weights):
        depth = max([arb(0)] + [(-w).upper() for w in weights])
        return low + depth * (low-high)

    # 单谷的每个高截止集合由不交的前缀和后缀组成。空前缀和允许整条非降。
    valleys = [sum(c[:i], arb(0)) + sum(c[j:], arb(0)) for i in range(20) for j in range(i, 20)]
    valley_bound = weak_bound(valleys)
    # 峰在前十七层，任一非空高截止区间的左端最多为第十七层。
    early_peaks = [sum(c[i:j]) for i in range(17) for j in range(i+1, 20)]
    early_bound = weak_bound(early_peaks)
    print("单谷、峰在前十七层时各自的弱差上界", valley_bound.str(14), early_bound.str(14), flush=True)
    assert valley_bound < arb("-.00027260") and early_bound < 0

    # 峰在第 18 层，令 X=x18。归一化后比较 S2−XS1 的线性下界与二次上界。
    k, w, err_k, err_w = thermal(fields[0])
    total = sum(abs(v) for v in c)
    lower = (1-total*err_k)/k - total*err_w/w
    other_negative = sum(max(arb(0), (-v).upper()) for i, v in enumerate(c) if i != 17)
    upper = other_negative/4
    print("第十八层最高时 (S2−XS1)/X 的必要下界和可达上界",
          lower.str(14), upper.str(14), flush=True)
    assert lower > upper
    # 峰在第十九层已包括在单谷类（空前缀）。两次转向的上界由已发表的三段根给出。

    # 唯一下降在 r 与 r+1 之间时，每个高截止集合是两个有序段各自的后缀。
    candidates = []
    for r in range(1, 19):
        weights = [sum(c[a:r], arb(0)) + sum(c[b:], arb(0))
                   for a in range(r+1) for b in range(r, 20)]
        if not weak_bound(weights) < 0:
            candidates.append(r)
    print("单次下降时弱差未排除的位置", candidates, flush=True)
    assert candidates == [7, 8, 18]

    # 对两层连续提升的显式见证，只在第八层至第九层下降。
    weak_at = native["weak_series"]()
    p0, p1 = (mz/mp)**2, (arb("1.176397e14")/mp)**2
    high_pair = run_path(str(root / "ch39-ordered-cross-gravity-roots.py"))["high_pair"]
    three = run_path(str(root / "ch39-adjacent-ordered-charge-six-cutoffs.py"))["solve_three"]
    middle = c[6] + c[7] + c[17]
    weights = (1-middle-c[18], middle, c[18])
    def solve_pair(abc, low):
        return high_pair(abc, low, fields[0], heat, moments)
    _, _, _, minimum = three(left, weights, ("1","1.001"), solve_pair, weak_at, high_weak, mz, p0, p1)
    assert mz < minimum


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

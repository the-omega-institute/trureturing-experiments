"""单调绕数权的前缀分解与固定竞争相位对参考味三窗的全阶排除。"""
from collections import Counter
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


BANDS = ((".00815", ".01280"), (".13670", ".14130"))
CUTOFF = 5000
# 两个整数对分别标记精确的原生几何相位。
LOCKS = {(1, 0): (29189, 41), (0, 1): (205479, 101)}


def verify_phases(root, shell, bands):
    previous = run_path(str(root / "ch39-twentyseven-coefficient-window.py"))
    inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
    model = run_path(str(root / "ch39-mincut-three-window-angle.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 192
    phases = shell["geo"]["phases"](acb.zeta_zero(11661).imag)
    d = (phases[724] / phases[804]).arg()
    a, b = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    assert a > 0 and abs(b) < arb("1e-45")
    distances = (d, d - 2 * (b / a).atan())
    assert all(arb(".1496557638107688") < v < arb(".1496557638107690") for v in distances)
    low, high = arb("18.94088212"), arb("18.94414547")
    bottom, top = model["sign"]["density"]["eigen_band"]()
    # 正根 (c+sqrt(c²+4m))/2 随 c 改变至多 |Δc|，上边加 4δ。
    cap = arb(23) / 25 * arb(2) ** -20 / 14**2
    slack = arb("4.48e-9") - cap
    assert 0 < slack and low < bottom - slack / 3 < top + 4 * slack < high
    aperture = arb.pi() / 2
    counts = [0, 0]
    for sign, distance in zip((1, -1), distances, strict=True):
        pending = [(-aperture, aperture)]
        while pending:
            left, right = pending.pop()
            t = left.union(right)
            x, y = sign * t, sign * (t - distance)
            if any(l < t < r for l, r in bands):
                counts[0] += 1
            elif (previous["signs"](inertia, 2, low, x, y)[0] >= 16
                  or previous["signs"](inertia, 2, high, x, y)[1] < 16):
                counts[1] += 1
            else:
                assert right - left > arb("1e-10"), (left, right, sign)
                mid = (left + right) / 2
                pending.extend(((left, mid), (mid, right)))
    print("两条等味线按荷二周期折叠后的窄带内与惯性排除末段数", *counts)
    return distances


def geometric_tail(rows, identically_zero=False):
    """共振面逐面给 1，其余面给 1/(N |sin(α/2)|)。"""
    near = 0
    error = arb(0)
    for angle, count in rows:
        denom = abs((angle / 2).sin())
        if denom.contains(0):
            if identically_zero:
                assert angle.is_zero()
            near += count
        else:
            error += count / denom.abs_lower()
    return near, error


def kernel_sum(rows, n):
    """余弦平均的闭式在可能跨零处用 sinc 连续延拓，并与绝对界 1 取交。"""
    result = arb(0)
    for angle, count in rows:
        if angle.is_zero():
            result += count
            continue
        if n == 1:
            value = angle.cos()
        else:
            half = angle / 2
            div = half.sin()
            if div.contains(0):
                ratio = (2 * n + 1) * ((2 * n + 1) * half).sinc() / half.sinc()
            else:
                ratio = ((2 * n + 1) * half).sin() / div
            value = ((ratio - 1) / (2 * n)).intersection(arb(0, 1))
        result += count * value
    return result


def local_example(root, lib, tri, negative_branch_distance):
    """同一组单调权在全部指定软化量上将同一个合格相位选成严格局部极小。"""
    x = arb("-.01046")
    y = x + negative_branch_distance
    orders = (1, 61, 119, 121, 204)
    cap = arb("4.48e-9")
    delta = arb(0).union(cap)
    columns = []
    hessians = []
    gradients = [arb(0) for _ in range(4)]
    hessian = [arb(0) for _ in range(3)]
    for q in range(1, max(orders) + 1):
        bare = lib["coeff"](tri, q, arb(0))
        soft = [sum(z**q for z, spoke in tri[k] if spoke)
                for k in ((1, 0), (0, 1), (1, -1))]
        first = lib["derivatives"]([(q, bare)], x, y)
        response = lib["derivatives"]([(q, soft)], x, y)
        gradients = [a + v for a, v in zip(gradients, (*first[:2], *response[:2]), strict=True)]
        hessian = [a + v - delta * w for a, v, w in zip(hessian, first[2:], response[2:], strict=True)]
        if q in orders:
            columns.append([arb(1), *(v / q for v in gradients)])
            hessians.append([v / q for v in hessian])
    matrix = arb_mat(columns).transpose()
    assert not matrix.det().contains(0)
    weights = matrix.solve(arb_mat([[1], [0], [0], [0], [0]])).entries()
    assert all(w > arb(".004") for w in weights)
    h, j, k = [sum(w * v for w, v in zip(weights, vs, strict=True)) for vs in zip(*hessians)]
    assert h > 30 and (h - 30) * (k - 30) > j*j
    print("局部见证的前缀权、Hessian 一致大于 30I", *[w.str(12) for w in weights])
    physical = run_path(str(root / "ch39-reference-flavor-global-winding.py"))
    physical["physical_all_heights"](root, x, y, cap)


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    exact = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 192
    bands = [(arb(l), arb(r)) for l, r in BANDS]
    distances = verify_phases(root, shell, bands)
    distance = distances[0].union(distances[1])
    all_faces = exact["enumerate_faces"]()
    multiplicities = Counter((k, z) for k, z in all_faces.values() if k != (0, 0))
    assert len(all_faces) == 529 and sum(multiplicities.values()) == 96
    assert set(multiplicities.values()) == {3}
    assert all(multiplicities[k, exact["conjugate"](z)] == n for (k, z), n in multiplicities.items())
    # 与软壳脚本逐频率逐面核对整数相位及辐边权。不同整数相位的包围互不相交。
    tri = lib["terms"](shell)
    groups = []
    soft = 0
    for k in ((1, 0), (0, 1), (1, -1)):
        units = {z: exact["unit"](z) for freq, z in multiplicities if freq == k}
        recovered = Counter()
        for native, spoke in tri[k]:
            matches = [z for z, unit in units.items() if (unit - native).contains(0)]
            assert len(matches) == 1 and spoke in (0, 1)
            recovered[matches[0]] += 1
            soft += spoke
        for z, unit in units.items():
            n = multiplicities[k, z]
            assert recovered[z] == n
            groups.append((k, z, unit.arg(), n))
    assert soft == 6
    correction = 2 * soft * arb("4.48e-9")
    assert acb.zeta_zero(1).imag > 14 and shell["soft"]["point"]["study"]["offset"] > 0
    assert arb(23) / 25 * arb(2) ** -20 / 14**2 < arb("4.48e-9")
    margin = arb(1) / 4 + correction
    u, v = [-exact["unit"](LOCKS[k]).arg() for k in ((1, 0), (0, 1))]
    rival = [(arb(0) if z == LOCKS.get(k) else a + k[0] * u + k[1] * v, n)
             for k, z, a, n in groups]
    aligned, rival_error = geometric_tail(rival, identically_zero=True)
    assert aligned == 6 and rival_error < 6085
    ref = [arb(0)] + [kernel_sum(rival, n) for n in range(1, CUTOFF)]
    for l, r, lift in ((l, r, k) for l, r in bands for k in (0, 1)):
        # 共轭配对使单位辐边的作用在同时取反时不变，软化统一由 correction 支付。
        def at(t):
            angles = []
            for (i, j), _, a, n in groups:
                angle = a + (i + j) * t - j * distance
                if lift and i + j:
                    angle += -arb.pi() if a > 0 else arb.pi()
                angles.append((angle, n))
            return angles

        near, target_error = geometric_tail(at(l.union(r)))
        assert near <= 3 and target_error < 6800, (near, target_error)
        bound = aligned - near - (rival_error + target_error) / CUTOFF
        assert bound > margin
        pending = [(l, r, range(1, CUTOFF))]
        checked = 0
        while pending:
            left, right, orders = pending.pop()
            cells = at(left.union(right))
            undecided = []
            for n in orders:
                if ref[n] - kernel_sum(cells, n) > margin:
                    checked += 1
                else:
                    undecided.append(n)
            if undecided:
                assert right - left > arb("1e-10"), (left, right, undecided[:5])
                mid = (left + right) / 2
                pending.extend(((left, mid, undecided), (mid, right, undecided)))
        print("窄带、有限前缀认证数、候选尾界常数、尾部平均差下界",
              l.str(6), r.str(6), lift, checked, target_error.abs_upper().str(9), bound.lower().str(9))
    print("固定竞争点、竞争尾界常数、软化校正上界", u.str(15), v.str(15),
          rival_error.abs_upper().str(12), correction.str(9))
    local_example(root, lib, tri, distances[1])


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)

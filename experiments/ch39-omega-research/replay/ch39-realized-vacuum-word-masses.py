"""以 213 个严格正的两面字质量实现指定的十四项真空权，认证全局最小总变差。"""

from itertools import product
from pathlib import Path
from runpy import run_path
import sys

import numpy as np
from flint import arb, arb_mat, ctx


SUPPORT = (1, 2, 3, 4, 28, 30, 35, 38, 40, 65, 67, 83, 104, 105)
WINNER = (1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1)
OFF_SUPPORT_DUAL = (
    -135059845, 23774560, 8210337, -4733902, -17349575, 84725972, -56522793,
    -19496217, 7271045, 32923263, -35336800, 74404806, 24197122, -5041300,
    -33989819, 7260198, -65114907, -17630907, 1361860, 23436577, 351247,
    36195963, 4984628, -11959973, -5070824, 3410268, 490131, 5537013, -16088506,
    -1221318, -2722583, 25526015, -9443138, -1652252, 1244461, -2010424,
    -26274168, 21439255, 1203306, -2315600, 4203297, 20241213, -28664608,
    3649396, 7337324, -3744350, -8745404, 28315820, -11160458, -13922577,
    1385815, -4451319, -19028215, 18057161, 18280123, 14215569, -23188056,
    -19966996, -111275, -18567832, 26049541, 25787996, 20562893, 967685,
    20267834, -50130686, -20691562, -17631630, -1691779, -22816890, 60556120,
    4808840, -2708629, 25830042, -52116940, 36259605, 18087191, 14240589,
    -25650020, 31885101, -72135515, -37075608, -25116970, 20184942, -13918567,
    83642444, 32479391, 22287964, -11341134, 6287275, -58495848,
)


def candidates(c, w, columns):
    # 双精度只建议两个支点和其余对偶符号。每个下界另以区间运算验证。
    c, w = np.array(list(map(float, c))), np.array(list(map(float, w)))
    matrix = np.array([[float(x) for x in col] for col in columns])
    signs = np.array([(1, *tail) for tail in product((1, -1), repeat=len(columns)-1)])
    for batch in np.array_split(signs, 32):
        v = batch @ matrix
        rows = np.arange(len(batch))

        def at(scale):
            offset = c + scale[:, None]*v
            roots = -offset/w
            order = np.argsort(roots, axis=1)
            cumulative = np.cumsum(np.abs(w)[order], axis=1)
            med = np.argmax(cumulative >= np.sum(np.abs(w))/2, axis=1)
            pivot = order[rows, med]
            free = roots[rows, pivot]
            residual = offset + free[:, None]*w
            h = np.sign(residual)
            h[rows, pivot] = 0
            h[rows, pivot] = -np.sum(h*w, axis=1)/w[pivot]
            return np.sum(h*v, axis=1), residual, pivot

        lower = np.full(len(batch), -.05)
        upper = -lower
        for _ in range(44):
            mid = (lower+upper)/2
            slope = at(mid)[0]
            lower, upper = np.where(slope < 0, mid, lower), np.where(slope < 0, upper, mid)
        _, residual, first = at((lower+upper)/2)
        abs_residual = np.abs(residual)
        abs_residual[rows, first] = np.inf
        second = np.argmin(abs_residual, axis=1)
        for i, choice in enumerate(batch):
            yield (tuple(map(int, choice)), int(first[i]), int(second[i]),
                   tuple(map(int, np.where(residual[i] > 0, 1, -1))))


def signed_sum(row, signs):
    return sum(x if h > 0 else -x for x, h in zip(row, signs) if h)


def certificate(c, v, w, i, j, selected):
    # 两个支点上的 h 使 h·v=h·w=0，因此任意 x,y 均有 ||c+xv+yw||₁ ≥ h·c/||h||∞。
    signs = list(selected)
    signs[i] = signs[j] = 0
    hv, hw = signed_sum(v, signs), signed_sum(w, signs)
    det = v[i]*w[j]-v[j]*w[i]
    assert not det.contains(0)
    hi = (-hv*w[j]+hw*v[j])/det
    hj = (hv*w[i]-hw*v[i])/det
    norm = arb(1).max(abs(hi)).max(abs(hj))
    lower = (signed_sum(c, signs)+hi*c[i]+hj*c[j])/(2*norm)
    return lower, hi, hj, det


def verify(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    local = load("mincut-reference-ckm-action")
    winding = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    target = load("height-independent-flavor-vacuum")
    global_model = load("reference-flavor-global-winding")
    geometry = load("mincut-real-triangle-tomography")
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    tri = winding["terms"](shell)
    weights = target["minimum_order_weights"](winding, tri, x, y)
    assert tuple(sorted(weights)) == SUPPORT and (sum(weights.values())-1).contains(0)
    theta = geometry["unit"]((108469, 109)).arg()
    A = arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)])
    inverse = A.inv()
    mult = [1]+[2]*106
    ref = [arb(n)/213 for n in mult]
    missing = [q for q in range(1, 106) if q not in SUPPORT]
    assert len(missing) == len(OFF_SUPPORT_DUAL) == 91
    dual = list(zip(missing, OFF_SUPPORT_DUAL))
    score = [sum(z*A[q, j] for q, z in dual) for j in range(107)]
    cutoff = abs(sum(h*r for h, r in zip(score, ref))) / (
        max(h.upper() for h in score)-min(h.lower() for h in score))
    assert cutoff > arb(".0031567504")
    c = [inverse[j, 0]-ref[j] for j in range(107)]
    w = [inverse[j, 106] for j in range(107)]
    columns = [[inverse[j, q]*weights[q].sqrt() for j in range(107)] for q in SUPPORT]
    previous = (0,)*len(SUPPORT)
    v = [arb(0)]*107
    other_lower = arb(1)
    planes = 0
    for signs, i, j, selected in candidates(c, w, columns):
        for old, new, column in zip(previous, signs, columns):
            if old != new:
                v = [a+(new-old)*b for a, b in zip(v, column)]
        # 增量更新仅涉及翻转的符号。v 始终包围精确的有符号列和。
        previous = signs
        planes += 1
        lower, hi, hj, det = certificate(c, v, w, i, j, selected)
        if signs != WINNER:
            other_lower = other_lower.min(lower)
            continue
        assert {i, j} == {3, 14} and abs(hi) < 1 and abs(hj) < 1
        scale = (-c[i]*w[j]+c[j]*w[i])/det
        free = (c[i]*v[j]-c[j]*v[i])/det
        change = [r+scale*a+free*b for r, a, b in zip(c, v, w)]
        assert all(selected[k]*z > 0 for k, z in enumerate(change) if k not in (i, j))
        masses = [u+z for u, z in zip(ref, change)]
        floor = min(m/n for m, n in zip(masses, mult))
        tv = sum(map(abs, change))/2
        assert (tv-lower).contains(0)
    assert planes == 8192 and other_lower > arb(".00334468")
    assert arb(".0033434696169") < tv < arb(".0033434696171")
    assert floor > arb(".0031330080774") and scale < 0
    moments = (A*arb_mat([[p] for p in masses])).entries()
    exact = {q: scale*h*weights[q].sqrt() for q, h in zip(SUPPORT, WINNER)}
    assert (moments[0]-1).contains(0) and (moments[106]-free).contains(0)
    assert all((moments[q]-exact.get(q, 0)).contains(0) for q in range(1, 106))
    norm = sum(m*m for m in moments[1:106])
    assert norm > 0 and (norm-scale**2).contains(0)
    assert all((moments[q]**2/norm-weights.get(q, 0)).contains(0) for q in range(1, 106))

    # 采用原证明覆盖全部正零点高度的共同软化区间，复核这一个权向量的全局比较。
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    modes = [(q, [a*z for z in winding["coeff"](tri, q, delta)]) for q, a in weights.items()]
    global_model["conjugation_symmetry"](root)
    drift = global_model["global_cover"](winding, tri, weights, modes, x, y)
    assert drift < arb("-.05")
    global_model["physical_all_heights"](root, x, y, cap)
    print("仅由 91 个零矩得到的总变差下界", cutoff.lower().str(16))
    print("符号平面数、其他平面最小对偶下界", planes, other_lower.lower().str(16))
    print("有符号尺度、第106矩、最小总变差、最小每字质量",
          *(a.str(18) for a in (scale, free, tv, floor)))
    print("归一前的矩平方总量、第105圈份额", norm.str(16), weights[105].str(16))
    print("依次为 j=0,…,106 的每字概率，共轭字使用相同概率")
    for start in range(0, 107, 10):
        print(start, *[(masses[j]/mult[j]).str(14) for j in range(start, min(start+10, 107))])


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)

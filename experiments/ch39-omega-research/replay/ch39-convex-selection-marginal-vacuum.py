"""三种严格凸选权在 105 圈共同真空的 100 维可行族上只能停在边缘真空。

模型对象沿用 105 圈共同真空：参考味负支上的位置 A(t)=(-t, D_- - t)，
0.00815 <= t <= 0.01280。P_t 为归一、非负且在全部高度驻于 A(t) 的权
（一行归一加四条梯度等式），F_t 为其中 A(t) 在全部有限高度与单位辐边极限
都是全局真空的权。已有证书给出 F_t 内 w_105 > 0.072。

本脚本对最大香农熵、最小二范数平方与最大对数体积（解析中心）三种泛函，
在 1024 个有理子段上同时给出
    P_t 上的最优值（可行点上界）  与  P_t ∩ {w_105 >= 0.072} 上的对偶下界，
从而证明每个泛函在 P_t 上的唯一最优权都不是共同真空，且 F_t 上的最优值
至少差出打印的余量。浮点只用于挑选可行点与乘子，所有不等式由 Arb 球算术认证。
"""
from pathlib import Path
from runpy import run_path
import json
import sys

import numpy as np
from flint import acb, arb, arb_mat, ctx

LO, HI = 815, 1280          # t 的端点（单位 1e-5）
PIECES = 1024
SHARE = "0.072"              # 已有证书给出的末圈必要份额
KINDS = ("entropy", "l2", "logvol")


def load(root, name):
    return run_path(str(root/("ch39-"+name+".py")))


def build_rows(root):
    """每圈五个约束行写成 K + q Im(h e^{-iqt})，第 0 行为归一。"""
    lib = load(root, "mincut-global-winding")
    shell = load(root, "three-ray-chiral-boundary")
    angles = load(root, "mincut-reference-ckm-action")
    flavor = load(root, "chiral-neutral-flavor-phase")
    fixed = load(root, "height-independent-flavor-vacuum")
    ctx.prec = 160
    tri = lib["terms"](shell)
    d = -angles["reference_angles"](flavor)[0]
    keys = ((1, 0), (0, 1), (1, -1))
    rows = []
    for q in range(1, 106):
        bare = lib["coeff"](tri, q, arb(0))
        soft = [sum((z**q for z, m in tri[k] if m), acb(0)) for k in keys]
        row = [(arb(1), acb(0))]
        for cx, cy, cd in (bare, soft):
            diag = cd*acb(0, -q*d).exp()
            row.append((q*diag.imag, cx))
            row.append((-q*diag.imag, cy*acb(0, q*d).exp()))
        rows.append(row)

    def entry(q, k, t, der=0):
        const, h = rows[q][k]
        z = h*acb(0, -(q+1)*t).exp()
        return const+(q+1)*z.imag if der == 0 else -(q+1)**2*z.real

    # 与原有梯度列逐项重合。
    for t in (arb("0.00815"), arb("0.01046"), arb("0.0128")):
        for q in (0, 39, 103, 104):
            ref = fixed["gradient_columns"](lib, tri, q+1, -t, d-t)
            assert all(entry(q, k+1, t).overlaps(ref[k]) for k in range(4))
    return entry, d


def columns(entry, t):
    return [[entry(q, k, t) for k in range(5)] for q in range(105)]


def float_matrix(entry, t):
    return np.array([[float(entry(q, k, arb(t)).mid()) for q in range(105)] for k in range(5)])


# ---------- 浮点：只用来选可行点与对偶乘子 ----------
LAST = (np.arange(105) == 104).astype(float)


def ascend(dual, p):
    for _ in range(400):
        v, g, h = dual(p)
        step = np.linalg.solve(h-1e-15*np.eye(len(p)), -g)
        s = 1.0
        while s > 1e-18:
            v2 = dual(p+s*step)[0]
            if np.isfinite(v2) and v2 >= v+1e-4*s*(g@step):
                break
            s /= 2
        p = p+s*step
        if np.abs(g).max() < 1e-14 or s <= 1e-18:
            break
    return p


def float_optimum(kind, m, share=None):
    """返回原始最优点与（未缩放）对偶乘子；share 给出时加 w_105 >= share。"""
    scale = np.abs(m).max(1)
    ms = m/scale[:, None]
    extra = [LAST] if share is not None else []
    if kind == "entropy":
        x = np.vstack([ms[1:], *extra])
        r = np.r_[np.zeros(4), [share] if share is not None else []]

        def dual(p):
            a = p@x
            top = a.max()
            e = np.exp(a-top)
            w = e/e.sum()
            mean = x@w
            return p@r-top-np.log(e.sum()), r-mean, -((x*w)@x.T-np.outer(mean, mean))
        start = np.zeros(len(r))
        primal = lambda p: np.exp(p@x-(p@x).max())/np.exp(p@x-(p@x).max()).sum()
        k, sc = 4, scale[1:]
    else:
        x = np.vstack([ms, *extra])
        r = np.r_[1.0, np.zeros(4), [share] if share is not None else []]
        if kind == "l2":
            def dual(p):
                a = p@x
                pos = np.maximum(a, 0)
                return p@r-(pos**2).sum()/4, r-x@pos/2, -(x*(a > 0))@x.T/2
            start = np.r_[2/105, np.zeros(len(r)-1)]
            primal = lambda p: np.maximum(p@x, 0)/2
        else:
            def dual(p):
                a = p@x
                if (a >= 0).any():
                    return -np.inf, None, None
                return p@r+(1+np.log(-a)).sum(), r+x@(1/a), -(x/a**2)@x.T
            start = np.r_[-105.0, np.zeros(len(r)-1)]
            primal = lambda p: -1/(p@x)
        k, sc = 5, scale
    p = ascend(dual, start)
    w = primal(p)
    p = p.copy()
    p[:k] = p[:k]/sc
    return w, p


# ---------- Arb 证书 ----------
def objective(kind, w):
    if kind == "entropy":
        return sum(v*v.log() for v in w if v != 0)
    if kind == "l2":
        return sum(v*v for v in w)
    return -sum(v.log() for v in w)


def primal_upper(kind, entry, mid, rad, m):
    """对子段内每个 t 构造 P_t 中同一族可行点，给出其目标值的上界。"""
    w, _ = float_optimum(kind, m)
    base = [arb(float(f"{v:.10e}")) if v > 1e-300 else arb(0) for v in w]
    c0, ct = columns(entry, mid), columns(entry, mid+arb(0, rad))
    # w(t) = base*(1 + C0^T y(t))，y 由 5x5 系统 M(t) diag(base) C0^T y = e1 - M(t) base 决定。
    gram = arb_mat([[sum(ct[q][i]*base[q]*c0[q][j] for q in range(105)) for j in range(5)]
                    for i in range(5)])
    rhs = arb_mat([[arb(int(i == 0))-sum(ct[q][i]*base[q] for q in range(105))] for i in range(5)])
    y = gram.solve(rhs).entries()
    shift = [sum(c0[q][j]*y[j] for j in range(5)) for q in range(105)]
    assert max(s.abs_upper() for s in shift) < arb("0.5")
    weights = [b*(1+s) for b, s in zip(base, shift)]
    assert all(v > 0 or v == 0 for v in weights)
    return objective(kind, weights).upper()


def dual_lower(kind, entry, mid, rad, m, share):
    """P_t ∩ {w_105 >= share} 上目标值的拉格朗日下界，对子段内全部 t 成立。"""
    _, p = float_optimum(kind, m, float(share.mid()))
    mult = [arb(float(v)) for v in p]
    nu = mult[-1]
    assert nu > 0
    col = columns(entry, mid+arb(0, rad))
    if kind == "entropy":
        a = [sum(mult[k]*col[q][k+1] for k in range(4))+(nu if q == 104 else 0) for q in range(105)]
        return (nu*share-sum(z.exp() for z in a).log()).lower()
    a = [sum(mult[k]*col[q][k] for k in range(5))+(nu if q == 104 else 0) for q in range(105)]
    if kind == "l2":
        return (mult[0]+nu*share-sum(arb(max(0, z.upper()))**2 for z in a)/4).lower()
    assert all(z < 0 for z in a)
    return (mult[0]+nu*share+sum(1+(-z).log() for z in a)).lower()


def flagship_bracket(root, entry):
    """x=-0.01046 处真实可行族上的两侧证书：竞争点对偶界与一个经全环面覆盖的共同真空权。"""
    lib = load(root, "mincut-global-winding")
    shell = load(root, "three-ray-chiral-boundary")
    angles = load(root, "mincut-reference-ckm-action")
    flavor = load(root, "chiral-neutral-flavor-phase")
    fixed = load(root, "height-independent-flavor-vacuum")
    cover = load(root, "reference-flavor-global-winding")
    ctx.prec = 224
    tri = lib["terms"](shell)
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    cols = {q: [arb(1), *fixed["gradient_columns"](lib, tri, q, x, y)] for q in range(1, 106)}
    bare = {q: lib["coeff"](tri, q, arb(0)) for q in range(1, 106)}
    home = {q: lib["energy"](q, bare[q], x, y) for q in range(1, 106)}
    data = json.loads(Path(__file__).with_suffix(".json").read_text())

    def complete(free, basis):
        rhs = [arb(int(k == 0))-sum(v*cols[q][k] for q, v in free.items()) for k in range(5)]
        mat = arb_mat([cols[q] for q in basis]).transpose()
        return {**free, **dict(zip(basis, mat.solve(arb_mat([[v] for v in rhs])).entries()))}

    # 已有 100 维族中心：全部 105 个权严格为正的共同真空权。
    original = fixed["minimum_order_weights"](lib, tri, x, y)
    solved = (30, 35, 38, 40, 105)
    interior = complete({q: original.get(q, arb("1e-8")) for q in range(1, 106) if q not in solved}, solved)
    assert all(v > 0 for v in interior.values())
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    result = {}
    for kind in ("entropy", "l2"):
        rec = data[kind]
        mult = [arb(v) for v in rec["multipliers"]]
        rows = []
        for bx, by, nu in rec["cuts"]:
            nu = arb(nu)
            assert nu > 0
            bx, by = arb(bx), arb(by)
            rows.append((nu, {q: lib["energy"](q, bare[q], bx, by)-home[q] for q in range(1, 106)}))
        # 对 F 中任一 w，每个竞争点给出 sum_q w_q [S_q(B)-S_q(A)] >= 0（单位辐边极限）。
        if kind == "entropy":
            a = [sum(m*c for m, c in zip(mult, cols[q][1:]))+sum(n*r[q] for n, r in rows)
                 for q in range(1, 106)]
            bound = sum(z.exp() for z in a).log()           # 熵的上界
        else:
            a = [sum(m*c for m, c in zip(mult, cols[q]))+sum(n*r[q] for n, r in rows)
                 for q in range(1, 106)]
            bound = mult[0]-sum(arb(max(0, z.upper()))**2 for z in a)/4   # 二范数平方的下界
        base = complete({int(q): arb(v) for q, v in rec["free"].items()}, rec["basis"])
        assert all(v > 0 or v == 0 for v in base.values())
        mix = arb(rec["mix"])
        weights = {q: (1-mix)*base[q]+mix*interior[q] for q in range(1, 106)}
        modes = [(q, [w*c for c in lib["coeff"](tri, q, delta)]) for q, w in weights.items()]
        drift = cover["global_cover"](lib, tri, weights, modes, x, y)
        assert drift < 0
        value = (-sum(v*v.log() for v in weights.values()) if kind == "entropy"
                 else sum(v*v for v in weights.values()))
        top = sorted(weights, key=lambda q: -weights[q].mid())[:6]
        result[kind] = (bound, value, weights[105])
        print("旗舰相位：泛函、可行族最优值的两侧界、可行见证的 w_105、最大六权的绕数",
              kind, bound.str(10), value.str(10), weights[105].str(6), top, flush=True)
    assert arb("2.6193") < result["entropy"][1] and result["entropy"][0] < arb("2.6291")
    assert arb("0.091428") < result["l2"][0] and result["l2"][1] < arb("0.092585")
    # 同一相位上 100 维切片本身的最优值（可行点）作对照。
    m = float_matrix(entry, 0.01046)
    slice_h = -primal_upper("entropy", entry, arb("0.01046"), arb(0), m)
    slice_l2 = arb(primal_upper("l2", entry, arb("0.01046"), arb(0), m))
    assert arb("4.2466") < slice_h and slice_l2 < arb("0.017062")
    print("同一相位切片：最大熵下界、最小二范数平方上界", slice_h.str(8), slice_l2.str(8))


def verify(root):
    load(root, "all-native-spectral-masses")["verify"](root)
    load(root, "fixed-vacuum-weight-dimension")["verify"]()
    entry, _ = build_rows(root)
    flagship_bracket(root, entry)
    ctx.prec = 160
    share = arb(SHARE)
    lo, hi = arb(LO)/10**5, arb(HI)/10**5
    rad = (hi-lo)/(2*PIECES)
    report = {}
    for kind in KINDS:
        gap, worst_slice, worst_vac = None, None, None
        for i in range(PIECES):
            mid = lo+(2*i+1)*rad
            m = float_matrix(entry, float(mid.mid()))
            upper = primal_upper(kind, entry, mid, rad, m)
            lower = dual_lower(kind, entry, mid, rad, m, share)
            g = lower-upper
            gap = g if gap is None else min(gap, g)
            worst_slice = upper if worst_slice is None else max(worst_slice, upper)
            worst_vac = lower if worst_vac is None else min(worst_vac, lower)
        report[kind] = (gap, worst_slice, worst_vac)
        print("泛函、最小认证余量、切片最优值的上界、真空权目标值的下界",
              kind, arb(gap).str(8), arb(worst_slice).str(8), arb(worst_vac).str(8), flush=True)
    assert report["entropy"][0] > arb("0.0379")
    assert report["l2"][0] > arb("0.00114")
    assert report["logvol"][0] > arb("7.51")
    # 全部共同真空权：香农熵 < 4.2260 奈特，二范数平方 > 0.018744，负对数体积 > 537.93。
    assert -report["entropy"][2] < arb("4.2260")
    assert report["l2"][2] > arb("0.018744")
    assert report["logvol"][2] > arb("537.93")
    # 数值参考（非证书）：t 端点处三种切片最优权的末圈份额都远低于 0.072。
    for t in (LO/1e5, HI/1e5):
        m = float_matrix(entry, t)
        shares = [float_optimum(kind, m)[0][104] for kind in KINDS]
        assert max(shares) < 0.034
        print("数值：t、三种切片最优权的 w_105", t, *[f"{s:.6f}" for s in shares])


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else here)

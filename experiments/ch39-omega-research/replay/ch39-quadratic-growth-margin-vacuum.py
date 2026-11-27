"""最大化二次增长裕度的选权在旗舰相位选出严格共同真空。

模型对象沿用 105 圈共同真空（ch39-convex-selection-marginal-vacuum.py）：参考味负支上的
相位 A=(x,y)=(-t, D_- - t)，旗舰 t0=0.01046。P 为归一、非负且四条梯度等式成立的权。

对 w in P 定义单位辐边极限的二次增长裕度
    g(w) = inf_{B != ±A} [S_w(B;0)-S_w(A;0)] / rho(B)^2,  rho(B)^2 = min(rho_A(B)^2, rho_{-A}(B)^2),
    rho_C(B)^2 = sum_{k in K} (2 - 2 cos k.(B-C)),
K 取三个钟面方向 {(1,0),(0,1),(1,-1)}（主度量）或两条钟轴 {(1,0),(0,1)}（对照度量）。
g 是一族关于 w 的线性函数的下确界，所以 max_P g 是半无限线性规划。

本脚本给出：
  1. 一个精确落在 P 中的有理见证权 W，用全环面二阶泰勒覆盖证明 g(W) >= L；
  2. 有限个竞争点与方向上的对偶证书，证明对一切 w in P，g(w) <= U；
  3. 对一切 g(w) >= L 的 w（包含全部最优权）给出 w_105、w_104、w_1、支撑外质量与镜像漂移的区间；
  4. 由 L、漂移上界与软辐边常数推出：每个这样的 w 在全部高度 delta in [0, delta_cap] 都是严格共同真空；
  5. 只取曲率的选权：一切最优权都使某个竞争点低于 A，因而不是共同真空；
  6. 只取能隙（半径 0.02 外）的选权：数值最优权经 global_cover 认证为严格共同真空（选权最优性仅数值）。
浮点只用于挑选见证、竞争点与乘子，全部不等式由 Arb 球算术认证（精度 160 位）。
"""
from functools import cache
from pathlib import Path
from runpy import run_path
import json
import sys

import numpy as np

from flint import acb, arb, arb_mat, ctx

PREC = 160
KEYS = ((1, 0), (0, 1), (1, -1))


def load(root, name):
    return run_path(str(root/("ch39-"+name+".py")))


def upper_max(values):
    """一组球的上端点的最大值（逐个比较精确端点）。"""
    best = None
    for v in values:
        u = v.upper()
        best = u if best is None else best.max(u)
    return best


def margin_cover(modes, x, y, L, metric, box=arb("0.00002"), limit=2**24):
    """证明 h(B)=S_0(B)-S_0(A)-L*min(rho_A^2,rho_{-A}^2) >= 0 于全环面。

    h >= h_C := S_0-S_0(A)-L rho_C^2（C=±A），h_C 沿三条轴可分。每条轴上的一维函数用格中心的值与
    导数、格上二阶导数的球界作二维二阶泰勒下界。叶格三类：某个 h_C 的下界为正；落在 ±A 的方框内；
    完全落在 rho_C < rho_{-C} 一侧且 grad h_C 不含零（此处 h=h_C，不可能是 h 的极小点）。
    方框内 Hess h_C >= Hess S_0 - 2L K 正定，且 h_C(C)=0、grad h_C(C)=0（P 的等式与偶性），故 h_C >= 0。"""
    centre = (x, y, x-y)
    naxes = 3 if metric == 3 else 2
    trial = arb(0)
    for q, cs in modes:
        for c, t in zip(cs, centre):
            trial += (c*(1-acb(0, q*t).exp())).real

    @cache
    def axis_data(axis, n, index):
        if axis == 2:
            tc, r = arb.pi()*arb(2*index)/n, 2*arb.pi()/n
        else:
            tc, r = arb.pi()*(arb(2*index+1)/n-1), arb.pi()/n
        t = tc+arb(0, r)
        f = f1 = f2 = arb(0)
        for q, cs in modes:
            c = cs[axis]
            rc, ri = c*acb(0, q*tc).exp(), c*acb(0, q*t).exp()
            f += c.real-rc.real
            f1 += q*rc.imag
            f2 += q*q*ri.real
        out = []
        for sg in (1, -1):
            if axis < naxes:
                a = sg*centre[axis]
                uc, ui = tc-a, t-a
                g, g1, g2, rho = f-L*(2-2*uc.cos()), f1-2*L*uc.sin(), f2-2*L*ui.cos(), 2-2*ui.cos()
            else:
                g, g1, g2, rho = f, f1, f2, arb(0)
            out.append((g, g1, arb(min(0, g2.lower()))*r*r/2, g1+arb(0, (abs(g2)*r).upper()), rho))
        return out

    def lower(hs, s, n):
        # 二维泰勒：中心值 + 合成梯度的线性项下界 + 各轴二阶项下界（第三轴 δ = dx-dy）。
        r = arb.pi()/n
        gx, gy = hs[0][s][1]+hs[2][s][1], hs[1][s][1]-hs[2][s][1]
        return sum(h[s][0]+h[s][2] for h in hs)-(abs(gx)+abs(gy))*r

    pending, counts = [(1, 0, 0)], [0, 0, 0]
    while pending:
        n, i, j = pending.pop()
        hs = (axis_data(0, n, i), axis_data(1, n, j), axis_data(2, n, i-j))
        if any(lower(hs, s, n)-trial > 0 for s in (0, 1)):
            counts[0] += 1
            continue
        a, b = [arb.pi()*(arb(2*k+1, 1)/n-1) for k in (i, j)]
        if any(abs(a-s*x) < box and abs(b-s*y) < box for s in (-1, 1)):
            counts[1] += 1
            continue
        ra, rm = sum(h[0][4] for h in hs), sum(h[1][4] for h in hs)
        excluded = False
        for s, side in ((0, ra < rm), (1, rm < ra)):
            if side:
                px, py = hs[0][s][3]+hs[2][s][3], hs[1][s][3]-hs[2][s][3]
                excluded |= not px.contains(0) or not py.contains(0)
        if excluded:
            counts[2] += 1
            continue
        assert n < limit, (n, i, j, a, b)
        pending.extend((2*n, 2*i+u, 2*j+v) for u in (0, 1) for v in (0, 1))
    kxx, kxy, kyy = (2, -1, 2) if metric == 3 else (1, 0, 1)
    for s in (-1, 1):
        hxx, hxy, hyy = hessian(modes, s*x+arb(0, box), s*y+arb(0, box))
        p, r, d = hxx-2*L*kxx, hyy-2*L*kyy, hxy-2*L*kxy
        assert p > 0 and p*r > d*d
    return counts


def hessian(modes, x, y):
    hxx = hxy = hyy = arb(0)
    for q, cs in modes:
        u, v, w = [c*acb(0, q*t).exp() for c, t in zip(cs, (x, y, x-y))]
        hxx += q*q*(u+w).real
        hxy -= q*q*w.real
        hyy += q*q*(v+w).real
    return hxx, hxy, hyy


class Model:
    def __init__(self, root):
        self.lib = load(root, "mincut-global-winding")
        shell = load(root, "three-ray-chiral-boundary")
        angles = load(root, "mincut-reference-ckm-action")
        flavor = load(root, "chiral-neutral-flavor-phase")
        self.fixed = load(root, "height-independent-flavor-vacuum")
        self.cover = load(root, "reference-flavor-global-winding")
        ctx.prec = PREC
        self.tri = self.lib["terms"](shell)
        assert sum(map(len, self.tri.values())) == 529
        self.d = angles["reference_angles"](flavor)[0]
        ctx.prec = PREC
        self.bare = {q: self.lib["coeff"](self.tri, q, arb(0)) for q in range(1, 106)}
        self.soft = {q: [sum((z**q for z, m in self.tri[k] if m), acb(0)) for k in KEYS] for q in range(1, 106)}
        self.cap = arb(23)/25*arb(2)**-20/14**2
        self.set_phase(arb("0.01046"))

    def set_phase(self, t):
        """A(t) = (-t, D_- - t)；t 可以是球，此时下面的量对球内每个 t 都成立。"""
        self.x = -t
        self.y = self.x-self.d
        self.cols = {q: [arb(1), *self.fixed["gradient_columns"](self.lib, self.tri, q, self.x, self.y)]
                     for q in range(1, 106)}
        self.home = {q: self.lib["energy"](q, self.bare[q], self.x, self.y) for q in range(1, 106)}
        # 镜像漂移系数：(S_q(A;δ)-S_q(-A;δ))/δ = -2 Σ_k Im(s_kq) sin(q k·A)。
        angle = (self.x, self.y, self.x-self.y)
        self.mu = {q: -2*sum(s.imag*(q*a).sin() for s, a in zip(self.soft[q], angle)) for q in range(1, 106)}

    def complete(self, free, basis):
        rhs = [arb(int(k == 0))-sum(v*self.cols[q][k] for q, v in free.items()) for k in range(5)]
        mat = arb_mat([self.cols[q] for q in basis]).transpose()
        assert not mat.det().contains(0)
        sol = mat.solve(arb_mat([[v] for v in rhs])).entries()
        weights = {**free, **dict(zip(basis, sol))}
        assert all(v > 0 for v in weights.values())
        # 见证精确落在 P：归一与四条梯度等式的残差球都含零。
        assert all((sum(w*self.cols[q][k] for q, w in weights.items())-int(k == 0)).contains(0)
                   for k in range(5))
        return weights

    def rho2(self, metric, bx, by, sign):
        keys = KEYS if metric == 3 else KEYS[:2]
        return sum(2-2*(i*(bx-sign*self.x)+j*(by-sign*self.y)).cos() for i, j in keys)

    def rows(self, cert, metric):
        """每行给出 (a_q, b)：对一切 w in P，g(w)·b <= Σ_q w_q a_q（竞争点或 A 处方向）。"""
        out = []
        kxx, kxy, kyy = (2, -1, 2) if metric == 3 else (1, 0, 1)
        for kind, u, v in cert["rows"]:
            u, v = arb(u), arb(v)
            if kind == "p":
                ra, rm = self.rho2(metric, u, v, 1), self.rho2(metric, u, v, -1)
                assert ra < rm      # 竞争点严格靠近 A，min(rho_A^2, rho_{-A}^2) = rho_A^2
                a = {q: self.lib["energy"](q, self.bare[q], u, v)-self.home[q] for q in range(1, 106)}
                assert ra > 0
                out.append((a, ra))
            else:
                a = {}
                for q in range(1, 106):
                    _, _, hxx, hxy, hyy = self.lib["derivatives"]([(q, self.bare[q])], self.x, self.y)
                    a[q] = (u*u*hxx+2*u*v*hxy+v*v*hyy)/2
                out.append((a, u*u*kxx+2*u*v*kxy+v*v*kyy))
                assert out[-1][1] > 0
        return out

    def dual_value(self, mult, rows, phi, level):
        """对一切 w in P 且对每行 Σ_q w_q a_q >= level·b：
        φ·w <= λ_0 + max(0, max_q(φ_q + Σ ν a_q - λ·M_q)) - level Σ ν b。"""
        lam = [arb(v) for v in mult["lam"]]
        nu = {int(j): arb(v) for j, v in mult["nu"].items()}
        assert all(v >= 0 for v in nu.values())
        coeff = [phi.get(q, arb(0))+sum(n*rows[j][0][q] for j, n in nu.items())
                 - sum(l*c for l, c in zip(lam, self.cols[q])) for q in range(1, 106)]
        top = upper_max(coeff).max(arb(0))
        return lam[0]+top-level*sum(n*rows[j][1] for j, n in nu.items()), sum(n*rows[j][1] for j, n in nu.items())


def ratio_certificate(model, cert, root, enclose=True):
    metric, L = int(cert["metric"]), arb(cert["L"])
    model.set_phase(arb(cert["t"]))
    W = model.complete({int(q): arb(v) for q, v in cert["free"].items()}, cert["basis"])
    modes = [(q, [W[q]*c for c in model.bare[q]]) for q in W]
    counts = margin_cover(modes, model.x, model.y, L, metric)
    rows = model.rows(cert, metric)
    # 上界：Σν b = N > 0，g(w) N <= Σ_q w_q Σν a_q <= λ_0 + max⁺(...)。
    numer, norm = model.dual_value(cert["upper"], rows, {}, arb(0))
    assert norm > 0
    U = numer/norm
    assert L < U
    enc = {}
    for key, phi in (("w105", {105: arb(1)}), ("-w105", {105: arb(-1)}), ("w104", {104: arb(1)}),
                     ("-w104", {104: arb(-1)}), ("w1", {1: arb(1)}), ("-w1", {1: arb(-1)}),
                     ("off", {q: arb(1) for q in range(1, 106) if q not in W}),
                     ("drift", model.mu)) if enclose else ():
        enc[key] = model.dual_value(cert["enclose"][key], rows, phi, L)[0]
    # 曲率读数：Hess S_W(A;0) - 2·200·K 正定，即 W 在 A 处的方向曲率比 κ(W) > 200。
    hxx, hxy, hyy = hessian(modes, model.x, model.y)
    kxx, kxy, kyy = (2, -1, 2) if metric == 3 else (1, 0, 1)
    p, r, d = hxx-400*kxx, hyy-400*kyy, hxy-400*kxy
    assert p > 0 and p*r > d*d
    drift_w = sum(W[q]*model.mu[q] for q in W)
    return dict(W=W, L=L, U=U, enc=enc, counts=counts, modes=modes, drift=drift_w)


def strict_at_all_heights(model, L, drift_max, metric):
    """书面引理的数值前提：对 w in P，若 g(w) >= L 且漂移 <= drift_max < 0，则对一切
    δ in [0, δ_cap]，A 是 S_w(·;δ) 的严格全局极小（δ=0 时除去 -A）且 Hessian 正定。"""
    hs = upper_max([arb(q*q)*(abs(s[0])+abs(s[1])+2*abs(s[2])) for q, s in model.soft.items()])
    grads = []
    for q in range(1, 106):
        gx, gy, *_ = model.lib["derivatives"]([(q, model.soft[q])], -model.x, -model.y)
        grads.append((gx*gx+gy*gy).sqrt())
    G = upper_max(grads)
    cap, pi = model.cap, arb.pi()
    Lp = L-cap*pi*pi*hs/8
    assert drift_max < 0 and Lp > 0
    assert cap*(pi/2)**2*G*G < 4*Lp*(-drift_max)
    assert 2*L-cap*hs > 0       # Hess S(A;δ) >= 2L·λ_min(K) - δ H_s，λ_min(K) = 1
    return hs, G


def curvature_certificate(model, cert):
    """只取曲率 κ(w)=min_u u'Hess S_w(A;0)u/(2u'Ku) 的选权：一切最优权在 B0 处低于 A。"""
    # 多面体界只对 A 处方向行成立：竞争点行约束的是 g 而不是 κ。
    assert all(r[0] != "p" for r in cert["rows"])
    rows = model.rows(cert, 3)
    W = model.complete({}, cert["support"])
    KL = arb(cert["KL"])
    hxx, hxy, hyy = hessian([(q, [W[q]*c for c in model.bare[q]]) for q in W], model.x, model.y)
    p, r, d = hxx/2-2*KL, hyy/2-2*KL, hxy/2+KL
    assert p > 0 and p*r > d*d          # κ(W) > KL
    numer, norm = model.dual_value(cert["upper"], rows, {}, arb(0))
    Uc = numer/norm
    bx, by = [arb(v) for v in cert["B0"]]
    phi = {q: model.lib["energy"](q, model.bare[q], bx, by)-model.home[q] for q in range(1, 106)}
    worst = model.dual_value(cert["enclose"], rows, phi, KL)[0]
    own = sum(W[q]*phi[q] for q in W)
    return KL, Uc, worst, own, W


def empty_pieces(model, pieces, lo, hi):
    """每段 t 球上：对一切 t 与一切 w in P_t，g_t(w) <= U < 0，故 A(t) 在单位辐边极限不是全局最低，F_t 为空。"""
    from fractions import Fraction
    pieces = sorted(pieces, key=lambda p: Fraction(p["lo"]))
    assert Fraction(pieces[0]["lo"]) == Fraction(lo) and Fraction(pieces[-1]["hi"]) == Fraction(hi)
    assert all(Fraction(a["hi"]) == Fraction(b["lo"]) for a, b in zip(pieces, pieces[1:]))
    worst, balls = None, 0
    for piece in pieces:
        mult = dict(lam=piece["lam"], nu={str(j): v for j, v in enumerate(piece["nu"])})
        pending = [(Fraction(piece["lo"]), Fraction(piece["hi"]))]
        while pending:
            a, b = pending.pop()
            # 同一组乘子对子球同样有效；球内高估过大时二分 t。
            model.set_phase((arb(a.numerator)/a.denominator).union(arb(b.numerator)/b.denominator))
            rows = model.rows(dict(rows=piece["rows"]), 3)
            numer, norm = model.dual_value(mult, rows, {}, arb(0))
            assert norm > 0
            U = numer/norm
            if not U < 0:
                assert b-a > Fraction(1, 10**9), (a, b, U)
                pending.extend([(a, (a+b)/2), ((a+b)/2, b)])
                continue
            balls += 1
            worst = U.upper() if worst is None else worst.max(U.upper())
    return len(pieces), balls, worst


def gap_witness(model, cert):
    model.set_phase(arb("0.01046"))
    W = model.complete({int(q): arb(v) for q, v in cert["free"].items()}, cert["basis"])
    delta = arb(0).union(model.cap)
    modes = [(q, [W[q]*c for c in model.lib["coeff"](model.tri, q, delta)]) for q in W]
    return W, model.cover["global_cover"](model.lib, model.tri, W, modes, model.x, model.y)


def gap_on_grid(model, W, radius, n=2048):
    """仅数值：网格上 Γ_r(W)=min_{rho(B)>=r}[S_W(B;0)-S_W(A;0)]（三方向度量）。"""
    grid = 2*np.pi*np.arange(n)/n-np.pi
    shift = (np.arange(n)[:, None]-np.arange(n)[None, :]+n//2) % n
    x, y = float(model.x.mid()), float(model.y.mid())
    axes = []
    for k in range(3):
        f = np.zeros(n)
        for q, w in W.items():
            c = complex(model.bare[q][k])
            f += float(w.mid())*(c*(1-np.exp(1j*q*grid))).real
        axes.append(f)
    energy = axes[0][:, None]+axes[1][None, :]+axes[2][shift]
    home = sum(float(w.mid())*float(model.home[q].mid()) for q, w in W.items())
    def rho(sx, sy):
        return ((2-2*np.cos(grid-sx))[:, None]+(2-2*np.cos(grid-sy))[None, :]
                + (2-2*np.cos(grid-(sx-sy)))[shift])
    far = np.minimum(rho(x, y), rho(-x, -y)) >= radius**2
    return (energy-home)[far].min()


def verify(root):
    load(root, "reference-flavor-global-winding")["conjugation_symmetry"](root)
    model = Model(root)
    data = json.loads(Path(__file__).with_suffix(".json").read_text())
    expect = {
        "ratio3": dict(L="3.08449", U="3.084492", w105=("0.111256", "0.111277"), w104=("0.023309", "0.023328"),
                       w1=("0.431379", "0.431387"), off="0.0000067", drift="-2.0538", support=11),
        "ratio2": dict(L="6.34872", U="6.348731", w105=("0.110870", "0.110874"), w104=("0.0247773", "0.0247786"),
                       w1=("0.431182", "0.431192"), off="0.0000106", drift="-2.0722", support=13),
    }
    results = {}
    for key in ("ratio3", "ratio2"):
        res = ratio_certificate(model, data[key], root)
        e, enc, metric = expect[key], res["enc"], int(data[key]["metric"])
        assert data[key]["L"] == e["L"]
        print("  U =", res["U"].str(12), flush=True)
        assert res["U"] < arb(e["U"])
        assert len(res["W"]) == e["support"]
        for name in ("w105", "w104", "w1"):
            lo, hi = e[name]
            assert arb(lo) < -enc["-"+name] and enc[name] < arb(hi)
        assert enc["off"] < arb(e["off"]) and enc["drift"] < arb(e["drift"])
        assert arb(e["w105"][0]) > arb("0.072")   # 与共同真空的必要份额相容
        hs, G = strict_at_all_heights(model, res["L"], enc["drift"], metric)
        assert hs < arb("66151") and G < arb("362.24")
        delta = arb(0).union(model.cap)
        dmodes = [(q, [res["W"][q]*c for c in model.lib["coeff"](model.tri, q, delta)]) for q in res["W"]]
        mirror = model.cover["global_cover"](model.lib, model.tri, res["W"], dmodes, model.x, model.y)
        assert mirror < 0 and mirror.overlaps(res["drift"])
        results[key] = res
        print("度量、裕度证书 L < g* <= U、覆盖叶格（能量、方框、梯度）", metric, res["L"].str(8),
              res["U"].str(10), res["counts"], flush=True)
        print("  全部 g>=L 的权：w_105、w_104、w_1 的区间，支撑外质量上界，镜像漂移上界",
              (-enc["-w105"]).str(7), enc["w105"].str(7), (-enc["-w104"]).str(7), enc["w104"].str(7),
              (-enc["-w1"]).str(7), enc["w1"].str(7), enc["off"].str(5), enc["drift"].str(6))
        print("  见证支撑大小、见证 w_105、w_104、全高度引理常数 H_s、G",
              len(res["W"]), res["W"][105].str(10), res["W"][104].str(10), hs.str(6), G.str(6), flush=True)
    assert arb("3.08449") < results["ratio3"]["U"]
    # 上界证书中 ν>0 的竞争点（靠近 A 的代表）落在六个类附近（两行同属一类）。
    classes = [(0.0315, 3.116), (0.164, 0.136), (0.0128, 0.0420), (0.0107, 0.1487), (-0.1576, -0.0698), (-1.063, 0.815)]
    active = [data["ratio3"]["rows"][int(j)] for j in data["ratio3"]["upper"]["nu"]]
    hit = [min(range(6), key=lambda c: abs(float(r[1])-classes[c][0])+abs(float(r[2])-classes[c][1])) for r in active]
    assert all(abs(float(r[1])-classes[c][0]) < 1e-3 and abs(float(r[2])-classes[c][1]) < 1e-3 for r, c in zip(active, hit))
    assert sorted(set(hit)) == list(range(6)) and len(active) == 7
    # κ(W) > 200 与 g* < 3.084492：曲率比裕度大六十倍以上。
    assert arb(200)/results["ratio3"]["U"] > 60
    # 最近的活跃竞争点类距 A 约 0.023（平坦距离，取靠近 A 的代表）。
    near = [r for r, c in zip(active, hit) if c == 3]
    model.set_phase(arb("0.01046"))
    for r in near:
        dist = ((arb(r[1])-model.x)**2+(arb(r[2])-model.y)**2).sqrt()
        assert arb("0.0231") < dist < arb("0.0233")
    w3, w2 = results["ratio3"]["W"], results["ratio2"]["W"]
    assert set(w3) == {1, 2, 30, 32, 35, 38, 40, 67, 81, 104, 105}
    assert abs(w3[105]-arb("0.1112662768")) < arb("1e-10") and abs(w3[104]-arb("0.0233186589")) < arb("1e-10")
    assert abs(w3[1]-arb("0.43138541")) < arb("1e-8")
    assert set(w2) == {1, 3, 5, 32, 35, 38, 40, 65, 67, 81, 83, 104, 105}
    # 旗舰以外的相位：严格共同真空见证与裕度上界。
    side = {"0.00925": ("0.19393", "0.193940"), "0.0095": ("1.76969", "1.769694"), "0.011": ("3.11111", "3.111121"),
            "0.012": ("1.24871", "1.248713"), "0.012125": ("0.18569", "0.185696")}
    for t, (low, up) in side.items():
        cert = data["phase"][t]
        assert cert["t"] == t and cert["L"] == low
        res = ratio_certificate(model, cert, root, enclose=False)
        assert res["U"] < arb(up)
        # 见证自身：g(W) >= L 由覆盖给出，漂移 m(W) 直接计算，代入全高度引理。
        strict_at_all_heights(model, res["L"], res["drift"], 3)
        assert len(res["W"]) == {"0.00925": 10, "0.0095": 12, "0.011": 12, "0.012": 11, "0.012125": 12}[t]
        results[t] = res
        print("相位 t、L < g*(t) <= U、见证支撑、w_105、w_104、见证漂移", t, res["L"].str(8), res["U"].str(10),
              len(res["W"]), res["W"][105].str(8), res["W"][104].str(8), res["drift"].str(6), flush=True)
    # t 作为变量时，旗舰相位不是最大裕度相位：g*(0.011) >= 3.11111 > U(0.01046)。
    assert results["0.011"]["L"] > results["ratio3"]["U"]
    for (lo, hi), key in ((("163/20000", "147/16000"), "low"), (("49/4000", "8/625"), "high")):
        count, balls, worst = empty_pieces(model, data["empty"][key], lo, hi)
        print("可行族为空的 t 段、乘子组数、t 球数、各球上界的最大值（自适应二分到为负即停）", lo, hi, count, balls, worst.str(6), flush=True)
        assert worst < 0 and (count, balls) == {"low": (42, 639), "high": (27, 170)}[key]
    model.set_phase(arb("0.01046"))
    KL, Uc, worst, own, wc = curvature_certificate(model, data["curvature"])
    assert data["curvature"]["KL"] == "632.2"
    print("  Uc =", Uc.str(12))
    assert Uc < arb("632.2016")
    assert worst < arb("-9.131") and own < arb("-9.1326")
    assert set(wc) == {22, 24, 39, 65, 105}
    assert abs(wc[105]-arb("0.1727213")) < arb("1e-7")
    print("只取曲率：632.2 < κ* <= U、一切最优权在 B0 处的能量差上界、见证自身的能量差、见证 w_105",
          Uc.str(10), worst.str(8), own.str(8), wc[105].str(8), flush=True)
    wg, drift = gap_witness(model, data["gap"])
    assert drift < 0 and len(wg) == 14
    # 不依赖 global_cover 的二阶步：g(W_gap) >= 0.004 由二次增长覆盖给出，再代入全高度引理。
    gap_modes = [(q, [wg[q]*c for c in model.bare[q]]) for q in wg]
    margin_cover(gap_modes, model.x, model.y, arb("0.004"), 3)
    drift_gap = sum(wg[q]*model.mu[q] for q in wg)
    strict_at_all_heights(model, arb("0.004"), drift_gap, 3)
    assert drift_gap < arb("-3.778")
    gamma = gap_on_grid(model, wg, 0.02)
    # 割平面线性规划的数值最优值是证书数据的记录值，这里不重算；网格值是见证的上估计。
    assert abs(data["gap"]["value"]-0.038444) < 1e-6 and 0.03894 < gamma < 0.03896
    print("仅数值：只取能隙的割平面最优值、见证在 2048 网格上的 Γ_0.02", round(data["gap"]["value"], 6), round(gamma, 6))
    assert abs(wg[105]-arb("0.193205")) < arb("1e-6") and abs(wg[104]-arb("0.044362")) < arb("1e-6")
    print("只取能隙（半径 0.02 外）：数值最优权的支撑、w_105、w_104，经全环面覆盖认证为严格共同真空",
          len(wg), wg[105].str(8), wg[104].str(8))


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else here)

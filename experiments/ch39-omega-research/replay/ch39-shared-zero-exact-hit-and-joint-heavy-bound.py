"""二十七项字典、三枚零点 #11661、#11934、#22027 共用同一组系数：精确完整命中的存在与局部维数，以及共用条件下的联合重块下界。

区间证书，448 位。
一、精确共用命中。完整命中写成 (w,u_n) 的 39 条实多项式方程（每零点 13 条：三部门质量平方的特征多项式系数、
    三个 CKM 模方、一个相位条件）。数据文件固定 72 个实坐标中的 33 个，脚本以前向自动微分在盒子上包围 39×39 雅可比，
    用 Krawczyk 判据认证余下 39 个坐标在盒内有唯一解；该解是精确的三零点共用完整命中。随后在盒上包围它的重—重上型块。
二、联合下界。对每个零点集合 N 与每个特征标选择 (k_n)_{n∈N}，认证 Σ_n θ_n H_{n,k_n} − Σ_{n,f} a_{n,f} F_{n,k_n,f} ≻ 0；
    由此任何在 N 上共用系数的命中（精确，或各质量相对误差 ≤ε）都满足 max_n ‖Y_hh^(n)‖_F² ≥ 下界。
"""
import json
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

HERE = Path(__file__).parent
prev = run_path(str(HERE/"ch39-per-leg-double-hop-transfer.py"))
single = run_path(str(HERE/"ch39-single-zero-heavy-block-double-hop-hit.py"))
study, adjoint, positive = prev["study"], prev["adjoint"], prev["positive"]
MIXING = run_path(str(HERE/"ch39-finite-orthogonal-ckm.py"))
PREV = json.loads((HERE/"ch39-per-leg-double-hop-transfer.json").read_text())
DATA = json.loads((HERE/"ch39-shared-zero-exact-hit-and-joint-heavy-bound.json").read_text())
ZEROS = (11661, 11934, 22027)

# 精确共用命中的重—重块：谱范数（Rayleigh 下界, Sylvester 上界）、Frobenius 区间、h² 系数 GeV² 区间。
HIT = {
    "11661": dict(spec=("66633.4", "66633.5"), frob=("86361.4", "86361.5"), quad=("5.899e40", "5.901e40")),
    "11934": dict(spec=("66488.6", "66488.7"), frob=("86193.4", "86193.5"), quad=("5.924e40", "5.926e40")),
    "22027": dict(spec=("65983.2", "65983.3"), frob=("85838.2", "85838.3"), quad=("5.748e40", "5.750e40")),
}
WNORM = ("6795043", "6795044")            # 精确共用命中的 ‖w‖
OLD_DIST = "1e-9"                          # 它与已有近似共用见证之 w 的相对距离上界
# 联合下界：max_n ‖Y_hh^(n)‖_F² 的下界（ε=0 与 ε=1e-3 两种），以及由此得到的谱范数下界。
JOINT = {
    "11661,11934,22027": ("6.216", "6.180", ".6233", ".6215"),
    "11661,11934": ("5.275", "5.244", ".5742", ".5725"),
    "11661,22027": ("5.653", "5.640", ".5944", ".5937"),
    "11934,22027": ("4.567", "4.553", ".5342", ".5334"),
}
EPS = arb("1e-3")


# ---------- 前向自动微分（值与对活动坐标的梯度，均为 acb 球） ----------
class DS:
    def __init__(s, v, g): s.v, s.g = v, g
    def __add__(a, b):
        return DS(a.v+b.v, [x+y for x, y in zip(a.g, b.g)]) if isinstance(b, DS) else DS(a.v+b, a.g)
    def __sub__(a, b):
        return DS(a.v-b.v, [x-y for x, y in zip(a.g, b.g)]) if isinstance(b, DS) else DS(a.v-b, a.g)
    def __mul__(a, b):
        if isinstance(b, DS):
            return DS(a.v*b.v, [x*b.v+a.v*y for x, y in zip(a.g, b.g)])
        return DS(a.v*b, [x*b for x in a.g])
    def conj(a): return DS(a.v.conjugate(), [x.conjugate() for x in a.g])
    def inv(a): return DS(1/a.v, [-x/(a.v*a.v) for x in a.g])
    def sqrt(a):
        v = a.v.real.sqrt()
        return DS(acb(v), [x/(2*v) for x in a.g])


class DM:
    def __init__(s, v, g): s.v, s.g = v, g
    def __mul__(a, b):
        if isinstance(b, DM):
            return DM(a.v*b.v, [ga*b.v+a.v*gb for ga, gb in zip(a.g, b.g)])
        return DM(a.v*b, [g*b for g in a.g])
    def __sub__(a, b): return DM(a.v-b, a.g)
    def H(a):
        t = lambda m: acb_mat([[m[j,i].conjugate() for j in range(3)] for i in range(3)])
        return DM(t(a.v), [t(g) for g in a.g])


def trace(m):
    return DS(sum((m.v[i,i] for i in range(3)), acb(0)), [sum((g[i,i] for i in range(3)), acb(0)) for g in m.g])


def det(m):
    e = lambda i, j: DS(m.v[i,j], [g[i,j] for g in m.g])
    return (e(0,0)*(e(1,1)*e(2,2)-e(1,2)*e(2,1)) - e(0,1)*(e(1,0)*e(2,2)-e(1,2)*e(2,0))
            + e(0,2)*(e(1,0)*e(2,1)-e(1,1)*e(2,0)))


def eye():
    return acb_mat([[int(i==j) for j in range(3)] for i in range(3)])


def targets():
    T = single["target_masses"]()
    lam = [[x*x for x in T[g]] for g in ("up", "down", "leptons")]
    s12, s23, s13 = [arb(a).sin() for a in (".22704", ".04275", ".00377")]
    c13 = (1-s13*s13).sqrt()
    return lam, [s12**2*c13**2, s23**2*c13**2, s13**2], s12**2*s23**2*s13**2*c13**2, arb("1.139")


def coefficient_matrix(blocks):
    # A(z)=D(z)·A0：第 (f,i,j) 行只含一个 Higgs 特征标 k（上型 k≡−i−j，下型与轻子 k≡i+j），D(z) 为对角的 z_k 或 z̄_k。
    rows = []
    for f in range(3):
        for i in range(3):
            for j in range(3):
                k = (-i-j) % 3 if f == 0 else (i+j) % 3
                rows.append([blocks[k][f][i,t] for t in range(27)])
    return acb_mat(rows)


def equations(A0s, x, active):
    """x 为 72 个实坐标：w 实部(27)、w 虚部(27)、每零点 u 实部(3)与虚部(3)；z=u/|u|。返回 (值, 梯度) 的 39 条实方程。"""
    m = len(active)
    slot = {a: i for i, a in enumerate(active)}
    def seed(idx, unit):
        g = [acb(0)]*m
        if idx in slot:
            g[slot[idx]] = unit
        return g
    w = [DS(acb(x[t], x[27+t]), [p+q for p, q in zip(seed(t, acb(1)), seed(27+t, acb(0,1)))]) for t in range(27)]
    lam, TV, kappa, delta = targets()
    vals, grads = [], []
    for a, n in enumerate(ZEROS):
        o = 54+6*a
        u = [DS(acb(x[o+k], x[o+3+k]), [p+q for p, q in zip(seed(o+k, acb(1)), seed(o+3+k, acb(0,1)))]) for k in range(3)]
        norm = (u[0]*u[0].conj() + u[1]*u[1].conj() + u[2]*u[2].conj()).sqrt()
        z = [uk*norm.inv() for uk in u]
        A0 = A0s[n]
        N = [DS(sum((A0[r,t]*w[t].v for t in range(27)), acb(0)),
                [sum((A0[r,t]*w[t].g[c] for t in range(27)), acb(0)) for c in range(m)]) for r in range(27)]
        mats = []
        for f in range(3):
            ent = [[None]*3 for _ in range(3)]
            for i in range(3):
                for j in range(3):
                    k = (-i-j) % 3 if f == 0 else (i+j) % 3
                    ent[i][j] = (z[k] if f == 0 else z[k].conj())*N[9*f+3*i+j]
            mats.append(DM(acb_mat([[ent[i][j].v for j in range(3)] for i in range(3)]),
                           [acb_mat([[ent[i][j].g[c] for j in range(3)] for i in range(3)]) for c in range(m)]))
        hs = [M*M.H() for M in mats]
        eqs = []
        for f, h in enumerate(hs):
            l = lam[f]
            t1, t2 = trace(h), trace(h*h)
            c2 = (t1*t1-t2)*acb(arb(1)/2)
            goal = (l[0]+l[1]+l[2], l[0]*l[1]+l[0]*l[2]+l[1]*l[2], l[0]*l[1]*l[2])
            eqs += [c*acb(1/g)-acb(1) for c, g in zip((t1, c2, det(h)), goal)]
        def proj(h, l, i):
            out = None
            for j in range(3):
                if j != i:
                    factor = (h - eye()*acb(l[j]))*acb(1/(l[i]-l[j]))
                    out = factor if out is None else out*factor
            return out
        P = [proj(hs[0], lam[0], i) for i in range(3)]
        Q = [proj(hs[1], lam[1], i) for i in range(3)]
        eqs += [trace(P[i]*Q[j])*acb(1/t)-acb(1) for (i, j), t in zip(((0,1), (1,2), (0,2)), TV)]
        phase = (trace(P[0]*Q[1]*P[1]*Q[2]) + acb(kappa))*acb((-delta).cos(), (-delta).sin())
        for e in eqs:
            vals.append(e.v.real); grads.append([g.real for g in e.g])
        vals.append(phase.v.imag); grads.append([g.imag for g in phase.g])
    return vals, grads


def krawczyk(A0s):
    x = [arb(s) for s in DATA["hit"]["x"]]
    active = DATA["hit"]["active"]
    assert len(active) == 39 and len(set(active)) == 39 and all(0 <= a < 72 for a in active)
    rad = arb(DATA["hit"]["radius"])
    vals, grads = equations(A0s, x, active)
    F = arb_mat([[v] for v in vals])
    Jm = arb_mat([[arb(g.mid()) for g in row] for row in grads])
    one = arb_mat([[int(i==j) for j in range(39)] for i in range(39)])
    Y = Jm.solve(one, algorithm="approx")
    Y = arb_mat([[arb(Y[i,j].mid()) for j in range(39)] for i in range(39)])
    box, half = list(x), []
    for a in active:
        r = rad*arb(max(arb(1), abs(x[a])).upper())
        half.append(arb(0, r.upper()))
        box[a] = x[a] + half[-1]
    _, gb = equations(A0s, box, active)
    # X−x̃ 取盒减中心（含中心球的半径），使结论对中心球内任一点成立。
    K = -(Y*F) + (one - Y*arb_mat(gb))*arb_mat([[box[a]-x[a]] for a in active])
    ratio = max((abs(K[i,0])/half[i].rad()).upper() for i in range(39))
    assert all(abs(K[i,0]).upper() < half[i].rad() for i in range(39))
    assert ratio < arb("1e-4")
    print("Krawczyk：39 维盒（半宽", DATA["hit"]["radius"], "·max(1,|x_a|)）内 K(X)⊂int X，分量比最大", arb(ratio).str(3))
    return box


def kernel_check(A0s):
    """数值核对（非证明）：中心处全 39×72 雅可比作用在 4 个相位方向与 3 个伸缩方向上为零。
    满秩 39 已由 Krawczyk 判据（39×39 活动雅可比非奇异）认证；这 7 个方向在核里由不变性书面论证给出。"""
    x = [arb(s) for s in DATA["hit"]["x"]]
    vals, grads = equations(A0s, x, list(range(72)))
    directions = []
    v = [arb(0)]*72                                  # w ↦ e^{iα}w
    for t in range(27):
        v[t], v[27+t] = -x[27+t], x[t]
    directions.append(v)
    for a in range(3):
        o = 54+6*a
        ph, sc = [arb(0)]*72, [arb(0)]*72            # u_n ↦ e^{iβ}u_n 与 u_n ↦ r u_n
        for k in range(3):
            ph[o+k], ph[o+3+k] = -x[o+3+k], x[o+k]
            sc[o+k], sc[o+3+k] = x[o+k], x[o+3+k]
        directions += [ph, sc]
    for v in directions:
        scale = max(sum((abs(g[j]*v[j]) for j in range(72)), arb(0)).upper() for g in grads)
        resid = max(abs(sum((g[j]*v[j] for j in range(72)), arb(0))).upper() for g in grads)
        assert arb(resid) < arb("1e-30")*arb(scale)
    print("7 个相位与伸缩方向上全雅可比的作用相对 <1e-30（数值核对）")


def hit_checks(models, box):
    lam, TV, kappa, delta = targets()
    T = single["target_masses"]()
    target = [y for g in ("up", "down", "leptons") for y in T[g]]
    goals = [arb(a).sin() for a in (".22704", ".04275", ".00377")] + [arb("1.139")]
    w = acb_mat([[acb(box[t], box[27+t])] for t in range(27)])
    old = single["column"](PREV["full"]["w"])
    wnorm = sum((abs(w[t,0])**2 for t in range(27)), arb(0)).sqrt()
    dist = sum((abs(w[t,0]-old[t,0])**2 for t in range(27)), arb(0)).sqrt()/wnorm
    assert arb(WNORM[0]) < wnorm < arb(WNORM[1]) and dist < arb(OLD_DIST)
    print("精确共用命中 ‖w‖ ∈ (", WNORM[0], ",", WNORM[1], ")，与已有近似共用见证的相对距离 <", OLD_DIST)
    four_pi = 4*arb.pi()
    for a, n in enumerate(ZEROS):
        o = 54+6*a
        u = acb_mat([[acb(box[o+k], box[o+3+k])] for k in range(3)])
        z = u/(adjoint(u)*u)[0,0].real.sqrt()
        blocks, heavy, cols = models[n]
        matrices = study["mass_matrices"](blocks, w, z)
        masses = []
        for y in matrices:
            sv = sorted((r.real.sqrt() for r in (y*adjoint(y)).eig(algorithm="rump")), key=float)
            assert all(sv[i].upper() < sv[i+1].lower() for i in range(2))
            masses += sv
        obs = MIXING["observables"](*matrices[:2])
        angles = obs[4:]
        # 盒内的精确解使每个包围都含其精确值；方程给出质量与 CKM 模方精确等于目标，相位方程给出 δ≡1.139 (mod π)。
        # δ 取值于 (−π,π]，另一候选为 1.139−π；δ 的包围含 1.139 且宽度 <1e-18，不含 1.139−π，于是 δ 精确等于 1.139。
        assert all((y-ref).contains(0) for y, ref in zip(masses, target))
        assert all((v-ref).contains(0) for v, ref in zip(angles[:4], goals))
        assert not (angles[3] - arb("1.139") + arb.pi()).contains(0)
        width = max(abs(y/ref-1).upper() for y, ref in zip(masses, target))
        mix = max(abs(v-ref).upper() for v, ref in zip(angles[:4], goals))
        assert width < arb("1e-21") and mix < arb("1e-18")
        assert abs(angles[4] - arb("3.20630e-5")) < arb("5e-11")
        flat = sum((z[k,0]*(cols[k]*w) for k in range(3)), acb_mat(361, 1))
        Yh = acb_mat([[flat[19*s+t,0] for t in range(19)] for s in range(19)])
        full, hh, mixed = single["vertex_blocks"](n, w, z)
        assert all((hh[i,j]-Yh[i,j]).contains(0) for i in range(19) for j in range(19))
        frob2 = sum((abs(Yh[i,j])**2 for i in range(19) for j in range(19)), arb(0))
        claim = HIT[str(n)]
        lo, hi = map(arb, claim["spec"])
        single["spectral_below"](Yh, claim["spec"][1])
        assert single["rayleigh"](Yh, DATA["hit"]["top"][str(n)]) > lo
        assert arb(claim["frob"][0]) < frob2.sqrt() < arb(claim["frob"][1])
        assert lo > four_pi
        (q_lo, q_hi), (u_lo, u_hi) = single["heavy_mass_balls"](n)
        quad = single["quad"]
        h_lo, h_hi = quad(q_lo**2, u_lo**2)*frob2, quad(q_hi**2, u_hi**2)*frob2
        assert arb(claim["quad"][0]) < h_lo and h_hi < arb(claim["quad"][1])
        print(n, "九质量与 CKM 四量包围均含目标（相对宽度 <", arb(width).str(3), "），J", angles[4].str(8))
        print("  重—重块谱范数 ∈ (", claim["spec"][0], ",", claim["spec"][1], ")，Frobenius", frob2.sqrt().str(8),
              "，h² 系数 ∈ (", h_lo.str(4), ",", h_hi.str(4), ") GeV²")


def rank_check(number):
    # 重—重块左乘 left=1−v̄_q v_q^T。按模型定义（三个轻模是正交归一的本征向量）它是秩 16 的精确正交投影，
    # 故 ‖Y‖₂ ≥ ‖Y‖_F/4；此处只以 1e-40 容差作一致性核对，不由它推出精确的秩。
    fields, ops = prev["legs"](number)
    vq = fields[1][0]
    left = ops[1][0] - vq.conjugate()*vq.transpose()
    assert all(abs((left*left - left)[i,j]) < arb("1e-40") for i in range(19) for j in range(19))
    assert all(abs(left[i,j] - left[j,i].conjugate()) < arb("1e-40") for i in range(19) for j in range(19))
    tr = sum((left[i,i] for i in range(19)), acb(0))
    assert abs(tr - 16) < arb("1e-40")


def joint_bounds(models, factor):
    T = single["target_masses"]()
    S = [sum((x*x for x in T[g]), arb(0)) for g in ("up", "down", "leptons")]
    for n in ZEROS:
        rank_check(n)
    grams = {}
    for n in ZEROS:
        blocks, heavy, _ = models[n]
        grams[n] = [(heavy[k], [adjoint(c)*c for c in blocks[k]]) for k in range(3)]
    results = {}
    for key, row in DATA["joint"].items():
        zs = tuple(int(s) for s in key.split(","))
        theta = [arb(s) for s in row["theta"]]
        a = [[arb(s) for s in r] for r in row["a"]]
        assert all(t >= 0 for t in theta) and sum(theta, arb(0)) > 0
        for ks in product(range(3), repeat=len(zs)):
            M = acb_mat(27, 27)
            for i, (n, k) in enumerate(zip(zs, ks)):
                H, F = grams[n][k]
                M += theta[i]*H
                for f in range(3):
                    M -= a[i][f]*F[f]
            positive(M)
        fac = None
        for n in zs:
            fac = factor[n] if fac is None else prev["lesser"](fac, factor[n])
        total = sum(theta, arb(0))
        exact = sum((a[i][f]*S[f] for i in range(len(zs)) for f in range(3)), arb(0))
        tolerant = sum((a[i][f]*S[f]*((1-EPS)**2 if a[i][f] > 0 else (1+EPS)**2)
                        for i in range(len(zs)) for f in range(3)), arb(0))
        assert all(a[i][f] > 0 or a[i][f] < 0 for i in range(len(zs)) for f in range(3))
        f0, f1 = fac*exact/total, fac*tolerant/total
        s0, s1 = (f0/16).sqrt(), (f1/16).sqrt()
        c = JOINT[key]
        assert f0 > arb(c[0]) and f1 > arb(c[1]) and s0 > arb(c[2]) and s1 > arb(c[3])
        results[key] = (f0, f1, s0, s1)
        print("共用零点", zs, "：max‖Y_hh‖_F² >", c[0], "(精确)，>", c[1], "(ε=1e-3)；谱范数 >", c[2], "，>", c[3])
    return results


def verify():
    ctx.prec = 448
    models = {n: prev["model"](n) for n in ZEROS}
    A0s = {n: coefficient_matrix(models[n][0]) for n in ZEROS}
    # A(z)=D(z)A0 的核对：对一个非纯真空 z 比较两种组装。
    z = acb_mat([[acb("0.3")], [acb("0.5", "0.2")], [acb("0.1", "-0.7")]])
    for n in ZEROS:
        blocks = models[n][0]
        for f in range(3):
            for i in range(3):
                for j in range(3):
                    k = (-i-j) % 3 if f == 0 else (i+j) % 3
                    fac = z[k,0] if f == 0 else z[k,0].conjugate()
                    assert all((fac*A0s[n][9*f+3*i+j,t] - fac*blocks[k][f][i,t]).contains(0) for t in range(27))
    box = krawczyk(A0s)
    kernel_check(A0s)
    hit_checks(models, box)
    factor = prev["leakage"](models)
    joint = joint_bounds(models, factor)
    # 与逐零点证书（#11661 的 0.4473 为三者最大）比较；两零点情形与相应单零点下界的最大者比较。
    # 逐零点下界引自 ch39-single-zero-heavy-block-double-hop-hit.py 的证书，本脚本不复算。
    single_lo = {11661: arb(".4473"), 11934: arb(".2976"), 22027: arb(".3852")}
    for key, (f0, f1, s0, s1) in joint.items():
        zs = [int(t) for t in key.split(",")]
        best = arb(max(single_lo[n].upper() for n in zs))          # 端点比较，不比较球
        assert s1 > best
    # 精确共用命中在 N 上即 N 共用命中；上端取其各零点 Sylvester 上界（十进制数）中的最大者。
    upper = {key: arb(max(arb(HIT[str(n)]["spec"][1]).upper() for n in map(int, key.split(",")))) for key in joint}
    assert upper["11661,11934,22027"] < arb("6.6634e4") and upper["11934,22027"] < arb("6.6489e4")
    for key in joint:
        print("共用零点", key, "：精确命中上 max_n 谱范数之下确界 ∈ (", JOINT[key][2], ",", upper[key].str(7), ")；ε=1e-3 容差集上 >", JOINT[key][3])


if __name__ == "__main__":
    verify()

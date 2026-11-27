"""每腿至多两次邻边输运的二十七项字典：三枚候选零点各自（不共用系数）存在精确完整命中，其重—重上型块谱范数低于 4π。

区间证书，448 位。每枚零点由数据文件给出一个真空 z 与五个酉矩阵的有限小数参数。脚本以精确目标质量与精确 CKM 参数
组装三部门质量矩阵 M，认证 27×27 线性方程组 A(z)w=vec(M) 可逆并求出解的包围球 w*；该方程组的唯一精确解即一个精确
完整命中。随后在 w* 球上认证重—重上型块谱范数的上下界、h² 系数的上下界，复算同一零点上任何命中被迫满足的谱范数下界，
并认证完整上型顶角与重双重态×轻单态混合块的谱范数都大于 4π；另认证已有跨零点共用见证的重块谱范数上界。
"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

HERE = Path(__file__).parent
prev = run_path(str(HERE/"ch39-per-leg-double-hop-transfer.py"))
study, adjoint, positive = prev["study"], prev["adjoint"], prev["positive"]
REFS = run_path(str(HERE/"ch39-shell-flavor-matching-windows.py"))
MIXING = run_path(str(HERE/"ch39-finite-orthogonal-ckm.py"))
DATA = json.loads((HERE/"ch39-single-zero-heavy-block-double-hop-hit.json").read_text())
PREV = json.loads((HERE/"ch39-per-leg-double-hop-transfer.json").read_text())
ZEROS = (11661, 11934, 22027)
MONOS = prev["MONOS"]
RD, RE = arb("61.36519"), arb("84.98147")

# 见证重—重块谱范数的 Rayleigh 下界与 Sylvester 上界、Frobenius 区间、h² 系数（GeV²）的上下界。
CLAIMS = {
    "11661": dict(spec=("4.4318", "4.4319"), frob=("9.7526", "9.7527"), quad=("7.524e32", "7.526e32")),
    "11934": dict(spec=("4.4287", "4.4289"), frob=("9.6185", "9.6186"), quad=("7.378e32", "7.380e32")),
    "22027": dict(spec=("4.1573", "4.1574"), frob=("9.6856", "9.6857"), quad=("7.319e32", "7.321e32")),
}
FORCED = {"11661": ".4473", "11934": ".2976", "22027": ".3852"}       # 同一零点上任何命中的谱范数下界（复算）
RATIO = {"11661": "9.91", "11934": "14.89", "22027": "10.80"}          # 见证上界 / 强制下界
FULL_LO = {"11661": "53.8", "11934": "59.0", "22027": "56.0"}          # 完整上型顶角谱范数下界
MIXED_LO = {"11661": "53.7", "11934": "58.9", "22027": "55.8"}         # 重双重态×轻单态块谱范数下界
SHARED_HI = {"11661": "6.6634e4", "11934": "6.6489e4", "22027": "6.5984e4"}   # 已有共用见证重块谱范数上界


def column(values):
    return acb_mat([[acb(arb(x), arb(y))] for x, y in values])


def unitary(v):
    # U=exp(iH)，H 由九个有限小数给出：对角 v[0:3]，上三角 (0,1),(0,2),(1,2) 为 v[3:6]+i v[6:9]。
    v = [arb(s) for s in v]
    H = acb_mat(3, 3)
    for i in range(3):
        H[i,i] = acb(v[i])
    for t, (i, j) in enumerate(((0,1), (0,2), (1,2))):
        H[i,j] = acb(v[3+t], v[6+t])
        H[j,i] = acb(v[3+t], -v[6+t])
    return (H*acb(0, 1)).exp()


TARGET_STRINGS = {"up": (".00000615", ".00311", ".8616"), "down": (".0000135", ".000268", ".01401"),
                  "leptons": (".0000028227", ".00059465", ".010101")}   # 与匹配窗脚本的中心值相同（脚本断言）


def target_masses():
    out = {g: [arb(x) for x in v] for g, v in TARGET_STRINGS.items()}
    for g in out:
        assert all(abs(a-b) < arb("1e-70") for a, b in zip(out[g], REFS[g]))
    return out


def targets():
    T = target_masses()
    dg = lambda m: acb_mat([[m[i] if i == j else 0 for j in range(3)] for i in range(3)])
    s12, s23, s13 = [arb(a).sin() for a in (".22704", ".04275", ".00377")]
    c12, c23, c13 = [(1-s*s).sqrt() for s in (s12, s23, s13)]
    e = (acb(0, -1)*arb("1.139")).exp()
    ckm = acb_mat([[c12*c13, s12*c13, s13*e],
                   [-s12*c23-c12*s23*s13/e, c12*c23-s12*s23*s13/e, s23*c13],
                   [s12*s23-c12*c23*s13/e, -c12*s23-s12*c23*s13/e, c23*c13]])
    return dg(T["up"]), dg(T["down"]), dg(T["leptons"]), ckm


def exact_hit(number, row, blocks):
    # 精确目标：M_u=VΣ_uW†，M_d=V V_CKM Σ_d X†，M_e=PΣ_eQ†；解 A(z)w=vec(M)。
    z = column([row["z"][k] for k in range(3)])
    z = z/(adjoint(z)*z)[0,0].real.sqrt()
    V, W, X, P, Q = [unitary(row["unitary"][9*i:9*i+9]) for i in range(5)]
    su, sd, se, ckm = targets()
    mats = (V*su*adjoint(W), V*ckm*sd*adjoint(X), P*se*adjoint(Q))
    rows, rhs = [], []
    for f, M in enumerate(mats):
        for i in range(3):
            for j in range(3):
                k = (-i-j) % 3 if f == 0 else (i+j) % 3
                fac = z[k,0] if f == 0 else z[k,0].conjugate()
                rows.append([fac*blocks[k][f][i,t] for t in range(27)])
                rhs.append([M[i,j]])
    w = acb_mat(rows).solve(acb_mat(rhs))      # 失败即抛异常；成功即认证可逆并包围唯一解
    radius = max(max(w[t,0].real.rad(), w[t,0].imag.rad()) for t in range(27))
    assert radius < arb("1e-15"), radius
    return w, z


def vertex_blocks(number, w, z):
    # 未投影的上型顶角 Y_full=Σ_k z_k Σ_t w_t (T^a)^T diag(h_c,k) T^b，以及它的重—重块与重双重态×轻单态块。
    fields, ops = prev["legs"](number)
    h, vq, vu = fields[3], fields[1][0], fields[-4][0]
    one = ops[1][0]
    left = one - vq.conjugate()*vq.transpose()
    right_heavy = one - vu*adjoint(vu)
    right_light = vu*adjoint(vu)
    full = acb_mat(19, 19)
    for k in range(3):
        for t, (a, b, c) in enumerate(MONOS):
            col = ops[1][a].transpose()*acb_mat([[h[c][x,k]*ops[-4][b][x,j] for j in range(19)] for x in range(19)])
            full += (z[k,0]*w[t,0])*col
    return full, left*full*right_heavy, left*full*right_light


def rayleigh(Y, v):
    v = column(v)
    gram = adjoint(Y)*Y
    return ((adjoint(v)*gram*v)[0,0].real/(adjoint(v)*v)[0,0].real).sqrt()


def spectral_below(Y, t):
    one = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    positive(arb(t)**2*one - adjoint(Y)*Y)


def heavy_mass_balls(number):
    shell = prev["shell"]
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19)) / 18
    out = []
    for q in (1, -4):
        levels = [abs(r.real*(r.real-center)) for r in shell["laplacian"](q, height).eig(algorithm="rump")]
        heavy = [v for v in levels if v > 9]
        assert len(heavy) == 16 and all(v > 9 or v < 9 for v in levels)
        lo, hi = min(v.lower() for v in heavy), min(v.upper() for v in heavy)
        scale = arb("2.435e18")/724**2*(32 if q == -4 else 1)
        out.append((arb(lo)*scale, arb(hi)*scale))
    return out


def quad(q2, u2):
    return q2 + u2 - q2*u2/(q2+u2)


def forced_lower(number, models, factor):
    cap, mu = map(arb, PREV["caps"]["all"][str(number)])
    assert cap > 0 and mu >= 0
    blocks, heavy, _ = models[number]
    for k in range(3):
        u, d, l = [adjoint(c)*c for c in blocks[k]]
        positive(cap*heavy[k] - u + mu*((RD**2*d + RE**2*l)/2 - u))
    return (factor[number]*arb(".86")**2/cap/16).sqrt()


def shared_witness_upper(models):
    # 已有跨零点共用见证（ch39-per-leg-double-hop-transfer.json）的重—重块谱范数上界。
    w = column(PREV["full"]["w"])
    for number in ZEROS:
        z = column(PREV["full"]["higgs"][str(number)])
        z = z/(adjoint(z)*z)[0,0].real.sqrt()
        cols = models[number][2]
        flat = sum((z[k,0]*(cols[k]*w) for k in range(3)), acb_mat(361, 1))
        Y = acb_mat([[flat[19*s+t,0] for t in range(19)] for s in range(19)])
        spectral_below(Y, SHARED_HI[str(number)])
        lo = rayleigh(Y, PREV["full"]["top"][str(number)])
        assert lo > arb(SHARED_HI[str(number)]) - 1
        print(number, "已有共用见证重块谱范数 ∈ (", lo.str(8), ",", SHARED_HI[str(number)], ")")


def verify():
    ctx.prec = 448
    four_pi = 4*arb.pi()
    T = target_masses()
    target = [y for g in ("up", "down", "leptons") for y in T[g]]
    goals = [arb(x).sin() for x in (".22704", ".04275", ".00377")] + [arb("1.139")]
    assert T["up"][-1] > arb(".86")
    # 逐零点下界证书的两条前提对容差集 H(ε)（各质量相对误差 ≤ε）在 ε=1e-3 时仍成立。
    S = {g: sum((x*x for x in T[g]), arb(0)) for g in T}
    gap = S["up"] - (RD**2*S["down"] + RE**2*S["leptons"])/2
    assert arb(".0029636") < gap < arb(".0029638")
    eps = arb("1e-3")
    assert S["up"]*(1-eps)**2 - (RD**2*S["down"] + RE**2*S["leptons"])/2*(1+eps)**2 > 0
    assert T["up"][-1]*(1-eps) > arb(".86")
    models = {n: prev["model"](n) for n in ZEROS}
    factor = prev["leakage"](models)
    shared_witness_upper(models)
    solutions = {}
    for number in ZEROS:
        key = str(number)
        row, claim = DATA[key], CLAIMS[key]
        blocks, heavy, cols = models[number]
        w, z = exact_hit(number, row, blocks)
        solutions[number] = w
        # 在 w* 球上复核命中：质量、CKM、J 的包围与目标一致（精确性由构造给出，见说明）。
        matrices = study["mass_matrices"](blocks, w, z)
        masses = []
        for y in matrices:
            sv = sorted((r.real.sqrt() for r in (y*adjoint(y)).eig(algorithm="rump")), key=float)
            assert all(sv[i].upper() < sv[i+1].lower() for i in range(2))      # 三个奇异值球互不重叠
            masses += sv
        error = max(abs(y/ref-1).upper() for y, ref in zip(masses, target))
        angles = MIXING["observables"](*matrices[:2])[4:]
        mix = max(abs(v-ref).upper() for v, ref in zip(angles[:4], goals))
        # 精确解 w* 在球内，故每个包围球都含精确目标；球宽只是区间放大，不是命中误差。
        assert all((y-ref).contains(0) for y, ref in zip(masses, target))
        assert all((v-ref).contains(0) for v, ref in zip(angles[:4], goals))
        assert error < arb("1e-12") and mix < arb("1e-10")
        assert abs(angles[4] - arb("3.20630e-5")) < arb("5e-11")
        # 重—重上型块：与 ch39-per-leg-double-hop-transfer.py 的列组装一致，再与未投影组装交叉核对。
        flat = sum((z[k,0]*(cols[k]*w) for k in range(3)), acb_mat(361, 1))
        Y = acb_mat([[flat[19*s+t,0] for t in range(19)] for s in range(19)])
        full, hh, mixed = vertex_blocks(number, w, z)
        assert all((hh[i,j]-Y[i,j]).contains(0) for i in range(19) for j in range(19))
        frob2 = sum((abs(Y[i,j])**2 for i in range(19) for j in range(19)), arb(0))
        lo, hi = map(arb, claim["spec"])
        spectral_below(Y, claim["spec"][1])
        assert rayleigh(Y, row["top_heavy"]) > lo
        f_lo, f_hi = map(arb, claim["frob"])
        assert f_lo < frob2.sqrt() < f_hi
        assert hi < four_pi and hi < arb("16.63") and lo > four_pi.sqrt()
        # 完整上型顶角与重双重态×轻单态块：谱范数下界均大于 4π。
        r_full, r_mixed = rayleigh(full, row["top_full"]), rayleigh(mixed, row["top_mixed"])
        assert r_full > arb(FULL_LO[key]) and r_mixed > arb(MIXED_LO[key]) and min(arb(FULL_LO[key]), arb(MIXED_LO[key])) > four_pi
        wnorm = sum((abs(w[t,0])**2 for t in range(27)), arb(0)).sqrt()
        w_lo, w_hi = {"11661": (5225, 5226), "11934": (4201, 4202), "22027": (3833, 3834)}[key]
        assert w_lo < wnorm < w_hi
        # h² 系数：f(a,b)=a+b-ab/(a+b) 对 a=m_Q²、b=m_U² 均单调增（∂f/∂a=1-(b/(a+b))²>0，对称），且 b<f<b(1+a²/b²)。
        (q_lo, q_hi), (u_lo, u_hi) = heavy_mass_balls(number)
        assert arb("5.08e13") < q_lo and q_hi < arb("5.10e13") and arb("2.79e15") < u_lo and u_hi < arb("2.83e15")
        h_lo, h_hi = quad(q_lo**2, u_lo**2)*frob2, quad(q_hi**2, u_hi**2)*frob2
        a, b = q_lo**2, u_lo**2
        assert b < quad(a, b) < b*(1 + a**2/b**2)
        qa, qb = map(arb, claim["quad"])
        assert qa < h_lo and h_hi < qb and h_lo/arb("125.25")**2 > arb("4.6e28")
        forced = forced_lower(number, models, factor)
        assert forced > arb(FORCED[key]) and hi/forced < arb(RATIO[key])
        print(number, "w* 包围半径 <1e-15；九质量包围含目标、相对宽度 <", error.str(3), "；CKM 四量包围含目标、宽度 <", mix.str(3), "J", angles[4].str(8))
        print("  重—重块 Frobenius", frob2.sqrt().str(8), "谱范数 ∈ (", lo.str(6), ",", hi.str(6), ")；h² 系数 ∈ (",
              h_lo.str(6), ",", h_hi.str(6), ") GeV²；‖w‖", wnorm.str(5))
        print("  完整上型顶角谱范数 >", r_full.str(5), "，重双重态×轻单态块谱范数 >", r_mixed.str(5), "，4π =", four_pi.str(6))
        print("  同一零点下确界区间 (", forced.str(5), ",", hi.str(6), ")，上界/下界 <", RATIO[key])
    # 三组系数两两不同：在第 0 个系数上球不含 0。
    for i in range(3):
        for j in range(i+1, 3):
            d = solutions[ZEROS[i]] - solutions[ZEROS[j]]
            assert any(not d[t,0].contains(0) for t in range(27))


if __name__ == "__main__":
    verify()

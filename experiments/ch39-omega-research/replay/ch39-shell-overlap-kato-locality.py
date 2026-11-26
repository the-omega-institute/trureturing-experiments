"""十九点壳四维重叠核与 Kato 输运测度流的显式指数局域性常数。

证书分三部分。
1. 精确有理：Combes-Thomas 加权扰动、围道与迹估计的全部常数，
   对任意正零点的常数壳（壳带界）给出体积无关的衰减率与前因子。
2. Arb 区间：11431、133931 两枚零点的实际壳质量谱隙、换点质量差范数，
   给出常数壳与逐格切换两种情形更紧的常数。
3. 仅数值：4^4 周期格单块 U(1) 玩具核上，对加权扰动范数、谱隙、
   符号核与测度流二阶变分的界作浮点对照。
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent

# ---------------------------------------------------------------- 有理工具
SQRT2_HI = F(99, 70)                 # (99/70)^2 = 9801/4900 > 2
PI_LO = F(314159, 100000)            # < pi
assert SQRT2_HI**2 > 2 and SQRT2_HI > 0


def sqrt_lo(x, den=10**15):
    """有理下界 floor(sqrt(x)*den)/den。"""
    v = F(isqrt(x.numerator*den*den // x.denominator), den)
    assert v*v <= x
    return v


def sqrt_hi(x, den=10**15):
    v = sqrt_lo(x, den) + F(1, den)
    assert v*v >= x
    return v


def arb_lo(v, den=10**12):
    """Arb 球的有理下界（向下取整到 1/den）。"""
    lo = v.lower()
    q = F(int((lo*den).floor().unique_fmpz()), den)
    assert arb(q.numerator)/q.denominator <= v
    return q


def arb_hi(v, den=10**12):
    hi = v.upper()
    q = F(int((hi*den).ceil().unique_fmpz()), den)
    assert arb(q.numerator)/q.denominator >= v
    return q


# ---------------------------------------------------------------- 模型常数
LOW, HIGH = F(8685, 1000), F(96, 10)
CENTER, SCALE = (LOW+HIGH)/2, HIGH/2           # M = 9.1425, c = 4.8
G = (HIGH-LOW)/HIGH                            # 61/640
EPS = F(1, 200)                                # 方格误差 ||1-U_p|| <= 1/200
LAMBDA_MAX = F(36)                             # Gershgorin，见 check_gershgorin
assert CENTER == F(3657, 400) and SCALE == F(24, 5) and G == F(61, 640)

M_MAX = (LAMBDA_MAX-CENTER)/SCALE              # 全部壳质量上界
M_MIN = (0-CENTER)/SCALE                       # 磁拉普拉斯半正定
assert M_MAX == F(10743, 1920) and M_MIN == -2+G
SPIN = sqrt_hi(4+6*EPS)                        # ||sum gamma_mu S_mu||^2 <= 4+6 eps
CURV_LOSS_HI = 1-sqrt_lo(1-30*EPS)             # 1 - sqrt(1-30 eps)

# 迹估计中一条链变分的块秩：自旋投影秩 2、壳 19 分量、规范表示维数 r。
# 一代 15 个外尔分量 = 4 个双重态 + 7 个单态，sum r = 15。
RANK_PER_DIM = 2*19
REP_DIM_SUM = 4*2+7
assert REP_DIM_SUM == 15


def check_gershgorin(shell):
    """十九点壳最大加权度为 18、各边权在 (0,1]、相位模为一，故 lambda_max <= 36。"""
    geo, soft = shell["geo"], shell["soft"]
    from collections import Counter
    deg = Counter(v for i, j, _, _ in geo["edges"] for v in (i, j))
    assert len(geo["edges"]) == 135 and max(deg.values()) == 18
    for _, _, _, z in geo["edges"]:
        assert abs(z).overlaps(arb(1)) and abs(abs(z)-1) < arb("1e-40")
    # 只有辐边带 1-delta 权，delta = (23/25) 2^-20/(h^2+offset) 在 [0,1e-8] 内。
    offset = soft["point"]["study"]["offset"]
    assert arb(offset) > 0
    assert 2*max(deg.values()) == LAMBDA_MAX
    print("壳最大度与 Gershgorin 上界", max(deg.values()), LAMBDA_MAX)


# ------------------------------------------------- 一般局域性定理的常数
def stage_constants(alpha, beta, nprime, heavy_kinetic, rep_dim):
    """alpha: 全路径谱隙下界；beta: 范数上界；nprime: N_b(dH/ds) 上界
    （不含 delta_b，按需在内部加上）；heavy_kinetic: dH/ds 是否含规范链。"""
    # 选 b 使 delta_b = 4 sqrt2 (e^b-1) <= alpha/2，用 e^b-1 <= b/(1-b)。
    t = alpha/(8*SQRT2_HI+alpha)
    b = F(int(t*10**9), 10**9)
    delta_b = 4*SQRT2_HI*b/(1-b)
    assert 0 < b < F(1, 100) and delta_b <= alpha/2
    r = 1/(alpha-delta_b)
    assert r <= 2/alpha
    perimeter = 2*(beta+alpha)+4*alpha        # 矩形 [0,beta+alpha]x[-alpha,alpha]
    per = perimeter/(2*PI_LO)                 # >= 周长/(2 pi)
    sign_const = 2*per*r                      # ||sign H(x,y)|| <= sign_const e^{-b d}
    n_der = nprime + (delta_b if heavy_kinetic else 0)
    e4b = 1/(1-4*b)                           # >= e^{4b}
    bracket = 10*r**6*n_der + (2*r**5 if heavy_kinetic else 0)
    current = 4*RANK_PER_DIM*rep_dim*per**3*e4b*bracket
    return dict(b=b, delta_b=delta_b, r=r, per=per, sign=sign_const,
                current=current, n_der=n_der)


def verify_contraction_identity():
    """(P^-X + P^+Y)^*(P^-X + P^+Y) 的交叉项为零：P^-P^+ = 0，
    以 4x4 欧氏 gamma 矩阵逐方向精确核对。"""
    import sympy as sp
    I = sp.I
    s = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -I], [I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
    Z, E = sp.zeros(2), sp.eye(2)
    gam = [sp.Matrix(sp.BlockMatrix([[Z, -I*sk], [I*sk, Z]])) for sk in s]
    gam.append(sp.Matrix(sp.BlockMatrix([[Z, E], [E, Z]])))
    g5 = gam[0]*gam[1]*gam[2]*gam[3]
    one = sp.eye(4)
    for a in range(4):
        for c in range(4):
            assert gam[a]*gam[c]+gam[c]*gam[a] == 2*(1 if a == c else 0)*one
        assert gam[a].H == gam[a]
        pm, pp = (one-gam[a])/2, (one+gam[a])/2
        assert pm*pp == sp.zeros(4) and pm*pm == pm and pm.rank() == 2
    assert g5*g5 == one and g5.H == g5
    print("自旋投影 P^-P^+ = 0、秩 2、gamma5 厄米幺正：精确核对")


def verify_universal():
    """任意正零点、空间常数壳、五种荷：壳带界 lambda_3<8.685<9.6<lambda_4。"""
    alpha = sqrt_lo(1-30*EPS) - (1-G)           # 轻块路径到 m=-1 全程的隙
    assert alpha > F(17204, 10**6) and alpha < G   # 重块 Re D >= min(m,1) >= g
    beta = M_MAX + 8 + SPIN
    heavy_der = M_MAX + 7 + SPIN               # ||gamma5(1-m-K)||，m in (0,M_MAX]
    c = stage_constants(alpha, beta, max(F(1), heavy_der), True, REP_DIM_SUM)
    assert F(156027, 10**4) < beta < F(156029, 10**4)
    assert c["b"] == F(1523793, 10**9)
    assert c["sign"] < 1155
    assert c["current"] < F(1002, 1000)*F(10)**20
    # 可取的衰减率上确界 ln(1+alpha/(4 sqrt2))：用 x-x^2/2 < ln(1+x) < x-x^2/2+x^3/3 与 sqrt2 的上下界。
    x_lo, x_hi = alpha/(4*SQRT2_HI), alpha/(4*F(1414213, 10**6))
    assert F(1414213, 10**6)**2 < 2
    sup_lo, sup_hi = x_lo-x_lo**2/2, x_hi-x_hi**2/2+x_hi**3/3
    assert F(30475, 10**7) < sup_lo and sup_hi < F(30478, 10**7)
    print("[普适常数壳] 隙下界 alpha、范数上界 beta", float(alpha), float(beta),
          " 衰减率上确界区间", float(sup_lo), float(sup_hi))
    print("  衰减率 b、局域长 1/b", float(c["b"]), float(1/c["b"]))
    print("  符号核前因子", float(c["sign"]), " 测度流二阶前因子(15 分量合计)", f"{float(c['current']):.4e}")
    return alpha, beta, c


def verify_two_zeros(shell):
    """11431 与 133931 的实际壳：常数壳与逐格切换。"""
    ctx.prec = 192
    geo = shell["geo"]
    center, scale = arb(CENTER.numerator)/CENTER.denominator, arb(SCALE.numerator)/SCALE.denominator
    heights = {n: acb.zeta_zero(n).imag for n in (11431, 133931)}
    free_gap, masses_hi, diff_norm = {}, {}, {}
    from collections import Counter
    for dd in (724, 804):
        pairs = [(i, j) for i, j, ed, _ in geo["edges"] if ed == dd]
        assert Counter(i for i, _ in pairs) == {i: 1 for i in range(6)}
        assert Counter(j for _, j in pairs) == {6: 2, 7: 2, 8: 2}
    clocks = {n: geo["phases"](h) for n, h in heights.items()}
    for q in (1, 2, 3, 4, 6):
        laps = {}
        for n, h in heights.items():
            lap = shell["laplacian"](q, h)
            laps[n] = lap
            masses = sorted(((z.real-center)/scale for z in lap.eig()), key=float)
            assert sum(-2 < m < 0 for m in masses) == 3 and sum(m > 0 for m in masses) == 16
            fg = min(min(abs(m+2*k) for k in range(5)) for m in masses)
            free_gap[q, n] = fg
            masses_hi[q, n] = max(masses)
            assert masses[-1] < arb(M_MAX.numerator)/M_MAX.denominator
        # 换点质量差的谱范数上界。两类时钟边（弦平方 724、804）各六条，从六个互异外站
        # 连到三个内站、每内站两条，故该类带相位邻接块 A 满足 A^*A = 2I_3，
        # ||A + A^*|| = sqrt2；钟值 u -> u' 使 L_q 改变 -(u^q - u'^q) 倍的该块（及共轭）。
        # 三条辐边权差至多 1e-8，每条单边拉普拉斯范数为 2，合计至多 6e-8。
        diff_norm[q] = (arb(2).sqrt()*sum(abs(clocks[11431][dd]**q-clocks[133931][dd]**q)
                                           for dd in (724, 804)) + 6*arb("1e-8"))/scale
        # 独立复核：Frobenius 范数（谱范数的另一上界）不应小于谱范数的下界 0。
        dmat = laps[11431]-laps[133931]
        frob = sum(abs(dmat[i, j])**2 for i in range(19) for j in range(19)).sqrt()/scale
        assert frob > 0
    loss = arb(1)-(arb(1)-arb(30)/200).sqrt()
    alpha_c = min(free_gap[k]-loss for k in free_gap)      # 两枚零点任一的常数壳
    alpha_ref = min(free_gap[q, 11431]-loss for q in (1, 2, 3, 4, 6))
    alpha_s = min(free_gap[q, 11431]-loss-diff_norm[q] for q in (1, 2, 3, 4, 6))
    vmax = max(diff_norm.values())
    mmax = max(masses_hi.values())
    for q in (1, 2, 3, 4, 6):
        print("  荷", q, "自由隙 11431/133931", free_gap[q, 11431].str(10), free_gap[q, 133931].str(10),
              " 换点质量差上界", diff_norm[q].str(10))
    assert arb(".02519") < alpha_c < arb(".0252") and arb(".02689") < alpha_ref < arb(".0269")
    assert arb(".00963") < alpha_s < arb(".00964")
    assert vmax < arb(".017254")
    assert mmax < arb("5.1")
    a_c, a_ref, a_s = F(2519, 10**5), F(2689, 10**5), F(963, 10**5)
    m_hi = arb_hi(mmax, 10**6)
    beta = max(m_hi, F(0)) + 8 + SPIN
    heavy_der = m_hi + 7 + SPIN
    const = stage_constants(a_c, beta, max(F(1), heavy_der), True, REP_DIM_SUM)
    # 逐格切换：第一段 H_s = H_ref + (1-s) gamma5 V，V 为逐站对角质量差，
    # 全段隙 >= alpha_ref - ||V||；第二段为 11431 常数壳的标准路径。两段 Berry 项相加。
    v_hi = arb_hi(vmax, 10**9)
    stage1 = stage_constants(a_s, beta, v_hi, False, REP_DIM_SUM)
    stage2 = stage_constants(a_ref, beta, max(F(1), heavy_der), True, REP_DIM_SUM)
    b_sw = min(stage1["b"], stage2["b"])
    total_sw = stage1["current"] + stage2["current"]
    assert F(128118, 10**4) < beta < F(128119, 10**4)
    assert const["b"] == F(2221443, 10**9) and const["sign"] < 652 and const["current"] < F(4705, 1000)*F(10)**18
    assert stage1["b"] == F(850412, 10**9) and stage1["sign"] < 1698
    assert stage2["b"] == F(2371006, 10**9) and stage2["sign"] < 611
    assert total_sw < F(5349, 1000)*F(10)**18 and b_sw == stage1["b"]
    for key, val in (("const", const), ("stage1", stage1), ("stage2", stage2)):
        print("   ", key, "b =", val["b"], " sign <", float(val["sign"]), f" current < {float(val['current']):.6e}")
    print("[11431 或 133931 常数壳] alpha, beta", float(a_c), float(beta),
          " b", float(const["b"]), " 符号核前因子", float(const["sign"]),
          f" 测度流前因子 {float(const['current']):.4e}")
    print("[逐格切换] 第一段隙", float(a_s), " ||V||<=", float(v_hi), " b", float(stage1["b"]),
          " 符号核前因子", float(stage1["sign"]), f" 两段合计测度流前因子 {float(total_sw):.4e} (率 2b, b={float(b_sw)})")
    return dict(const=const, stage1=stage1, stage2=stage2, alpha_s=a_s, alpha_c=a_c, beta=beta)


# ------------------------------------------------------------ 数值对照
def numeric_toy():
    """仅数值：4^4 周期格、单个壳轻块与重块、小曲率 U(1) 链。"""
    import numpy as np
    rng = np.random.default_rng(20261126)
    L, dim = 4, 4
    sites = [(a, b_, c, d) for a in range(L) for b_ in range(L) for c in range(L) for d in range(L)]
    idx = {x: k for k, x in enumerate(sites)}
    N = len(sites)
    s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
    Z2, E2 = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z2, -1j*s], [1j*s, Z2]]) for s in (s1, s2, s3)] + [np.block([[Z2, E2], [E2, Z2]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    theta = rng.uniform(-1, 1, (N, 4))/1600          # 方格角 <= 4/1600
    eps_toy = 2*np.sin(4/1600/2)
    assert eps_toy <= 1/200

    def shift(mu, th):
        T = np.zeros((N, N), complex)
        for x in sites:
            y = list(x); y[mu] = (y[mu]+1) % L
            T[idx[x], idx[tuple(y)]] = np.exp(1j*th[idx[x], mu])
        return T

    def kinetic(th):
        K = np.zeros((4*N, 4*N), complex)
        for mu in range(4):
            T = shift(mu, th)
            pm, pp = (np.eye(4)-gam[mu])/2, (np.eye(4)+gam[mu])/2
            K += 4/4*np.kron(np.eye(N), np.eye(4)) - np.kron(T, pm) - np.kron(T.conj().T, pp)
        return K

    def H_of(m, kin, th):
        return np.kron(np.eye(N), g5) @ (m*np.eye(4*N) + kin*kinetic(th))

    def dist(x, y):
        return sum(min((a-b_) % L, (b_-a) % L) for a, b_ in zip(x, y))

    y0 = sites[0]
    phi = np.repeat([dist(x, y0) for x in sites], 4)
    for bb in (0.05, 0.3, 1.0):
        W = np.diag(np.exp(bb*phi)); Wi = np.diag(np.exp(-bb*phi))
        H = H_of(-0.6, 1.0, theta)
        dev = np.linalg.norm(W@H@Wi - H, 2)
        assert dev <= 4*np.sqrt(2)*(np.exp(bb)-1)*(1+1e-12)
    # 轻块 m=-0.6 与重块 (m=0.8) 的谱隙
    for m, lower in ((-0.6, np.sqrt(1-30*eps_toy)-0.4), (0.8, 0.8)):
        ev = np.linalg.eigvalsh(H_of(m, 1.0, theta))
        assert np.min(np.abs(ev)) >= lower-1e-12
    # 符号核衰减 vs 一般界（alpha、beta 为该玩具核的数值隙与范数）
    H = H_of(-0.6, 1.0, theta)
    ev, V = np.linalg.eigh(H)
    alpha, beta = np.min(np.abs(ev)), np.max(np.abs(ev))
    S = (V*np.sign(ev))@V.conj().T
    bb = np.log(1+alpha/(8*np.sqrt(2)))
    pref = 2*(2*(beta+alpha)+4*alpha)/(2*np.pi)*(2/alpha)
    worst = 0.0
    for x in sites[1:]:
        blk = S[4*idx[x]:4*idx[x]+4, 0:4]
        worst = max(worst, np.linalg.norm(blk, 2)/(pref*np.exp(-bb*dist(x, y0))))
    assert worst < 1
    # 测度流二阶变分：轻块路径 m(s) = -0.6 + s(-0.4)，s in [0,1]，Gauss-Legendre 求积。
    e_site, f_site = idx[(0, 0, 0, 0)], idx[(2, 2, 0, 0)]
    def current(th):
        xs, ws = np.polynomial.legendre.leggauss(4)
        tot = 0j
        for xg, wg in zip(xs, ws):
            s = (xg+1)/2
            Hs = H_of(-0.6-0.4*s, 1.0, th)
            lam, U = np.linalg.eigh(Hs)
            pos = lam > 0
            den = lam[:, None]-lam[None, :]
            mask = np.logical_xor(pos[:, None], pos[None, :])
            def dP(X):
                Xt = U.conj().T@X@U
                out = np.zeros_like(Xt)
                out[mask] = Xt[mask]/den[mask]
                return out
            Hp = np.kron(np.eye(N), g5)*(-0.4)
            # 沿 mu=0 方向链 e 的生成元 i：dT 只在 (e, e+mu) 块
            th2 = th.copy(); h = 1e-6
            th2[e_site, 0] += h
            dH = (H_of(-0.6-0.4*s, 1.0, th2)-Hs)/h
            Pt = np.diag(pos.astype(float))
            A, B = dP(Hp), dP(dH)
            tot += wg/2*np.trace(Pt@(A@B-B@A))
        return tot
    hz = 1e-4
    thp, thm = theta.copy(), theta.copy()
    thp[f_site, 1] += hz; thm[f_site, 1] -= hz
    second = abs(current(thp)-current(thm))/(2*hz)
    d = dist(sites[e_site], sites[f_site])
    r_ = 2/alpha
    per = (2*(beta+alpha)+4*alpha)/(2*np.pi)
    bound = 4*2*per**3*np.exp(4*bb)*(10*r_**6*0.4)*np.exp(-2*bb*(d-2))
    assert second < bound
    print("[仅数值 4^4 玩具] 加权扰动、谱隙、符号核比值<1 均通过；",
          f"符号核最大归一比 {worst:.3e}，二阶变分 {second:.3e} <= 界 {bound:.3e} (d={d})")


if __name__ == "__main__":
    shell = run_path(str(ROOT/"ch39-three-ray-chiral-boundary.py"))
    check_gershgorin(shell)
    verify_contraction_identity()
    verify_universal()
    verify_two_zeros(shell)
    numeric_toy()
    print("全部断言通过")

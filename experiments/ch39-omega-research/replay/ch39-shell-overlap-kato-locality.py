"""十九点壳四维重叠核与 Kato 输运项的显式指数局域性常数。

脚本只复算书面推导给出的常数公式（项数 10 与 2、块秩 4*38*r、l^3、e^{4b}
等由书面证明得到，在 stage_constants 中按公式写入），不复算推导本身。
1. 精确有理：普适常数壳的全部常数，以及界开始小于一的距离门槛。
2. Arb 区间：11431、133931 两枚零点的壳谱隙、最大壳质量、换点质量差
   （结构上界并与 Arb 实算谱范数比较）。
3. 仅数值：4^4 周期格上对两条引理（加权扰动范数、Combes-Thomas 预解式）
   的浮点检验，取一般随机链与非平凡质量；它检验引理而非最终界。
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
    # 边权：laplacian() 中只有接极点且另一端在调权轨道上的辐边取 1-delta，其余为 1。
    pole, tuned = soft["pole"], soft["star"]["orbit"]
    reweighted = [(i, j) for i, j, _, _ in geo["edges"]
                  if pole in (i, j) and (i in tuned or j in tuned)]
    assert len(reweighted) == 3
    # delta = (23/25) 2^-20/(h^2+offset)，offset>0；对 h >= 14（第一枚零点 14.13...）有 0 < delta <= 1e-8。
    offset = arb(soft["point"]["study"]["offset"])
    assert offset > 0
    delta_max = arb(23)/25*arb(2)**-20/(arb(14)**2+offset)
    assert 0 < delta_max < arb("3.9e-9")
    print("调权辐边数、权 1-delta 中 delta 上界", len(reweighted), delta_max.str(5))
    assert 2*max(deg.values()) == LAMBDA_MAX
    print("壳最大度与 Gershgorin 上界", max(deg.values()), LAMBDA_MAX)


# ------------------------------------------------- 一般局域性定理的常数
def stage_constants(alpha, beta, nprime, heavy_kinetic, rep_dim):
    """按书面定理 A、B 的公式计算常数。
    alpha: 全路径谱隙下界；beta: 范数上界；nprime: dH/ds 的范数上界
    （heavy_kinetic 时内部再加 delta_b 得加权范数）；rep_dim: sum_c r_c。
    测度流前因子对应 ||eta||, ||zeta|| 取十五个分量表示中像范数的最大值。"""
    t = alpha/(8*SQRT2_HI+alpha)
    b = F(int(t*10**9), 10**9)
    delta_b = 4*SQRT2_HI*b/(1-b)               # >= 4 sqrt2 (e^b-1)，因 e^b-1 <= b/(1-b)
    assert 0 < b < F(1, 100) and delta_b <= alpha/2
    r = 1/(alpha-delta_b)
    assert r <= 2/alpha
    perimeter = 2*(beta+alpha)+4*alpha        # 矩形 [0,beta+alpha]x[-alpha,alpha]
    per = perimeter/(2*PI_LO)                 # >= 周长/(2 pi)
    sign_const = 2*per*r                      # ||sign H(x,y)|| <= sign_const e^{-b d}
    n_der = nprime + (delta_b if heavy_kinetic else 0)
    e4b = 1/(1-4*b)                           # >= e^{4b}
    # 书面证明：10 项（6 个 R、1 个 H'）+ 2 项（5 个 R、1 个 delta H'，仅重块动能缩放时）。
    bracket = 10*r**6*n_der + (2*r**5 if heavy_kinetic else 0)
    current = 4*RANK_PER_DIM*rep_dim*per**3*e4b*bracket
    return dict(b=b, delta_b=delta_b, r=r, per=per, sign=sign_const,
                current=current, n_der=n_der)


def threshold(const, rate):
    """界 const*e^{-rate d} 小于一的距离 ln(const)/rate（Arb 区间）。"""
    return (arb(const.numerator)/const.denominator).log()/(arb(rate.numerator)/rate.denominator)


def check_integer_threshold(table_const, exact_const, rate, n):
    """说明中所写整数门槛：以表中取整后的前因子 C 计，C e^{-rate n} < 1 <= C e^{-rate (n-1)}；
    并核对 C 不小于脚本中的精确前因子，故 d >= n 时精确界也小于一。"""
    assert table_const >= exact_const
    c = arb(table_const.numerator)/table_const.denominator
    r = arb(rate.numerator)/rate.denominator
    assert c*(-r*n).exp() < 1 and c*(-r*(n-1)).exp() >= 1
    return n


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
    assert F(172669, 10**7) < alpha < F(172670, 10**7) and alpha < G
    beta = M_MAX + 8 + SPIN
    heavy_der = M_MAX + 7 + SPIN               # ||gamma5(1-m-K)||，m in (0,M_MAX]
    c = stage_constants(alpha, beta, max(F(1), heavy_der), True, REP_DIM_SUM)
    assert F(156027, 10**4) < beta < F(156029, 10**4)
    assert c["b"] == F(1523793, 10**9) and 656 < 1/c["b"] < 657
    assert c["sign"] < 1155
    assert F(4983, 1000) < c["per"] < F(4984, 1000) and c["r"] < F(11583, 100) and c["n_der"] < F(146115, 10**4)
    assert c["current"] < F(1002, 1000)*F(10)**20
    # 本方法可取的 b 的上确界 ln(1+alpha/(4 sqrt2))（不是最优衰减率）：
    # 用 x-x^2/2 < ln(1+x) < x-x^2/2+x^3/3 与 sqrt2 的上下界。
    x_lo, x_hi = alpha/(4*SQRT2_HI), alpha/(4*F(1414213, 10**6))
    assert F(1414213, 10**6)**2 < 2
    sup_lo, sup_hi = x_lo-x_lo**2/2, x_hi-x_hi**2/2+x_hi**3/3
    assert F(30475, 10**7) < sup_lo and sup_hi < F(30478, 10**7)
    d_sign, d_cur = threshold(c["sign"], c["b"]), threshold(c["current"], 2*c["b"])
    assert 4620 < d_sign < 4640 and 15100 < d_cur < 15120
    check_integer_threshold(F(1155), c["sign"], c["b"], 4628)
    check_integer_threshold(F(1002, 1000)*F(10)**20, c["current"], 2*c["b"], 15112)
    print("[普适常数壳] alpha, beta", float(alpha), float(beta), " b 上确界区间", float(sup_lo), float(sup_hi))
    print("  b、1/b", float(c["b"]), float(1/c["b"]), " 符号核前因子", float(c["sign"]),
          f" 输运项前因子 {float(c['current']):.4e}", " 界小于一的门槛", d_sign.str(6), d_cur.str(6))
    # 两钟绕行情形：外部证书 ch39-shell-winding-clocks 给核本身奇异值 > 0.008
    # （前提：周期格每边 >= 4097、相邻格点质量差 <= 2^-22、各格 delta in [0,1e-8]），
    # Gershgorin 范数界对任意相位成立，故符号核同样有显式界。
    w = stage_constants(F(8, 1000), beta, F(1), False, REP_DIM_SUM)
    assert w["b"] == F(706571, 10**9) and w["sign"] < 2488
    d_w = threshold(w["sign"], w["b"])
    assert 11060 < d_w < 11070
    check_integer_threshold(F(2488), w["sign"], w["b"], 11067)
    print("[两钟绕行，隙 0.008] b", float(w["b"]), " 符号核前因子", float(w["sign"]), " 门槛", d_w.str(6))
    return alpha, beta, c


def verify_two_zeros(shell):
    """11431 与 133931 的实际壳：常数壳与逐格切换。"""
    ctx.prec = 192
    geo = shell["geo"]
    center, scale = arb(CENTER.numerator)/CENTER.denominator, arb(SCALE.numerator)/SCALE.denominator
    heights = {n: acb.zeta_zero(n).imag for n in (11431, 133931)}
    free_gap, masses_hi, diff_norm, spec_norm = {}, {}, {}, {}
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
            # 轻块 m in (-2,0)：min_k|m+2k| = 1-|m+1|；重块 m>0：min_k|m+2k| = m。
            fg = min(min(abs(m+2*k) for k in range(5)) for m in masses)
            for m in masses:
                if m < 0:
                    assert abs(min(abs(m+2*k) for k in range(5))-(1-abs(m+1))) < arb("1e-40")
                else:
                    assert abs(min(abs(m+2*k) for k in range(5))-m) < arb("1e-40")
            free_gap[q, n] = fg
            masses_hi[q, n] = max(masses)
            assert masses[-1] < arb(M_MAX.numerator)/M_MAX.denominator
        # 换点质量差的谱范数结构上界：两类时钟边各六条，从六个互异外站连到三个内站、
        # 每内站两条，故带相位邻接块 A 满足 A^*A = 2I_3，||A + A^*|| = sqrt2；
        # 三条辐边权差至多 1e-8，每条单边拉普拉斯范数为 2，合计至多 6e-8。
        diff_norm[q] = (arb(2).sqrt()*sum(abs(clocks[11431][dd]**q-clocks[133931][dd]**q)
                                           for dd in (724, 804)) + 6*arb("1e-8"))/scale
        # 与 Arb 实算谱范数比较（厄米差矩阵的本征值，含重根簇的严格包络）。
        ev = (laps[11431]-laps[133931]).eig(multiple=True)
        spec_norm[q] = max(abs(e) for e in ev)/scale
        assert spec_norm[q] < diff_norm[q]
    table = {1: ("0.2969081983", "0.2961239655", "0.002875828178", "0.002336310825"),
             2: ("0.1980393569", "0.1990920853", "0.005751615049", "0.005694052288"),
             3: ("0.1353842326", "0.1354330226", "0.008627344306", "0.007301916158"),
             4: ("0.1104950524", "0.1123251024", "0.01150298714", "0.01086964161"),
             6: ("0.1049393951", "0.1032375969", "0.01725389833", "0.01585502834")}
    for q, (g1, g2, vq, sq) in table.items():
        for val, ref in ((free_gap[q, 11431], g1), (free_gap[q, 133931], g2),
                         (diff_norm[q], vq), (spec_norm[q], sq)):
            assert abs(val-arb(ref)) < arb("1e-10"), (q, val, ref)
        print("  荷", q, "自由隙 11431/133931", free_gap[q, 11431].str(10), free_gap[q, 133931].str(10),
              " 换点差结构上界/Arb 谱范数", diff_norm[q].str(10), spec_norm[q].str(10))
    loss = arb(1)-(arb(1)-arb(30)/200).sqrt()
    alpha_c = min(free_gap[k]-loss for k in free_gap)      # 两枚零点任一的常数壳
    alpha_ref = min(free_gap[q, 11431]-loss for q in (1, 2, 3, 4, 6))
    alpha_s = min(free_gap[q, 11431]-loss-diff_norm[q] for q in (1, 2, 3, 4, 6))
    vmax = max(diff_norm.values())
    mmax = max(masses_hi.values())
    assert arb(".025192") < alpha_c < arb(".025193") and arb(".026893") < alpha_ref < arb(".026894")
    assert arb(".0096399") < alpha_s < arb(".0096400")
    assert vmax < arb(".017254")
    a_c, a_ref, a_s = F(2519, 10**5), F(2689, 10**5), F(963, 10**5)
    m_hi = arb_hi(mmax, 10**6)
    assert m_hi < F(28044, 10**4)
    beta = max(m_hi, F(0)) + 8 + SPIN
    heavy_der = m_hi + 7 + SPIN
    const = stage_constants(a_c, beta, max(F(1), heavy_der), True, REP_DIM_SUM)
    # 逐格切换：第一段 H_s = H_ref + (1-s) gamma5 V，V 为逐站对角质量差，
    # 全段隙 >= alpha_ref - ||V||；第二段为 11431 常数壳的标准路径。两段相加，取较小率。
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
    assert 450 < 1/const["b"] < 451 and 1175 < 1/stage1["b"] < 1176
    th = {"const": (threshold(const["sign"], const["b"]), threshold(const["current"], 2*const["b"])),
          "switch": (threshold(stage1["sign"], stage1["b"]), threshold(total_sw, 2*b_sw))}
    assert 2910 < th["const"][0] < 2930 and 9670 < th["const"][1] < 9690
    assert 8740 < th["switch"][0] < 8760 and 25340 < th["switch"][1] < 25370
    check_integer_threshold(F(652), const["sign"], const["b"], 2918)
    check_integer_threshold(F(4705, 1000)*F(10)**18, const["current"], 2*const["b"], 9678)
    check_integer_threshold(F(1698), stage1["sign"], stage1["b"], 8746)
    check_integer_threshold(F(5349, 1000)*F(10)**18, total_sw, 2*b_sw, 25355)
    print("  表中整数门槛（取整前因子，d >= n 时界 < 1）：4628 15112 2918 9678 8746 25355 绕行 11067")
    for key, val in (("const", const), ("stage1", stage1), ("stage2", stage2)):
        print("   ", key, "b =", val["b"], " sign <", float(val["sign"]), f" current < {float(val['current']):.6e}")
    print("[11431 或 133931 常数壳] alpha, beta", float(a_c), float(beta), " b", float(const["b"]),
          " 1/b", float(1/const["b"]), " 符号核前因子", float(const["sign"]),
          f" 输运项前因子 {float(const['current']):.4e}", " 门槛", th["const"][0].str(6), th["const"][1].str(6))
    print("[逐格切换] 第一段隙", float(a_s), " ||V||<=", float(v_hi), " b", float(stage1["b"]),
          " 1/b", float(1/stage1["b"]), " 符号核前因子", float(stage1["sign"]),
          f" 两段合计输运项前因子 {float(total_sw):.4e}", " 门槛", th["switch"][0].str(6), th["switch"][1].str(6))


# ------------------------------------------------------------ 数值检验引理
def numeric_lemmas():
    """仅数值：4^4 周期格、单壳分量、U(1) 链。检验两条可能失败的引理：
    (i) 任意幺正链（这里取均匀随机相位）上 ||W H W^{-1} - H|| <= 4 sqrt2 (e^b-1)；
    (ii) 小曲率链、轻块与重块质量下，围道点 z 上
         ||W (z-H)^{-1} W^{-1}|| <= 1/(dist(z,spec H) - delta_b)。"""
    import numpy as np
    rng = np.random.default_rng(20261126)
    L = 4
    sites = [(a, b_, c, d) for a in range(L) for b_ in range(L) for c in range(L) for d in range(L)]
    idx = {x: k for k, x in enumerate(sites)}
    N = len(sites)
    s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
    Z2, E2 = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z2, -1j*s], [1j*s, Z2]]) for s in (s1, s2, s3)] + [np.block([[Z2, E2], [E2, Z2]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    G5 = np.kron(np.eye(N), g5)

    def H_of(m, th):
        K = 4*np.eye(4*N, dtype=complex)
        for mu in range(4):
            T = np.zeros((N, N), complex)
            for x in sites:
                y = list(x); y[mu] = (y[mu]+1) % L
                T[idx[x], idx[tuple(y)]] = np.exp(1j*th[idx[x], mu])
            K -= np.kron(T, (np.eye(4)-gam[mu])/2) + np.kron(T.conj().T, (np.eye(4)+gam[mu])/2)
        return G5@(m*np.eye(4*N)+K)

    def dist(x, y):
        return sum(min((a-b_) % L, (b_-a) % L) for a, b_ in zip(x, y))

    phis = [np.repeat([dist(x, sites[0]) for x in sites], 4),
            np.repeat([min(dist(x, sites[0]), dist(x, sites[85])+1) for x in sites], 4)]
    worst_i = 0.0
    th_rand = rng.uniform(-np.pi, np.pi, (N, 4))
    H = H_of(-0.6, th_rand)
    for phi in phis:
        for bb in (0.01, 0.2, 1.0):
            W, Wi = np.exp(bb*phi)[:, None], np.exp(-bb*phi)[None, :]
            dev = np.linalg.norm(W*H*Wi - H, 2)
            ratio = dev/(4*np.sqrt(2)*(np.exp(bb)-1))
            worst_i = max(worst_i, ratio)
            assert ratio <= 1+1e-12
    worst_ii = 0.0
    th_small = rng.uniform(-1, 1, (N, 4))/1600
    for m in (-0.6, 0.8):
        H = H_of(m, th_small)
        ev = np.linalg.eigvalsh(H)
        alpha, beta = np.min(np.abs(ev)), np.max(np.abs(ev))
        bb = np.log(1+alpha/(8*np.sqrt(2)))
        dlt = 4*np.sqrt(2)*(np.exp(bb)-1)
        W, Wi = np.exp(bb*phis[0])[:, None], np.exp(-bb*phis[0])[None, :]
        for z in (0j, 0.5*alpha+1j*alpha, alpha*1j, beta+alpha, 0.5*(alpha+beta)+1j*alpha):
            R = np.linalg.inv(z*np.eye(4*N)-H)
            dz = np.min(np.abs(z-ev))
            ratio = np.linalg.norm(W*R*Wi, 2)*(dz-dlt)
            worst_ii = max(worst_ii, ratio)
            assert ratio <= 1+1e-10
    assert 0.60 < worst_i < 0.61 and 0.66 < worst_ii < 0.67
    print(f"[仅数值 4^4] 加权扰动比 max {worst_i:.4f} <= 1；Combes-Thomas 预解式比 max {worst_ii:.4f} <= 1")
    return worst_i, worst_ii


if __name__ == "__main__":
    shell = run_path(str(ROOT/"ch39-three-ray-chiral-boundary.py"))
    check_gershgorin(shell)
    verify_contraction_identity()
    verify_universal()
    verify_two_zeros(shell)
    numeric_lemmas()
    print("全部断言通过")

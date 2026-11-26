"""外尔测度流的曲率部分显式指数局域，而对称部分可随参考选择任意放大、衰减任意变慢。

说明：外尔测度流的曲率部分显式指数局域而对称部分可随参考选择任意放大且衰减任意变慢.md

1. 书面证明（说明文）：基联络的曲率恒等式即 Lüscher 测度项的可积性条件（文献结果，
   说明中为约定一致重证）；新的是 Kato 输运项曲率 = F_{P_1} - F_{P_0}（一步 Stokes）
   与曲率的显式常数。curvature_constants 按书面推导的公式（2 项、5 个预解式、
   块秩 4*rank*r、l^3、e^{4b}）计算常数。
2. 精确有理：普适常数壳、两枚零点常数壳、逐格切换、参考端点四种情形的率与前因子，
   界开始小于一的整数门槛（对数用 Arb），以及与上一篇定理 B 精确前因子之比（只比较界）。
   谱隙与范数输入取自 ch39-shell-overlap-kato-locality.py 的证书，本脚本重新运行其断言。
3. 精确整数：四维中含一条链的方格数为 6；两链基点距离 d >= 3 时无公共方格，
   所含方格中心的距离 <= d+2。
4. 仅数值：(i) 有限维玩具族上核对两条曲率恒等式；(ii) 3^4 周期格小随机容许 U(1) 场上，
   十五分量端点曲率之和不为零、按场幅三次方标度，并给出规范变换与纯规范场的噪声；
   (iii) 4^4 格上相位泛函 Phi 的格点对称与二阶变分闭式。
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import itertools
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb

ROOT = Path(__file__).parent
PREV = run_path(str(ROOT/"ch39-shell-overlap-kato-locality.py"))
SQRT2_HI, PI_LO = PREV["SQRT2_HI"], PREV["PI_LO"]
EPS, SPIN, sqrt_lo = PREV["EPS"], PREV["SPIN"], PREV["sqrt_lo"]
REP_DIM_SUM = PREV["REP_DIM_SUM"]
assert REP_DIM_SUM == 15 and EPS == F(1, 200)


# ------------------------------------------------------------ 曲率常数
def curvature_constants(alpha, beta, rank_per_dim, rep_dim):
    """F_P(eta,zeta) = Tr(P[delta_eta P, delta_zeta P])，对两条相距 d 的链：
    |F| <= 2 * 4*rank*r * l^3 * e^{4b} * rho^5 * ||eta|| ||zeta|| e^{-2bd}。
    b 的取法与上一篇相同（delta_b <= alpha/2）。rank_per_dim：一条链变分块的秩除以 r。"""
    t = alpha/(8*SQRT2_HI+alpha)
    b = F(int(t*10**9), 10**9)
    delta_b = 4*SQRT2_HI*b/(1-b)                # >= 4 sqrt2 (e^b - 1)
    assert 0 < b < F(1, 10) and delta_b <= alpha/2
    rho = 1/(alpha-delta_b)
    assert rho <= 2/alpha
    per = (2*(beta+alpha)+4*alpha)/(2*PI_LO)    # >= |Gamma|/(2 pi)
    e4b = 1/(1-4*b)                             # >= e^{4b}
    const = 2*4*rank_per_dim*rep_dim*per**3*e4b*rho**5
    return dict(b=b, rho=rho, per=per, const=const)


# 说明中表列的取整前因子（均不小于脚本中的精确前因子）。
PREV_EXACT = {}

TABLE = {"shell_universal": F(1184, 1000)*F(10)**16, "shell_two_zero": F(1003, 1000)*F(10)**15,
         "shell_switch": F(1208, 1000)*F(10)**17, "reference": F(2606, 1000)*F(10)**6}


def threshold(const, rate):
    return (arb(const.numerator)/const.denominator).log()/(arb(rate.numerator)/rate.denominator)


def check_integer_threshold(table_const, exact_const, rate, n):
    assert table_const >= exact_const
    c = arb(table_const.numerator)/table_const.denominator
    r = arb(rate.numerator)/rate.denominator
    assert c*(-r*n).exp() < 1 and c*(-r*(n-1)).exp() >= 1


def certify_constants():
    shell = run_path(str(ROOT/"ch39-three-ray-chiral-boundary.py"))
    PREV["check_gershgorin"](shell)
    alpha_u, beta_u, kato_u = PREV["verify_universal"]()
    # 证书：alpha 2519e-5, 2689e-5, 963e-5，beta < 12.8119。记录其中三次 stage_constants
    # 调用（常数壳、切换第一段、第二段）的精确前因子，供比较。
    g = PREV["verify_two_zeros"].__globals__
    orig, calls = g["stage_constants"], []

    def rec(*a, **k):
        out = orig(*a, **k); calls.append(out); return out
    g["stage_constants"] = rec
    try:
        PREV["verify_two_zeros"](shell)
    finally:
        g["stage_constants"] = orig
    assert len(calls) == 3
    PREV_EXACT.update(const=calls[0], stage1=calls[1], stage2=calls[2])
    beta_z = F(128119, 10**4)
    # 端点参考：三个标准重叠块，Wilson 质量 -1；H(-1)^2 >= 1-30 eps = 17/20；
    # ||M + sum A|| <= 7（sum A 的谱在 [0,8]），||sum gamma S|| <= sqrt(4.03)。
    alpha_0 = sqrt_lo(1-30*EPS)
    beta_0 = 7 + SPIN
    assert 1-30*EPS == F(17, 20) and F(921954, 10**6) < alpha_0 < F(921955, 10**6)
    assert beta_0 < F(90075, 10**4)
    cases = {
        "shell_universal": curvature_constants(alpha_u, beta_u, 2*19, REP_DIM_SUM),
        "shell_two_zero": curvature_constants(F(2519, 10**5), beta_z, 2*19, REP_DIM_SUM),
        "shell_switch": curvature_constants(F(963, 10**5), beta_z, 2*19, REP_DIM_SUM),
        "reference": curvature_constants(alpha_0, beta_0, 2*3, REP_DIM_SUM),
    }
    su, sz, ss, rf = (cases[k] for k in ("shell_universal", "shell_two_zero", "shell_switch", "reference"))
    # 壳的三种情形与上一篇的符号核、输运项取同一 b；参考端点的率快约五十倍。
    assert su["b"] == F(1523793, 10**9) and sz["b"] == F(2221443, 10**9) and ss["b"] == F(850412, 10**9)
    assert rf["b"] == F(75346220, 10**9) and 49 < rf["b"]/su["b"] < 50
    assert F(11582, 100) < su["rho"] < F(11583, 100) and F(4983, 1000) < su["per"] < F(4984, 1000)
    assert F(2169, 1000) < rf["rho"] < F(2170, 1000) and F(3747, 1000) < rf["per"] < F(3748, 1000)
    for key, bound in TABLE.items():
        assert cases[key]["const"] < bound, key
    assert cases["shell_universal"]["const"] > F(1183, 1000)*F(10)**16
    # 输运项的曲率 = F_{P_1} - F_{P_0}：两前因子相加，率取较小者（壳的 2b）。
    kato_curv = {k: cases[k]["const"]+rf["const"] for k in ("shell_universal", "shell_two_zero", "shell_switch")}
    for key, val in kato_curv.items():
        assert val < TABLE[key], key
    # 与上一篇定理 B 的精确前因子（同一率 2b）之比。定理 B 界的是单个 delta_zeta j_e(eta)，
    # 由它推出的曲率界为 2 C_B，故相对曲率界的缩小倍数是下列比值的两倍。只比较界。
    ratios = {"shell_universal": kato_u["current"]/su["const"]}
    assert 8460 < ratios["shell_universal"] < 8470
    for key, val in (("shell_two_zero", PREV_EXACT["const"]["current"]),
                     ("shell_switch", PREV_EXACT["stage1"]["current"]+PREV_EXACT["stage2"]["current"])):
        ratios[key] = val/cases[key]["const"]
    assert PREV_EXACT["const"]["b"] == sz["b"] and PREV_EXACT["stage1"]["b"] == ss["b"]
    assert 4690 < ratios["shell_two_zero"] < 4700 and 44 < ratios["shell_switch"] < 45
    th = {k: threshold(cases[k]["const"], 2*cases[k]["b"]) for k in cases}
    assert 12140 < th["shell_universal"] < 12150 and 7770 < th["shell_two_zero"] < 7780
    assert 23120 < th["shell_switch"] < 23130 and 98 < th["reference"] < 99
    return cases, kato_curv, th, kato_u, ratios


def integer_thresholds(cases):
    """表中整数 n：取整前因子 C 下 C e^{-2b n} < 1 <= C e^{-2b (n-1)}。"""
    out = {}
    for key, c in TABLE.items():
        rate = 2*cases[key]["b"]
        n = int((threshold(c, rate)).upper().ceil().unique_fmpz())
        if not (arb(c.numerator)/c.denominator*(-arb(rate.numerator)/rate.denominator*n).exp() < 1):
            n += 1
        check_integer_threshold(c, cases[key]["const"], rate, n)
        out[key] = n
    assert out == {"shell_universal": 12145, "shell_two_zero": 7775, "shell_switch": 23126, "reference": 99}, out
    return out


# ------------------------------------------------------------ 精确组合
def plaquette_geometry(L=6):
    """四维周期格：含一条链的方格恰 6 个；两链基点距离 d >= 3 时无公共方格，
    且所含方格中心（半整数坐标）的周期 l^1 距离 <= d+2（对 L=6 全部链对逐一核对）。"""
    sites = list(itertools.product(range(L), repeat=4))

    def add(x, mu, s=1):
        y = list(x); y[mu] = (y[mu]+s) % L; return tuple(y)

    def dist(x, y):
        return sum(min((a-c) % L, (c-a) % L) for a, c in zip(x, y))

    def plaqs_of(x, mu):
        out = []
        for nu in range(4):
            if nu != mu:
                out.append((x, min(mu, nu), max(mu, nu)))          # 基点 x
                out.append((add(x, nu, -1), min(mu, nu), max(mu, nu)))  # 基点 x - nu
        return out

    def links_of(p):
        x, m, n = p
        return {(x, m), (add(x, m), n), (add(x, n), m), (x, n)}

    for mu in range(4):
        ps = plaqs_of(sites[0], mu)
        assert len(set(ps)) == 6 and all((sites[0], mu) in links_of(p) for p in ps)
    def center2(p):                      # 方格中心坐标的两倍（模 2L）
        x, m, n = p
        return tuple(2*x[k]+(k == m)+(k == n) for k in range(4))

    def cdist(p, q):
        a, c = center2(p), center2(q)
        return F(sum(min((u-v) % (2*L), (v-u) % (2*L)) for u, v in zip(a, c)), 2)

    x0 = sites[0]
    for mu in range(4):                  # 中心到所含链基点的距离恰为 1
        for p in plaqs_of(x0, mu):
            c = center2(p)
            assert sum(min(u % (2*L), (-u) % (2*L)) for u in c) == 2
    worst = 0
    for y in sites:
        d = dist(x0, y)
        if d < 3:
            continue
        for mu in range(4):
            for nu in range(4):
                pe, pf = plaqs_of(x0, mu), plaqs_of(y, nu)
                assert not (set(pe) & set(pf))
                assert all((y, nu) not in links_of(p) for p in pe)
                worst = max(worst, max(cdist(p, q) - d for p in pe for q in pf))
    assert worst == 2
    print("[精确] 含一链的方格数 6；d>=3 时无公共方格，方格中心距离 - d 的最大值", worst)
    return worst


# ------------------------------------------------------------ 仅数值
def numeric_identities():
    """玩具族 H(s,t1,t2) = H0 + s H1 + t1 K1 + t2 K2（8 维，秩 4）：
    (a) j_k = int_0^1 Tr(P[d_s P, d_k P]) ds 满足 d_1 j_2 - d_2 j_1 = F_{P_1} - F_{P_0}；
    (b) 光滑正交基 V 的联络 L_k = Tr(V^* d_k V) 满足 d_1 L_2 - d_2 L_1 = F_P(1,2)。"""
    import numpy as np
    rng = np.random.default_rng(1)
    n = 8

    def herm(s):
        A = rng.standard_normal((n, n))+1j*rng.standard_normal((n, n)); return s*(A+A.conj().T)/2
    H0 = np.diag([3, 2.5, 2, 1.5, -1.5, -2, -2.5, -3]).astype(complex)+herm(0.2)
    H1, K1, K2 = herm(0.4), herm(0.3), herm(0.3)

    def Hf(s, t1, t2):
        return H0+s*H1+t1*K1+t2*K2

    def pd(H, Ks):
        lam, W = np.linalg.eigh(H); pos = lam > 0
        assert pos.sum() == 4 and np.min(np.abs(lam)) > 0.5
        Wp, Wm = W[:, pos], W[:, ~pos]; den = lam[pos][:, None]-lam[~pos][None, :]
        P = Wp@Wp.conj().T
        out = []
        for K in Ks:
            X = Wp@((Wp.conj().T@K@Wm)/den)@Wm.conj().T
            out.append(X+X.conj().T)
        return P, out

    def curv(P, a, b):
        return np.trace(P@(a@b-b@a))

    xs, ws = np.polynomial.legendre.leggauss(60); xs, ws = (xs+1)/2, ws/2

    def j(k, t1, t2):
        tot = 0
        for x, w in zip(xs, ws):
            P, (dS, dk) = pd(Hf(x, t1, t2), [H1, (K1, K2)[k]])
            tot += w*np.trace(P@(dS@dk-dk@dS))
        return tot
    h = 1e-5
    lhs = (j(1, h, 0)-j(1, -h, 0))/(2*h)-(j(0, 0, h)-j(0, 0, -h))/(2*h)
    P0, (a0, b0) = pd(Hf(0, 0, 0), [K1, K2]); P1, (a1, b1) = pd(Hf(1, 0, 0), [K1, K2])
    F0, F1 = curv(P0, a0, b0), curv(P1, a1, b1)
    assert abs(lhs-(F1-F0)) < 1e-8*abs(F1-F0) + 1e-12 and abs(F1-F0) > 1e-2
    V0 = np.linalg.eigh(Hf(1, 0, 0))[1][:, 4:]

    def basis(t1, t2):
        P, _ = pd(Hf(1, t1, t2), []); X = P@V0
        w, U = np.linalg.eigh(X.conj().T@X); return X@U@np.diag(w**-0.5)@U.conj().T

    def Lc(k, t1, t2):
        e = ((h, 0), (0, h))[k]
        return np.trace(basis(t1, t2).conj().T@(basis(t1+e[0], t2+e[1])-basis(t1-e[0], t2-e[1]))/(2*h))
    H2 = 1e-4
    lhsL = (Lc(1, H2, 0)-Lc(1, -H2, 0))/(2*H2)-(Lc(0, 0, H2)-Lc(0, 0, -H2))/(2*H2)
    assert abs(lhsL-F1) < 1e-3*abs(F1) and abs(F1) > 1e-3
    print(f"[仅数值 玩具族] 输运项曲率 {lhs.imag:.10f}i vs F1-F0 {(F1-F0).imag:.10f}i；"
          f"基联络曲率 {lhsL.imag:.7f}i vs F1 {F1.imag:.7f}i")


def numeric_endpoint_curvature():
    """3^4 周期格、Wilson 质量 -1 的标准重叠核、U(1) 链。十五分量（右手单态按模型约定
    取 1-P，等价地改写为荷 -y 的左手共轭）端点曲率之和 sum_c sigma_c y_c^2 g(U^{y_c})。"""
    import numpy as np
    L = 3
    V = L**4
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    Z, E = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z, -1j*sk], [1j*sk, Z]]) for sk in s]+[np.block([[Z, E], [E, Z]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    sites = list(itertools.product(range(L), repeat=4))
    idx = {x: k for k, x in enumerate(sites)}
    nb = [[idx[tuple((x[m]+(mu == m)) % L for m in range(4))] for mu in range(4)] for x in sites]
    Pm = [(np.eye(4)-g)/2 for g in gam]; Pp = [(np.eye(4)+g)/2 for g in gam]
    rng = np.random.default_rng(20261126)
    base = rng.standard_normal((V, 4))

    def H_of(theta, y):
        D = 3*np.eye(4*V, dtype=complex)            # m + 4 = 3
        for i in range(V):
            for mu in range(4):
                j = nb[i][mu]; u = np.exp(1j*y*theta[i, mu])
                D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*u
                D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*np.conj(u)
        return np.kron(np.eye(V), g5)@D

    def K_of(theta, y, e):
        i, mu = e; j = nb[i][mu]; u = np.exp(1j*y*theta[i, mu])
        D = np.zeros((4*V, 4*V), complex)
        D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*1j*u
        D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*(-1j)*np.conj(u)
        return np.kron(np.eye(V), g5)@D

    def g(theta, y, e, f):
        lam, W = np.linalg.eigh(H_of(theta, y)); pos = lam > 0
        assert pos.sum() == 2*V and np.min(np.abs(lam)) > np.sqrt(0.85)
        Wp, Wm = W[:, pos], W[:, ~pos]; den = lam[pos][:, None]-lam[~pos][None, :]
        A = (Wp.conj().T@K_of(theta, y, e)@Wm)/den; B = (Wp.conj().T@K_of(theta, y, f)@Wm)/den
        return np.trace(A@B.conj().T)-np.trace(B@A.conj().T)

    def max_plaq(theta):
        mx = 0
        for i in range(V):
            for mu in range(4):
                for nu in range(mu+1, 4):
                    ang = theta[i, mu]+theta[nb[i][mu], nu]-theta[nb[i][nu], mu]-theta[i, nu]
                    mx = max(mx, abs(ang))
        return mx
    e, f = (0, 0), (idx[(0, 1, 0, 0)], 2)
    left = {1: 6, -3: 2}; right = {4: 3, -2: 3, -6: 1}
    res = {}
    for amp in (5e-5, 1e-4):
        theta = amp*base
        # 容许：|q|=6 表示中 ||1-U_p|| = 2|sin(3 theta_p)| <= 1/200。
        assert 2*np.sin(3*max_plaq(theta)) <= 1/200
        gs = {y: g(theta, y, e, f) for y in (1, -3, 4, -2, -6, 3, -4, 2, 6)}
        model = sum(m*y*y*gs[y] for y, m in left.items())-sum(m*y*y*gs[y] for y, m in right.items())
        conj = sum(m*y*y*gs[y] for y, m in left.items())+sum(m*y*y*gs[-y] for y, m in right.items())
        res[amp] = (model, conj, gs[1])
    (m1, c1, g1), (m2, c2, g2) = res[5e-5], res[1e-4]
    assert abs(m2.real) < 1e-15 and abs(c2.real) < 1e-15
    # 合计约 1.6e-12，按幅的三次方标度；零场浮点噪声约 1e-16。
    assert 1.5e-12 < m2.imag < 1.7e-12 and 7.9 < m2.imag/m1.imag < 8.1
    assert abs(c2-m2) < 1e-14 and abs(c2-m2) < 1e-3*abs(m2)
    assert 4.5e-9 < -g2.imag < 4.7e-9      # 荷一单分量（未乘 y^2 与重数）
    theta0 = 0*base
    z0 = sum(m*y*y*g(theta0, y, e, f) for y, m in left.items())-sum(m*y*y*g(theta0, y, e, f) for y, m in right.items())
    assert abs(z0) < 1e-14
    # 噪声估计：O(1) 随机规范变换 theta_mu(x) + lam(x) - lam(x+mu) 前后之差，以及纯规范场上的值。
    lam = np.random.default_rng(7).uniform(-np.pi, np.pi, V)
    grad = np.array([[lam[i]-lam[nb[i][mu]] for mu in range(4)] for i in range(V)])

    def model_total(theta):
        return (sum(m*y*y*g(theta, y, e, f) for y, m in left.items())
                - sum(m*y*y*g(theta, y, e, f) for y, m in right.items()))
    gauge_diff = abs(model_total(1e-4*base+grad)-m2)
    pure = abs(model_total(grad))
    assert gauge_diff < 1e-14 and pure < 1e-14 and abs(m2) > 100*max(gauge_diff, pure, abs(z0))
    print(f"[仅数值 3^4] 十五分量端点曲率合计 {m2.imag:.4e}i（幅 1e-4），幅减半比 {m2.imag/m1.imag:.3f}；"
          f"单分量 g(1) = {g2.imag:.3e}i；规范变换前后差 {gauge_diff:.1e}，纯规范 {pure:.1e}，零场 {abs(z0):.1e}")


def numeric_phase_functional():
    """4^4 格：Phi = kappa sum_{p,p'} e^{-b' D(p,p')} w(p) w(p')，w = 1 - cos theta_p，
    D 为方格中心的周期 l^1 距离（Phi 在平移与超立方转动、反射下不变）。
    背景只在链 e、f 上取角 a、a'，二阶变分 = 2 kappa sin a sin a' sum e^{-b'D}
    >= 72 kappa sin a sin a' e^{-b'(d+2)}，以有限差分核对。"""
    import numpy as np
    L = 4
    sites = list(itertools.product(range(L), repeat=4))
    idx = {x: k for k, x in enumerate(sites)}
    nb = [[idx[tuple((x[m]+(mu == m)) % L for m in range(4))] for mu in range(4)] for x in sites]
    plaqs = [(i, m, n) for i in range(len(sites)) for m in range(4) for n in range(m+1, 4)]
    # D：方格中心（半整数坐标）的周期 l^1 距离；用两倍坐标模 2L 计算。
    coords = np.array([[2*sites[i][k]+(k == m)+(k == n) for k in range(4)] for i, m, n in plaqs])
    diff = np.abs(coords[:, None, :]-coords[None, :, :]); diff = np.minimum(diff, 2*L-diff)
    D = diff.sum(-1)/2
    kappa, bprime = 1.0, 0.3
    Wt = np.exp(-bprime*D)

    def phi(theta):
        th = np.array([theta[i, m]+theta[nb[i][m], n]-theta[nb[i][n], m]-theta[i, n] for i, m, n in plaqs])
        w = 1-np.cos(th)
        return kappa*w@Wt@w
    e, f = (0, 0), (idx[(2, 2, 0, 0)], 1)
    d = 4
    a, ap = 1/2000, 1/2000
    assert 2*np.sin(3*a) <= 1/200                 # |q|=6 表示的方格误差
    base = np.zeros((len(sites), 4)); base[e] = a; base[f] = ap
    h = 1e-4

    def at(de, df):
        t = base.copy(); t[e] += de; t[f] += df; return phi(t)
    hess = (at(h, h)-at(h, -h)-at(-h, h)+at(-h, -h))/(4*h*h)
    pe = [k for k, (i, m, n) in enumerate(plaqs) if e in ((i, m), (nb[i][m], n), (nb[i][n], m), (i, n))]
    pf = [k for k, (i, m, n) in enumerate(plaqs) if f in ((i, m), (nb[i][m], n), (nb[i][n], m), (i, n))]
    assert len(pe) == 6 and len(pf) == 6 and not set(pe) & set(pf)
    closed = 2*kappa*np.sin(a)*np.sin(ap)*sum(Wt[p, q] for p in pe for q in pf)
    lower = 72*kappa*np.sin(a)*np.sin(ap)*np.exp(-bprime*(d+2))
    assert abs(hess-closed) < 1e-4*closed and closed >= lower > 0
    # Phi 的格点对称：随机场在平移、0<->1 轴交换与 x0 反射下取值不变。
    rnd = np.random.default_rng(3).uniform(-1, 1, (len(sites), 4))

    def transform(theta, kind):
        out = np.zeros_like(theta)
        for i, x in enumerate(sites):
            for mu in range(4):
                if kind == "shift":
                    out[idx[tuple((x[k]+(k == 2)) % L for k in range(4))], mu] = theta[i, mu]
                elif kind == "swap":
                    y = (x[1], x[0], x[2], x[3]); nu = {0: 1, 1: 0}.get(mu, mu)
                    out[idx[y], nu] = theta[i, mu]
                else:                                   # x0 -> -x0
                    y = ((-x[0]) % L,)+x[1:]
                    if mu == 0:
                        out[idx[((-x[0]-1) % L,)+x[1:]], 0] = -theta[i, 0]
                    else:
                        out[idx[y], mu] = theta[i, mu]
        return out
    p0 = phi(rnd)
    for kind in ("shift", "swap", "reflect"):
        assert abs(phi(transform(rnd, kind))-p0) < 1e-9*abs(p0)
    print(f"[仅数值 4^4] Phi 的二阶变分 {hess:.6e}，闭式 {closed:.6e}，下界 {lower:.6e}")


if __name__ == "__main__":
    cases, kato_curv, th, kato_u, ratios = certify_constants()
    ints = integer_thresholds(cases)
    for k, v in cases.items():
        print(f"[{k}] b = {v['b']}  2b = {float(2*v['b']):.9f}  rho <= {float(v['rho']):.4f}"
              f"  l <= {float(v['per']):.4f}  曲率前因子 < {float(v['const']):.4e}  门槛 {th[k].str(6)}  整数门槛 {ints[k]}")
    for k, v in kato_curv.items():
        print(f"[输运项曲率 {k}] 前因子 < {float(v):.4e}（率取壳的 2b）")
    print("上一篇定理 B 精确前因子 / 本文壳曲率前因子（同率，只比较界）",
          {k: round(float(v), 1) for k, v in ratios.items()})
    plaquette_geometry()
    numeric_identities()
    numeric_endpoint_curvature()
    numeric_phase_functional()
    print("全部断言通过")

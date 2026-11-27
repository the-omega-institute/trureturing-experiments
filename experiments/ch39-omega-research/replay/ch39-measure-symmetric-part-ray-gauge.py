"""外尔测度流对称部分：与相位无关的部分恰为曲率导数，射线平行输运定相下有显式指数界。

说明：凸片上外尔测度流对称部分模相位Hessian由曲率导数决定且纯超荷凸片上的射线定相参考部分以率零点一一三指数局域.md

1. 书面证明（说明文）：
   (a) 旋度恒等式 d_g F(e,f) = d_e S(f,g) - d_f S(e,g)（三条两两可交换的变分）；
   (b) 凸片上的完全分类：对称张量 S 是曲率为 F 的某个连接的对称部分 <=> 旋度 = dF；
       两个这样的连接之差为 -i dPhi，对称部分之差为 -2i Hess Phi；
   (c) 曲率导数的三链界：10 项六预解式 + 至多 2 项五预解式（同链二阶导），
       |d_g F_P(e,f)| <= 88 k r l^3 rho^6 e^{6b} ||.||^3 e^{-b P}，P 为三基点周长；
   (d) 射线平行输运定相 Phi(A) = -i int_0^1 L^ref_{tA}(A) dt：L^rad(A)(A) = 0，
       L^rad(A)(eta) = int_0^1 t F_{tA}(A,eta) dt，对称部分
       S^rad = int_0^1 t^2 [d_zeta F(A,eta) + d_eta F(A,zeta)](tA) dt。
   (e) 规范判据：F(grad w, eta) = i sum w delta_eta A，射线定相规范不变 <=> A = 0 于 K_a；
       射线上有隙时由解析性，K_a 外一点 A != 0 即可判定。
   本脚本只按这些公式计算常数，不复算推导本身。
2. 精确有理：K_a 容许（a = 1/4800）；三链前因子、射线定相下参考对称部分与完整测度
   对称部分的前因子与率；格点求和 sum_{Z^4} e^{-c|n|_1} <= (2(1+c)/c)^4。
3. Arb：界开始小于一的整数门槛。
4. 仅数值：(i) 有限维玩具族上旋度恒等式、射线定相公式与 L^rad(A)(A) = 0；
   (ii) 3^4 周期格小随机容许 U(1) 场上十五分量端点曲率的三链导数不为零；
   (iii) 3^4 格上十五分量手征密度 A 按 eps^4 非零（射线定相规范判据，场不在 K_a 内）。
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
CURV = run_path(str(ROOT/"ch39-shell-measure-curvature-locality.py"))
EPS, REP_DIM_SUM = CURV["EPS"], CURV["REP_DIM_SUM"]
assert EPS == F(1, 200) and REP_DIM_SUM == 15

QMAX = 6                      # 超荷 |q| 的最大值；坐标方向 d/dA_g 在各表示中的像范数之最大值
A_RADIUS = F(1, 4800)         # K_a = {A : max_g |A_g| <= a}
LAMBDA = F(1, 2)              # 格点求和中把周长拆为 lambda*d + (1-lambda)*d(z,x)
TERMS_SIX, TERMS_FIVE = 10, 2  # 曲率导数：十项六预解式，同链时至多两项五预解式


def admissible_patch():
    """K_a 中每点（及射线上每点）在五种超荷表示中满足 ||1-U_p|| <= eps：
    |theta_p| <= 4 |q| a，||1-e^{i theta}|| = 2|sin(theta/2)| <= |theta|。"""
    assert 4*QMAX*A_RADIUS == EPS
    return 4*QMAX*A_RADIUS


def three_link_constant(case, rank_per_dim):
    """|sum_c sigma_c d_theta F_{P^c}(eta,zeta)| <= C3 ||theta|| ||eta|| ||zeta|| e^{-b P}。
    六预解式项：8 组块，每组 <= k r ||.||^3 rho^6 e^{-b(P-6)}；
    五预解式项：4 组块，每组 <= k r ||.||^3 rho^5 e^{-b(P-4)} <= 同上（rho >= 1）。"""
    b, rho, per = case["b"], case["rho"], case["per"]
    assert rho >= 1 and 6*b < 1
    e6b = 1/(1-6*b)                               # >= e^{6b}
    const = (TERMS_SIX*8+TERMS_FIVE*4)*rank_per_dim*REP_DIM_SUM*per**3*rho**6*e6b
    assert TERMS_SIX*8+TERMS_FIVE*4 == 88
    return const


def lattice_sum_hi(c):
    """sum_{n in Z^4} e^{-c |n|_1} = coth(c/2)^4 <= (2(1+c)/c)^4（e^{-c} <= 1/(1+c)）。
    周期格：每点取各坐标在 (-L/2, L/2] 的代表，单射且保持周期 l^1 距离，故同一上界成立。"""
    return (2*(1+c)/c)**4


def ray_symmetric_constant(case, c3):
    """|S^rad(eta,zeta)| <= (1/3) * a * 2 * QMAX * C3 * e^{-b d} * sum_g e^{-b(d(z,x)+d(z,y))}
    <= (2/3) a QMAX C3 * 4 * lattice_sum((1-lambda) b) * e^{-(1+lambda) b d}。"""
    b = case["b"]
    const = F(2, 3)*A_RADIUS*QMAX*c3*4*lattice_sum_hi((1-LAMBDA)*b)
    return const, (1+LAMBDA)*b


def first_threshold(const, rate):
    c = arb(const.numerator)/const.denominator
    r = arb(rate.numerator)/rate.denominator
    n = int(((c.log())/r).upper().ceil().unique_fmpz())
    while not (c*(-r*n).exp() < 1):
        n += 1
    while n > 1 and c*(-r*(n-1)).exp() < 1:
        n -= 1
    assert c*(-r*n).exp() < 1 and c*(-r*(n-1)).exp() >= 1
    return n


def certify():
    cases, kato_curv, th, kato_u, ratios = CURV["certify_constants"]()
    prev_exact = CURV["PREV_EXACT"]
    rf, su = cases["reference"], cases["shell_universal"]
    assert rf["b"] == F(75346220, 10**9) and su["b"] == F(1523793, 10**9)
    admissible_patch()
    out = {}
    # 三链界：参考端点（三个标准块，k = 2*3）与完整测度（十九点壳，k = 2*19）。
    c3_ref = three_link_constant(rf, 2*3)
    c3_sh = three_link_constant(su, 2*19)
    assert F(7928, 1000)*10**7 < c3_ref < F(7929, 1000)*10**7
    assert F(1512, 1000)*10**19 < c3_sh < F(1513, 1000)*10**19
    out["c3_ref"], out["c3_sh"] = c3_ref, c3_sh
    # 周长 P >= 2 d(x,y)，故用 e^{-b P} <= e^{-2 b d}：门槛按 P 报。
    out["n3_ref"] = first_threshold(c3_ref, rf["b"])
    out["n3_sh"] = first_threshold(c3_sh, su["b"])
    assert out["n3_ref"] == 242 and out["n3_sh"] == 28983, (out["n3_ref"], out["n3_sh"])
    # 射线定相下参考连接的对称部分。
    cs_ref, rate_ref = ray_symmetric_constant(rf, c3_ref)
    assert rate_ref == F(3, 2)*rf["b"] == F(11301933, 10**8)
    assert F(2433, 1000)*10**12 < cs_ref < F(2434, 1000)*10**12
    out["cs_ref"], out["rate_ref"] = cs_ref, rate_ref
    out["ns_ref"] = first_threshold(cs_ref, rate_ref)
    # 有理上界非最优：精确 e^{6b} 与精确 coth 下的同一表达式（仍用 rho、l 的有理上界）。
    ab = arb(rf["b"].numerator)/rf["b"].denominator
    tight3 = (arb(c3_ref.numerator)/c3_ref.denominator)*(1-6*ab)*(6*ab).exp()
    assert tight3 > 682*10**5 and tight3 < 684*10**5
    cth = (1/(ab/4).tanh())**4
    tight_s = (arb(cs_ref.numerator)/cs_ref.denominator)/(arb(lattice_sum_hi(rf["b"]/2).numerator)
                                                           / lattice_sum_hi(rf["b"]/2).denominator)*cth
    assert tight_s > 209*10**10 and tight_s < 211*10**10
    out["tight3"], out["tight_s"] = tight3, tight_s
    assert out["ns_ref"] == 253, out["ns_ref"]
    # 完整测度流对称部分 = 参考部分（射线定相）+ Kato 输运项的对称部分；
    # 后者 = delta_zeta j_e(eta) + delta_eta j_f(zeta)，各由输运项定理 B 界定（率 2b，15 分量）。
    totals = {}
    for key, cur, b in (("shell_universal", kato_u["current"], su["b"]),
                        ("shell_two_zero", prev_exact["const"]["current"], prev_exact["const"]["b"]),
                        ("shell_switch", prev_exact["stage1"]["current"]+prev_exact["stage2"]["current"],
                         prev_exact["stage1"]["b"])):
        rate = 2*b
        assert rate < rate_ref                        # 参考部分衰减更快，可并入较慢的率
        const = cs_ref+2*cur
        totals[key] = dict(const=const, rate=rate, n=first_threshold(const, rate), kato=2*cur)
    tu, tz, ts = totals["shell_universal"], totals["shell_two_zero"], totals["shell_switch"]
    assert tu["const"] < F(2004, 1000)*10**20 and tz["const"] < F(9410, 1000)*10**18 \
        and ts["const"] < F(1070, 100)*10**18
    assert tu["const"] > F(2003, 1000)*10**20
    assert (tu["n"], tz["n"], ts["n"]) == (15339, 9834, 25762), (tu["n"], tz["n"], ts["n"])
    # 参考部分在总前因子中所占比例极小。
    assert cs_ref/tu["const"] < F(1, 10**7)
    out["totals"] = totals
    out["cases"] = cases
    return out


# ------------------------------------------------------------ 仅数值
def numeric_toy():
    """8 维、秩 4 的投影族 H(A) = H0 + sum_i A_i K_i（A in R^3）。参考基：固定标架投影后
    正交化；连接 L_i = Tr(V^* d_i V)。核对 (a) 旋度恒等式；(b) 射线定相 L^rad = L - i dPhi
    满足 L^rad(A)(A) = 0 与 S^rad = int t^2 [d_j F(A,e_i) + d_i F(A,e_j)](tA) dt。"""
    import numpy as np
    rng = np.random.default_rng(11)
    n = 8

    def herm(s):
        M = rng.standard_normal((n, n))+1j*rng.standard_normal((n, n)); return s*(M+M.conj().T)/2
    H0 = np.diag([3, 2.5, 2, 1.5, -1.5, -2, -2.5, -3]).astype(complex)+herm(0.2)
    Ks = [herm(0.3) for _ in range(3)]
    V0 = np.linalg.eigh(H0)[1][:, 4:]+0.3*(rng.standard_normal((n, 4))+1j*rng.standard_normal((n, 4)))

    def proj(A):
        H = H0+sum(a*K for a, K in zip(A, Ks))
        lam, W = np.linalg.eigh(H); assert np.min(np.abs(lam)) > 0.5 and (lam > 0).sum() == 4
        Wp = W[:, lam > 0]; return Wp@Wp.conj().T

    def basis(A):
        X = proj(A)@V0
        w, U = np.linalg.eigh(X.conj().T@X); return X@U@np.diag(w**-0.5)@U.conj().T

    h = 1e-4
    E = np.eye(3)

    def Lc(A, i):
        A = np.asarray(A, float)
        return np.trace(basis(A).conj().T@(basis(A+h*E[i])-basis(A-h*E[i]))/(2*h))

    def Fc(A, i, j):
        P = proj(A)
        dPi = (proj(A+h*E[i])-proj(A-h*E[i]))/(2*h); dPj = (proj(A+h*E[j])-proj(A-h*E[j]))/(2*h)
        return np.trace(P@(dPi@dPj-dPj@dPi))
    A0 = np.array([0.31, -0.22, 0.17])
    H2 = 1e-3

    def dd(fun, A, k):
        return (fun(A+H2*E[k])-fun(A-H2*E[k]))/(2*H2)

    def S(A, i, j):
        return dd(lambda B: Lc(B, i), A, j)+dd(lambda B: Lc(B, j), A, i)
    # (a) d_3 F_12 = d_1 S_23 - d_2 S_13
    lhs = dd(lambda B: Fc(B, 0, 1), A0, 2)
    rhs = dd(lambda B: S(B, 1, 2), A0, 0)-dd(lambda B: S(B, 0, 2), A0, 1)
    assert abs(lhs-rhs) < 2e-3*abs(lhs) and abs(lhs) > 1e-3
    # (b) 射线定相：Phi(A) = -i int_0^1 L_{tA}(A) dt；L^rad_i = L_i - i d_i Phi。
    xs, ws = np.polynomial.legendre.leggauss(24); xs, ws = (xs+1)/2, ws/2

    def phi(A):
        return -1j*sum(w*sum(A[i]*Lc(x*A, i) for i in range(3)) for x, w in zip(xs, ws))

    def Lrad(A, i):
        return Lc(A, i)-1j*dd(phi, A, i)
    radial = sum(A0[i]*Lrad(A0, i) for i in range(3))
    assert abs(radial) < 1e-5
    # L^rad 只由曲率决定：L^rad_i(A) = int t F_{tA}(A, e_i) dt
    for i in range(3):
        pred = sum(w*x*sum(A0[k]*Fc(x*A0, k, i) for k in range(3)) for x, w in zip(xs, ws))
        assert abs(Lrad(A0, i)-pred) < 1e-5
    i, j = 0, 1
    Srad = dd(lambda B: Lrad(B, i), A0, j)+dd(lambda B: Lrad(B, j), A0, i)
    pred = sum(w*x*x*sum(A0[k]*(dd(lambda B: Fc(B, k, i), x*A0, j)+dd(lambda B: Fc(B, k, j), x*A0, i))
                         for k in range(3)) for x, w in zip(xs, ws))
    assert abs(Srad-pred) < 2e-3*abs(pred) and abs(pred) > 1e-4
    print(f"[仅数值 玩具族] 旋度 {lhs.imag:.6f}i vs {rhs.imag:.6f}i；射线分量 {abs(radial):.1e}；"
          f"S^rad {Srad.imag:.6f}i vs 公式 {pred.imag:.6f}i")


def numeric_three_link_derivative():
    """3^4 周期格、Wilson 质量 -1、U(1) 链：十五分量端点曲率
    G(e,f) = sum_c sigma_c y_c^2 g(U^{y_c}) 对第三条链角的中心差分，与零场、规范噪声比较。"""
    import numpy as np
    L = 3
    V = L**4
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    Z, Id = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z, -1j*sk], [1j*sk, Z]]) for sk in s]+[np.block([[Z, Id], [Id, Z]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    sites = list(itertools.product(range(L), repeat=4))
    idx = {x: k for k, x in enumerate(sites)}
    nb = [[idx[tuple((x[m]+(mu == m)) % L for m in range(4))] for mu in range(4)] for x in sites]
    Pm = [(np.eye(4)-g)/2 for g in gam]; Pp = [(np.eye(4)+g)/2 for g in gam]
    G5 = np.kron(np.eye(V), g5)

    def H_of(theta, y):
        D = 3*np.eye(4*V, dtype=complex)
        for i in range(V):
            for mu in range(4):
                j = nb[i][mu]; u = np.exp(1j*y*theta[i, mu])
                D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*u
                D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*np.conj(u)
        return G5@D

    def K_of(theta, y, e):
        i, mu = e; j = nb[i][mu]; u = np.exp(1j*y*theta[i, mu])
        D = np.zeros((4*V, 4*V), complex)
        D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*1j*u
        D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*(-1j)*np.conj(u)
        return G5@D

    def gfun(theta, y, e, f):
        lam, W = np.linalg.eigh(H_of(theta, y)); pos = lam > 0
        assert pos.sum() == 2*V and np.min(np.abs(lam)) > np.sqrt(0.85)
        Wp, Wm = W[:, pos], W[:, ~pos]; den = lam[pos][:, None]-lam[~pos][None, :]
        A = (Wp.conj().T@K_of(theta, y, e)@Wm)/den; B = (Wp.conj().T@K_of(theta, y, f)@Wm)/den
        return np.trace(A@B.conj().T)-np.trace(B@A.conj().T)
    left = {1: 6, -3: 2}; right = {4: 3, -2: 3, -6: 1}

    def total(theta, e, f):
        return (sum(m*y*y*gfun(theta, y, e, f) for y, m in left.items())
                - sum(m*y*y*gfun(theta, y, e, f) for y, m in right.items()))

    def max_plaq(theta):
        mx = 0
        for i in range(V):
            for mu in range(4):
                for nu in range(mu+1, 4):
                    mx = max(mx, abs(theta[i, mu]+theta[nb[i][mu], nu]-theta[nb[i][nu], mu]-theta[i, nu]))
        return mx
    rng = np.random.default_rng(20261127)
    theta = 1e-4*rng.standard_normal((V, 4))
    e, f, g = (0, 0), (idx[(0, 1, 0, 0)], 2), (idx[(0, 0, 1, 0)], 1)
    hstep = 6e-5
    tp, tm = theta.copy(), theta.copy(); tp[g] += hstep; tm[g] -= hstep
    for t in (theta, tp, tm):
        assert 2*np.sin(3*max_plaq(t)) <= 1/200
    def dg(th, hs):
        a, b = th.copy(), th.copy(); a[g] += hs; b[g] -= hs
        for t in (a, b):
            assert 2*np.sin(3*max_plaq(t)) <= 1/200
        return (total(a, e, f)-total(b, e, f))/(2*hs)
    deriv = dg(theta, hstep)
    half_step = dg(theta, hstep/2)
    half_amp = dg(theta/2, hstep)
    zero = np.zeros((V, 4)); zp, zm = zero.copy(), zero.copy(); zp[g] += hstep; zm[g] -= hstep
    noise = abs((total(zp, e, f)-total(zm, e, f))/(2*hstep))
    assert abs(deriv.real) < 1e-10 and abs(deriv) > 1e-9 and abs(deriv) > 100*noise
    assert 2.2e-9 < deriv.imag < 2.3e-9
    assert abs(half_step-deriv) < 1e-2*abs(deriv)
    assert 3.9 < deriv.imag/half_amp.imag < 4.1
    print(f"[仅数值 3^4] 十五分量端点曲率的三链导数 {deriv.imag:.4e}i（幅 1e-4），步长减半相对差 "
          f"{abs(half_step-deriv)/abs(deriv):.1e}，幅减半比 {deriv.imag/half_amp.imag:.3f}，零场同一差分 {noise:.1e}")
    return deriv


def numeric_anomaly_density():
    """3^4 周期格、Wilson 质量 -1、U(1) 链，场 A = eps*T（T 为固定随机场）。
    A(x) = sum_c sigma_c q_c m_c tr_spin P_+^{(q_c)}(x,x) 减去零场值（三个标准块的内部因子 3 略去）。
    核对：A 不恒为零且按 eps^4 标度；sum_x A 近零；随机规范变换下不变；
    射线 {t eps T : 0<=t<=1} 上 H 有隙（浮点本征值加 Weyl 扰动的格点插值），
    以便用解析性把非零结论传回 K_a 中的射线初段。"""
    import itertools
    import numpy as np
    L = 3
    V = L**4
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]])]
    Z, Id = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z, -1j*sk], [1j*sk, Z]]) for sk in s]+[np.block([[Z, Id], [Id, Z]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    sites = list(itertools.product(range(L), repeat=4))
    idx = {x: k for k, x in enumerate(sites)}
    nb = [[idx[tuple((x[m]+(mu == m)) % L for m in range(4))] for mu in range(4)] for x in sites]
    Pm = [(np.eye(4)-g)/2 for g in gam]; Pp = [(np.eye(4)+g)/2 for g in gam]
    G5 = np.kron(np.eye(V), g5)

    def H_of(theta, y):
        D = 3*np.eye(4*V, dtype=complex)
        for i in range(V):
            for mu in range(4):
                j = nb[i][mu]; u = np.exp(1j*y*theta[i, mu])
                D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*u
                D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*np.conj(u)
        return G5@D

    def diagP(theta, y):
        lam, W = np.linalg.eigh(H_of(theta, y)); pos = lam > 0
        assert pos.sum() == 2*V
        Wp = W[:, pos]
        d = np.einsum("ij,ij->i", Wp, Wp.conj()).real
        return d.reshape(V, 4).sum(1), np.min(np.abs(lam))
    weights = {1: 6, -3: 2, 4: -3, -2: -3, -6: -1}      # sigma_c * 重数
    assert sum(w*q**3 for q, w in weights.items()) == 0 and sum(w*q for q, w in weights.items()) == 0
    zero = np.zeros((V, 4))
    base = {q: diagP(zero, q)[0] for q in weights}

    def anomaly(theta):
        out, gap = np.zeros(V), np.inf
        for q, w in weights.items():
            d, g = diagP(theta, q); out += w*q*(d-base[q]); gap = min(gap, g)
        return out, gap
    T = np.random.default_rng(20261128).standard_normal((V, 4))
    res = {}
    for eps in (1e-2, 5e-3, 2.5e-3):
        a, gap = anomaly(eps*T)
        res[eps] = (np.max(np.abs(a)), a.sum(), gap, a)
    (m1, s1, g1, a1), (m2, s2, g2, a2), (m3, s3, g3, a3) = res[1e-2], res[5e-3], res[2.5e-3]
    zero_val = np.max(np.abs(anomaly(zero)[0]))
    lam = np.random.default_rng(7).uniform(-np.pi, np.pi, V)
    grad = np.array([[lam[i]-lam[nb[i][mu]] for mu in range(4)] for i in range(V)])
    gauge_diff = np.max(np.abs(anomaly(1e-2*T+grad)[0]-a1))
    # 射线 t*eps*T（eps = 1e-2）上的隙：N+1 个格点上的浮点最小 |lambda|，相邻格点间
    # ||H(t)-H(s)|| <= 8 |q| max|A| |t-s|，由 Weyl 不等式插值。
    N, A1 = 20, 1e-2*T
    lip = {q: 8*abs(q)*np.max(np.abs(A1)) for q in weights}
    seg_gap = min(min(diagP(k/N*A1, q)[1] for k in range(N+1))-lip[q]/(2*N) for q in weights)
    assert 2.10e-6 < m1 < 2.11e-6 and 1.32e-7 < m2 < 1.33e-7 and 8.26e-9 < m3 < 8.28e-9
    assert 15.5 < m1/m2 < 16.5 and 15.5 < m2/m3 < 16.5
    assert max(abs(s1), abs(s2), abs(s3)) < 1e-12 and gauge_diff < 1e-12 and zero_val < 1e-14
    assert m3 > 1000*max(abs(s3), gauge_diff, zero_val, 1e-12)
    assert seg_gap > 0.9
    assert 2.5e-3*np.max(np.abs(T)) > float(A_RADIUS)           # 所用场均不在 K_a 内
    print(f"[仅数值 3^4] 手征密度 max|A| = {m1:.4e}, {m2:.4e}, {m3:.4e}（eps = 1e-2, 5e-3, 2.5e-3），"
          f"相邻比 {m1/m2:.3f}, {m2/m3:.3f}；|sum A| <= {max(abs(s1), abs(s2), abs(s3)):.1e}；"
          f"规范变换差 {gauge_diff:.1e}；零场 {zero_val:.1e}；射线隙下界 {seg_gap:.4f}")


if __name__ == "__main__":
    out = certify()
    rf = out["cases"]["reference"]
    print(f"[三链界 参考端点] C3 < {float(out['c3_ref']):.4e}  率 b = {float(rf['b']):.8f}（按周长）  周长门槛 {out['n3_ref']}")
    print(f"[三链界 普适壳]   C3 < {float(out['c3_sh']):.4e}  率 b = 0.001523793（按周长）  周长门槛 {out['n3_sh']}")
    print(f"[非最优性] 精确 e^(6b) 下参考 C3 = {out['tight3'].str(5)}；精确 coth 下 C_S = {out['tight_s'].str(5)}")
    print(f"[射线定相 参考对称部分] 前因子 < {float(out['cs_ref']):.4e}  率 {float(out['rate_ref']):.8f}  门槛 d >= {out['ns_ref']}")
    for k, v in out["totals"].items():
        print(f"[完整测度对称部分 {k}] 前因子 < {float(v['const']):.4e}（其中输运项 {float(v['kato']):.4e}）"
              f"  率 {float(v['rate']):.9f}  门槛 d >= {v['n']}")
    numeric_toy()
    numeric_three_link_derivative()
    numeric_anomaly_density()
    print("全部断言通过")

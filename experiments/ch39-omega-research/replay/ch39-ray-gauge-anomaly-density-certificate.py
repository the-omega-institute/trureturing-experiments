"""纯超荷凸片 K_a 内一点的十五分量手征密度 A 的 Arb 区间证书。

对象（与《凸片上外尔测度流对称部分模相位Hessian由曲率导数决定……》命题 5 同一记号）：
  A(x)(A) = sum_c sigma_c q_c tr P_{0,+}^{(c)}(A)(x,x) - (同式在 A=0 的值)，
终点 P_{0,+} = (1+sign H)/2，H = gamma_5 D_W，D_W 为 Wilson 质量 -1（对角 3）的标准 Wilson 算子，
每个分量的终点在内部十九维上含三个相同的标准块，故 A = 3 * Atil，
  Atil(x) = sum_q w_q q [tr_spin P_+^{(q)}(x,x) - tr_spin P_+^{(0)}(x,x)]，
w_q = sigma*重数 = {1:6, -3:2, 4:-3, -2:-3, -6:-1}。

证书方法（全部为 Arb 球算术，精度 900 位（Bauer-Fike 一步用 128 位））：
  1. 谱区间：浮点本征基 W 换成精确二进制数后，以 Arb 计算 R = HW - W Lambda 与 G = W^*W - I，
     由 Bauer-Fike 得到 H 的全部本征值满足 lo <= |mu| <= hi。
  2. sign(H) = H (H^2)^(-1/2)；取有理 a <= lo^2 与有理 s，使 b = a((1+s)/(1-s))^2 >= hi^2，
     t^(-1/2) = kappa^(-1/2) sum_n (-s)^n P_n(u)（Legendre 生成函数，u = (t-c)/d，kappa = c/(1+s^2)），
     截断到 N 项的误差在 [a,b] 上一致不超过 kappa^(-1/2) s^(N+1)/(1-s)，乘 ||H|| <= sqrt(b) 即 sign 的算子范数误差。
  3. 对角元 e^* sign(H) e 由矩阵-向量 Legendre 递推在 Arb 中计算，再加截断误差半径。
"""
import itertools
import time
from fractions import Fraction as Fr

import numpy as np
from flint import acb, acb_mat, arb, ctx

PREC = 900
TOL_LOG10 = -30
WEIGHTS = {1: 6, -3: 2, 4: -3, -2: -3, -6: -1}      # sigma_c * 重数，十五个分量
INTERNAL_STANDARD_BLOCKS = 3
A_HALF_WIDTH = Fr(1, 4800)                           # K_a 的 a
FIELD_SCALE = Fr(1, 5000)                            # 证书场 A = T/5000，T 取 +-1，在 K_a 内部
PATTERNS = {
    (2, 2, 2, 2): "+++---+--+++-+++-+----+--+--+++-+++--++--++--++--++--++--++-+++-",
    (2, 2, 1, 1): "+-+++-+---+++--+",
}


def gammas():
    I = 1j
    s = [np.array([[0, 1], [1, 0]]), np.array([[0, -I], [I, 0]]), np.array([[1, 0], [0, -1]])]
    Z, I2 = np.zeros((2, 2)), np.eye(2)
    gam = [np.block([[Z, -I*sk], [I*sk, Z]]) for sk in s]+[np.block([[Z, I2], [I2, Z]]).astype(complex)]
    g5 = gam[0]@gam[1]@gam[2]@gam[3]
    return gam, g5


def lattice(Ls):
    sites = list(itertools.product(*[range(l) for l in Ls]))
    idx = {x: k for k, x in enumerate(sites)}
    nb = [[idx[tuple((x[m]+(mu == m)) % Ls[m] for m in range(4))] for mu in range(4)] for x in sites]
    return len(sites), nb


def half_gauss(z):
    re, im = Fr(round(2*z.real), 2), Fr(round(2*z.imag), 2)
    assert abs(complex(float(re), float(im))-z) < 1e-12
    return acb(arb(re.numerator)/re.denominator, arb(im.numerator)/im.denominator), (re, im) != (0, 0)


def build_H(Ls, field, q):
    """H_q = gamma_5 D_W(U^q)，U = exp(i A)，A 为精确有理角。返回 Arb 矩阵与其浮点中点。"""
    V, nb = lattice(Ls)
    n = 4*V
    gam, g5 = gammas()
    Pm = [(np.eye(4)-g)/2 for g in gam]
    Pp = [(np.eye(4)+g)/2 for g in gam]
    D = [[acb(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        D[i][i] = acb(3)
    for i in range(V):
        for mu in range(4):
            j = nb[i][mu]
            th = field.get((i, mu), Fr(0))*q
            ang = arb(th.numerator)/th.denominator
            u, uc = acb(ang.cos(), ang.sin()), acb(ang.cos(), -ang.sin())
            for r in range(4):
                for c in range(4):
                    pm, nz = half_gauss(Pm[mu][r, c])
                    if nz:
                        D[4*i+r][4*j+c] -= pm*u
                    pp, nz = half_gauss(Pp[mu][r, c])
                    if nz:
                        D[4*j+r][4*i+c] -= pp*uc
    # gamma_5 在此基中为对角 diag(1,1,-1,-1)
    assert np.allclose(g5, np.diag(np.diag(g5))) and set(np.round(np.diag(g5).real).astype(int)) == {1, -1}
    g5d = [int(round(v.real)) for v in np.diag(g5)]
    H = [[D[i][j]*g5d[i % 4] for j in range(n)] for i in range(n)]
    Hm = acb_mat(H)
    Hn = np.array([[complex(float(H[i][j].real.mid()), float(H[i][j].imag.mid())) for j in range(n)] for i in range(n)])
    return Hm, Hn


def norm2_upper(M):
    """谱范数上界 sqrt(||M||_1 ||M||_inf)，逐个以 arb.max 归约。"""
    n, m = M.nrows(), M.ncols()
    rows, cols = [arb(0)]*n, [arb(0)]*m
    for i in range(n):
        for j in range(m):
            v = abs(M[i, j]).upper()
            rows[i] += v
            cols[j] += v
    r = rows[0]
    for x in rows[1:]:
        r = r.max(x)
    c = cols[0]
    for x in cols[1:]:
        c = c.max(x)
    return (r*c).sqrt().upper()


def identity(n):
    return acb_mat([[acb(1) if i == j else acb(0) for j in range(n)] for i in range(n)])


def spectrum_bounds(Hm, Hn):
    """Bauer-Fike：H = W Lambda W^{-1} + R W^{-1}，||W^{-1}|| <= (1-E)^(-1/2)，kappa(W) <= ((1+E)/(1-E))^(1/2)。"""
    lam, W = np.linalg.eigh(Hn)
    n = len(lam)
    saved = ctx.prec
    ctx.prec = 128
    Wm = acb_mat([[acb(float(W[i, j].real), float(W[i, j].imag)) for j in range(n)] for i in range(n)])
    Lm = acb_mat([[acb(float(lam[i])) if i == j else acb(0) for j in range(n)] for i in range(n)])
    R = Hm*Wm - Wm*Lm
    E = norm2_upper(Wm.conjugate().transpose()*Wm - identity(n))
    assert E < 1
    rho = (norm2_upper(R)*((1+E)/(1-E)).sqrt()/(1-E).sqrt()).upper()
    lo = arb(float(np.min(np.abs(lam)))) - rho
    hi = arb(float(np.max(np.abs(lam)))) + rho
    ctx.prec = saved
    assert lo > 0
    return lo, hi, rho


def legendre_window(lo, hi):
    """取有理 a <= lo^2 与有理 s，使 b = a((1+s)/(1-s))^2 >= hi^2。"""
    lo2 = float((lo*lo).lower())
    a = Fr(int(lo2*10**6) - 1, 10**6)
    assert arb(a.numerator)/a.denominator < (lo*lo).lower()
    hi2 = float((hi*hi).upper())
    r = (hi2/float(a))**0.5
    s = Fr(int((r-1)/(r+1)*10**6) + 2, 10**6)
    b = a*((1+s)/(1-s))**2
    assert arb(b.numerator)/b.denominator > (hi*hi).upper()
    return a, s, b


def sign_diagonal(Hm, cols, a_q, s_q, tol_log10=None):
    if tol_log10 is None:
        tol_log10 = TOL_LOG10
    """对 cols 中的标准基向量 e_k，返回 e_k^* sign(H) e_k 的 Arb 球列表（已含截断误差半径）。"""
    n = Hm.nrows()
    a = arb(a_q.numerator)/a_q.denominator
    s = arb(s_q.numerator)/s_q.denominator
    bq = a_q*((1+s_q)/(1-s_q))**2
    b = arb(bq.numerator)/bq.denominator
    c, d = (a+b)/2, (b-a)/2
    kap = c/(1+s*s)
    N = 1
    while True:
        err = (b.sqrt()/kap.sqrt()*s**(N+1)/(1-s)).upper()
        if err < arb(10)**tol_log10:
            break
        N += 1
    H2 = Hm*Hm
    X = acb_mat([[(H2[i, j]-(c if i == j else 0))/d for j in range(n)] for i in range(n)])
    Ecol = acb_mat([[acb(1) if i == cols[k] else acb(0) for k in range(len(cols))] for i in range(n)])
    Q0, Q1 = Ecol, X*Ecol
    S = Q0 + Q1*(-s)
    pw = -s
    for k in range(1, N):
        Q2 = (X*Q1*(2*k+1) - Q0*k)*(arb(1)/(k+1))
        pw = pw*(-s)
        S = S + Q2*pw
        Q0, Q1 = Q1, Q2
    Y = Hm*S
    out = []
    for k, col in enumerate(cols):
        v = Y[col, k]/kap.sqrt()
        assert abs(v.imag) < arb(10)**-30       # 对角元为实数
        out.append(v.real + arb(0, err))
    return out, N, err


def field_from_pattern(Ls, pattern, scale):
    V, _ = lattice(Ls)
    assert len(pattern) == 4*V
    return {(i, mu): scale*(1 if pattern[4*i+mu] == "+" else -1) for i in range(V) for mu in range(4)}


def plaquette_max(Ls, field):
    V, nb = lattice(Ls)
    m = Fr(0)
    for i in range(V):
        for mu in range(4):
            for nu in range(mu+1, 4):
                th = field[(i, mu)] + field[(nb[i][mu], nu)] - field[(nb[i][nu], mu)] - field[(i, nu)]
                m = max(m, abs(th))
    return m


def spin_trace_P_plus(Ls, field, q, sites, log):
    """tr_spin P_+(x,x) = 2 + tr_spin sign(H)(x,x)/2 的 Arb 球，x 取 sites。"""
    Hm, Hn = build_H(Ls, field, q)
    lo, hi, rho = spectrum_bounds(Hm, Hn)
    a_q, s_q, b_q = legendre_window(lo, hi)
    cols = [4*x+k for x in sites for k in range(4)]
    diag, N, err = sign_diagonal(Hm, cols, a_q, s_q)
    tr = []
    for t, x in enumerate(sites):
        ssum = diag[4*t]+diag[4*t+1]+diag[4*t+2]+diag[4*t+3]
        tr.append(2 + ssum/2)
    log.append((q, lo, hi, rho, a_q, s_q, N, err))
    return tr, lo


def anomaly_density(Ls, field, sites, log):
    zero = {k: Fr(0) for k in field}
    base, _ = spin_trace_P_plus(Ls, zero, 1, sites, log)
    At = [arb(0) for _ in sites]
    gap = None
    for q, w in WEIGHTS.items():
        tr, lo = spin_trace_P_plus(Ls, field, q, sites, log)
        gap = lo if gap is None else gap.min(lo)
        for t in range(len(sites)):
            At[t] += w*q*(tr[t]-base[t])
    return [INTERNAL_STANDARD_BLOCKS*x for x in At], gap


def run():
    ctx.prec = PREC
    t0 = time.time()
    assert sum(w for w in WEIGHTS.values()) == 1 and sum(abs(w) for w in WEIGHTS.values()) == 15
    assert sum(w*q for q, w in WEIGHTS.items()) == 0 and sum(w*q**3 for q, w in WEIGHTS.items()) == 0
    moments = {m: sum(w*q**m for q, w in WEIGHTS.items()) for m in range(1, 6)}
    assert moments == {1: 0, 2: -72, 3: 0, 4: -1944, 5: 4320}
    print("超荷矩 sum w q^m (m=1..5):", moments)
    assert FIELD_SCALE < A_HALF_WIDTH
    # 电荷共轭 C = gamma_2 gamma_4：C gamma_mu C^{-1} = -gamma_mu^T，[C, gamma_5] = 0，C 为实正交矩阵（精确整数）
    gam, g5 = gammas()
    C = gam[1]@gam[3]
    assert np.array_equal(C.imag, np.zeros((4, 4))) and np.array_equal(C@C.T, np.eye(4))
    Cinv = C.T
    assert all(np.array_equal(C@g@Cinv, -g.T) for g in gam) and np.array_equal(C@g5, g5@C)
    # 隙引理：|spec H_q(A)| >= 1 - 4*max|q|*a = 199/200 于 K_a
    GAP = 1 - 4*max(abs(q) for q in WEIGHTS)*A_HALF_WIDTH
    assert GAP == Fr(199, 200)
    results = {}
    for Ls, pattern in PATTERNS.items():
        V, nb = lattice(Ls)
        field = field_from_pattern(Ls, pattern, FIELD_SCALE)
        pm = plaquette_max(Ls, field)
        assert max(abs(v) for v in field.values()) == FIELD_SCALE < A_HALF_WIDTH   # K_a 内部
        assert 6*pm <= Fr(1, 200)
        log = []
        sites = list(range(V))
        A, gap = anomaly_density(Ls, field, sites, log)
        # 总和恒为零（秩不变）：证书球须含 0
        tot = arb(0)
        for x in A:
            tot += x
        assert tot.contains(0)
        A0 = A[0]
        assert A0.rad() < arb(10)**-25
        assert gap > arb(GAP.numerator)/GAP.denominator
        results[Ls] = (A0, A, pm, gap, log)
        Nmax = max(e[6] for e in log)
        errmax = log[0][7]
        for e in log:
            errmax = errmax.max(e[7])
        print(f"[格 {Ls}，V={V}] 场 A=T/5000（max|A_g|=1/5000<1/4800），最大方格角 {pm}，6*方格角 <= 1/200")
        print(f"  谱区间：min|mu| >= {gap.lower().str(8)}；Legendre 项数 <= {Nmax}，截断误差 <= {errmax.str(3)}")
        print(f"  A(0) = {A0.str(15)}")
        print(f"  sum_x A(x) 球含 0：{tot.str(5)}")
    # 断言给出的数字
    A0_24 = results[(2, 2, 2, 2)][0]
    A0_22 = results[(2, 2, 1, 1)][0]
    assert results[(2, 2, 2, 2)][2] == Fr(1, 2500) and results[(2, 2, 1, 1)][2] == Fr(1, 2500)
    assert results[(2, 2, 2, 2)][3] > arb("0.999100") and results[(2, 2, 1, 1)][3] > arb("0.998330")
    assert max(e[6] for Ls_ in results for e in results[Ls_][4]) <= 247
    assert A0_24 > arb("2.0362339e-11") and A0_24 < arb("2.0362340e-11")
    assert A0_22 > arb("-2.6783766e-11") and A0_22 < arb("-2.6783765e-11")
    # 2^4 上 A 在全部 16 个格点上的最大值与最小值
    A24 = results[(2, 2, 2, 2)][1]
    nonzero = sum(1 for x in A24 if not x.contains(0))
    print(f"[2^4] 16 个格点中 A(x) 证书非零的格点数：{nonzero}；各点 A(x)：", [x.str(6) for x in A24])
    assert nonzero == 16
    assert all(abs(x) > arb("2.44e-12") for x in A24)
    # 规范协变核对：对一个规范变换 omega（值 +-1/5000、0），A 的证书球重叠
    Ls = (2, 2, 2, 2)
    V, nb = lattice(Ls)
    field = field_from_pattern(Ls, PATTERNS[Ls], FIELD_SCALE)
    om = [FIELD_SCALE*((i % 3) - 1) for i in range(V)]
    gfield = {(i, mu): field[(i, mu)] + om[nb[i][mu]] - om[i] for i in range(V) for mu in range(4)}
    assert gfield != field and 6*plaquette_max(Ls, gfield) <= Fr(1, 200)
    log = []
    Ag, _ = anomaly_density(Ls, gfield, [0], log)
    assert Ag[0].overlaps(A0_24)
    print(f"[2^4] 规范变换后 A(0) = {Ag[0].str(15)}（与原值重叠）")
    # 半场：A(T/10000) 与比值
    hfield = field_from_pattern(Ls, PATTERNS[Ls], FIELD_SCALE/2)
    log = []
    Ah, _ = anomaly_density(Ls, hfield, [0], log)
    ratio = A0_24/Ah[0]
    assert ratio > arb("15.9997") and ratio < arb("15.9998")
    print(f"[2^4] A(T/10000)(0) = {Ah[0].str(15)}，比值 A(T/5000)/A(T/10000) = {ratio.str(10)}")
    assert Ah[0] > arb("1.2726669e-12") and Ah[0] < arb("1.2726670e-12")
    # 电荷共轭的相容性核对：f(A_*) 与 f(-A_*)（单个超荷 1，格点 0）的球重叠
    nfield = {k: -v for k, v in field.items()}
    fp, _ = spin_trace_P_plus(Ls, field, 1, [0], [])
    fm, _ = spin_trace_P_plus(Ls, nfield, 1, [0], [])
    assert fp[0].overlaps(fm[0]) and fp[0].rad() < arb(10)**-25
    print(f"[2^4] f(A_*)(0) = {fp[0].str(20)}，f(-A_*)(0) = {fm[0].str(20)}（重叠）")
    float_thin_lattices()
    print(f"耗时 {time.time()-t0:.1f} s")
    return results


def float_thin_lattices():
    """仅数值：只有一个方向长度大于一的长方格上，浮点 max_x|A|/3 的读数。"""
    gam, g5 = gammas()
    Pm = [(np.eye(4)-g)/2 for g in gam]
    Pp = [(np.eye(4)+g)/2 for g in gam]
    rng = np.random.default_rng(2)
    for Ls in [(2, 1, 1, 1), (3, 1, 1, 1), (4, 1, 1, 1)]:
        V, nb = lattice(Ls)
        G5 = np.kron(np.eye(V), g5)

        def Hq(th, q):
            D = 3*np.eye(4*V, dtype=complex)
            for i in range(V):
                for mu in range(4):
                    j = nb[i][mu]
                    u = np.exp(1j*q*th[i, mu])
                    D[4*i:4*i+4, 4*j:4*j+4] -= Pm[mu]*u
                    D[4*j:4*j+4, 4*i:4*i+4] -= Pp[mu]*np.conj(u)
            return G5@D

        def at(th):
            out = np.zeros(V)
            for q, w in WEIGHTS.items():
                lam, W = np.linalg.eigh(Hq(th, q))
                Wp = W[:, lam > 0]
                out += w*q*np.einsum("ij,ij->i", Wp, Wp.conj()).real.reshape(V, 4).sum(1)
            return out
        T = rng.standard_normal((V, 4))
        a0 = at(0*T)
        vals = [float(np.max(np.abs(at(e*T)-a0))) for e in (2e-2, 1e-2)]
        print(f"[仅数值 {Ls}] 随机高斯场 eps=2e-2, 1e-2：max_x|A|/3 = {vals[0]:.2e}, {vals[1]:.2e}")
        assert max(vals) < 1e-13


if __name__ == "__main__":
    run()
    print("全部断言通过")

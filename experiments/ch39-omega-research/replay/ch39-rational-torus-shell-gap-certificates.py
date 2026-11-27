"""十九点球壳单步与两步全环面结论的精确有理证书（R2）。

说明：钟相位取圆上自由点时两步分界只剩秩15且单步三代窗包含八点二二一至十点六并含于八点二二零六四八至十点六三五六一七均有精确有理证书.md
数据：ch39-rational-torus-shell-gap-certificates.json（盒划分、阈值与有理见证点；只作为待核验的候选输入）

全程只用整数与 fractions.Fraction。不计算任何余弦、正弦或平方根的浮点值。

A. 有理化模型。原壳固定边相位 z_e=s_e(A_e-i*17*sqrt3*D_e)/sqrt(A_e^2+867 D_e^2)
   是代数数（约二十个平方类，不是有理数）。对每个荷 q，取分母 2^30 的高斯有理数 zeta_{e,q}
   代替 z_e^q，|z_e^q-zeta_{e,q}| 由 sqrt3 与 sqrt(N_e) 的整数平方根上下界精确控制。
   钟相位 U=(x,y)、V=(u,v) 是单位圆上的自由点，荷 q 的联合相位取 (u^q,v^q)。
   得到的 Lhat_q(U,V) 的系数全在 Q(i) 中；真实模型 L_{q,delta} 的每个有序本征值落在
   [lambda_j(Lhat)-eta-4 delta, lambda_j(Lhat)+eta] 内（Weyl，||E||<=||E||_F，0<=K<=4）。
B. 全环面惯性。Lhat 与三循环置换精确对易；有理基 [T|N] 把 Lhat-tau 合同为 7+12 两块；
   每块只有中心（内轨道）行列含钟相位。Haynsworth：N(块)=N(F-tau)+N(S)，S 为 1x1 或 2x2
   Schur 补。S 的元素乘 (1+t^2)(1+s^2) 后是 Cayley 参数 (t,s) 的多项式；两张图
   U=±(1-t^2+2it)/(1+t^2)，t in [-1,1] 覆盖整个圆。在盒上把 s、det S、S_11 展开成
   伯恩斯坦系数，系数同号即该多项式在闭盒上定号（逐点代数不等式，覆盖盒内全部实点）。
C. 两步秩 15、其余秩的有理见证点、单步窗端点。
D. 雷吉向量扇区不变量的平方和恒等式。
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, factorial, isqrt, lcm
from pathlib import Path
import json
import time

HERE = Path(__file__).with_name
DATA = json.loads(HERE("ch39-rational-torus-shell-gap-certificates.json").read_text())
CHARGES = (1, 2, 3, 4, 6)

# ---------------- A. 几何与有理化 ----------------
R = 867
OUTER = sorted(set(permutations((1, 5, 29))))
INNER = sorted(set(permutations((13, 13, 23))))
SITES = OUTER+INNER+sorted(z for z in product(range(30), repeat=3)
                           if sum(a*a for a in z) == R and z not in OUTER+INNER)
POLE_PT = (17, 17, 17)
POLE = SITES.index(POLE_PT)
TUNED = (0, 3, 4)                      # 软辐边的三个外端
P_BITS = 30
DMAX = Q(23, 25)/2**20/(14**2+30)      # 辐边亏额上界：delta=(23/25)2^-20/(gamma^2+30.3833...)，gamma>14，
                                       # 把偏移 30.3833... 换成 30 只会增大分式，故仍是合法上界
CMIN, CMAX = 15-DMAX/3, Q(15)          # c = tr L_0/18 = 15 - delta/3
EPS = Q(5, 10**8)                      # 大于 eta + 4 delta_max 的统一余量
ETA = Q(1, 10**8)                      # 有理化误差的算子范数上界（下面断言）


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def det3(p, x, y):
    return sum(p[i]*(x[(i+1) % 3]*y[(i+2) % 3]-x[(i+2) % 3]*y[(i+1) % 3]) for i in range(3))


def edges():
    out = []
    for i, j in combinations(range(19), 2):
        d = 2*(R-dot(SITES[i], SITES[j]))
        if d <= 804:
            A = R*(R+dot(POLE_PT, SITES[i])+dot(SITES[i], SITES[j])+dot(SITES[j], POLE_PT))
            D = det3(POLE_PT, SITES[i], SITES[j])
            out.append((i, j, d, A, D, -1 if d == 200 else 1))
    return out


EDGES = edges()


def cadd(a, b): return (a[0]+b[0], a[1]+b[1])
def csub(a, b): return (a[0]-b[0], a[1]-b[1])
def cmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cconj(a): return (a[0], -a[1])
def cneg(a): return (-a[0], -a[1])


def cinv(a):
    n = a[0]*a[0]+a[1]*a[1]
    return (a[0]/n, -a[1]/n)


ZERO, ONE = (Q(0), Q(0)), (Q(1), Q(0))


def isqrt_bounds(n, k=200):
    s = isqrt(n << (2*k))
    lo, hi = Q(s, 1 << k), Q(s+1, 1 << k)
    assert lo*lo <= n <= hi*hi
    return lo, hi


def phase_power(A, D, sgn, q):
    """z^q 的实部、虚部的有理区间，z = sgn (A - i 17 sqrt3 D)/sqrt(A^2+867 D^2)。"""
    a, b = 1, 0                        # a + b i sqrt3
    for _ in range(q):
        a, b = a*A+3*b*17*D, -a*17*D+b*A
    N2 = A*A+R*D*D
    s3 = isqrt_bounds(3)
    if q % 2 == 0:
        inv = (Q(1, N2**(q//2)),)*2
    else:
        nl, nh = isqrt_bounds(N2)
        inv = (1/(N2**((q-1)//2)*nh), 1/(N2**((q-1)//2)*nl))

    def scale(x, lohi):
        v = (x*lohi[0], x*lohi[1])
        return min(v), max(v)
    sg = sgn**q
    return scale(sg*a, inv), scale(sg*b, (s3[0]*inv[0], s3[1]*inv[1]))


def rational_phase(A, D, sgn, q):
    re, im = phase_power(A, D, sgn, q)
    zr, zi = (Q(round((iv[0]+iv[1])/2*(1 << P_BITS)), 1 << P_BITS) for iv in (re, im))
    er = max(abs(re[0]-zr), abs(re[1]-zr))
    ei = max(abs(im[0]-zi), abs(im[1]-zi))
    return (zr, zi), er*er+ei*ei


def model(q):
    """Lhat_q(U,V) = B + U C1 + V C2 + conj(U) C1^H + conj(V) C2^H；返回 2 sum_e |z_e^q-zeta|^2 的上界。"""
    B = [[ZERO]*19 for _ in range(19)]
    C1 = [[ZERO]*19 for _ in range(19)]
    C2 = [[ZERO]*19 for _ in range(19)]
    err2 = Q(0)
    for i, j, d, A, D, sgn in EDGES:
        B[i][i] = cadd(B[i][i], ONE)
        B[j][j] = cadd(B[j][j], ONE)
        if POLE in (i, j):
            val = (Q(-1), Q(0))
        else:
            z, e2 = rational_phase(A, D, sgn, q)
            err2 += e2
            val = cneg(z)
        if d in (724, 804):
            assert i < 6 and j in (6, 7, 8)
            (C1 if d == 724 else C2)[i][j] = val
        else:
            B[i][j], B[j][i] = val, cconj(val)
    return B, C1, C2, 2*err2


def rotation():
    return [SITES.index((y, z, x)) for x, y, z in SITES]


def bases():
    rot = rotation()
    seen, orbs = set(), []
    for i in range(19):
        if i not in seen:
            o = [i]
            while rot[o[-1]] != i:
                o.append(rot[o[-1]])
            seen.update(o)
            orbs.append(o)
    T, N, cT, cN = [], [], None, []
    for o in orbs:
        v = [Q(int(i in o)) for i in range(19)]
        if 6 in o:
            cT = len(T)
        T.append(v)
        if len(o) == 3:
            for a, b in ((0, 1), (1, 2)):
                w = [Q(0)]*19
                w[o[a]], w[o[b]] = Q(1), Q(-1)
                if 6 in o:
                    cN.append(len(N))
                N.append(w)
    assert len(T) == 7 and len(N) == 12 and len(orbs) == 7
    return T, N, cT, cN, orbs


def congr(Bv, M, Cv=None):
    Cv = Bv if Cv is None else Cv
    nz = [[(i, x) for i, x in enumerate(v) if x] for v in Bv]
    mz = [[(i, x) for i, x in enumerate(v) if x] for v in Cv]
    out = [[ZERO]*len(Cv) for _ in Bv]
    for a in range(len(Bv)):
        for b in range(len(Cv)):
            s = ZERO
            for i, x in nz[a]:
                for j, y in mz[b]:
                    m = M[i][j]
                    if m != ZERO:
                        s = (s[0]+x*y*m[0], s[1]+x*y*m[1])
            out[a][b] = s
    return out


def rank_rational(rows):
    A = [r[:] for r in rows]
    rk, col = 0, 0
    n = len(A[0])
    for col in range(n):
        piv = next((i for i in range(rk, len(A)) if A[i][col] != 0), None)
        if piv is None:
            continue
        A[rk], A[piv] = A[piv], A[rk]
        for i in range(len(A)):
            if i != rk and A[i][col] != 0:
                f = A[i][col]/A[rk][col]
                A[i] = [x-f*y for x, y in zip(A[i], A[rk])]
        rk += 1
    return rk


def det_rational(M):
    A = [r[:] for r in M]
    n, d = len(A), Q(1)
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k] != 0), None)
        if piv is None:
            return Q(0)
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            d = -d
        d *= A[k][k]
        for i in range(k+1, n):
            f = A[i][k]/A[k][k]
            A[i] = [a-f*b for a, b in zip(A[i], A[k])]
    return d


def block_data(q):
    B, C1, C2, e2 = model(q)
    rot = rotation()
    for M in (B, C1, C2):
        assert all(M[rot[i]][rot[j]] == M[i][j] for i in range(19) for j in range(19))
    T, N, cT, cN, orbs = bases()
    assert rank_rational(T+N) == 19
    for M in (B, C1, C2):
        assert all(x == ZERO for row in congr(T, M, N) for x in row)
        assert all(x == ZERO for row in congr(N, M, T) for x in row)
    out = []
    for Bv, cidx in ((T, [cT]), (N, cN)):
        XB, X1, X2 = congr(Bv, B), congr(Bv, C1), congr(Bv, C2)
        n = len(Bv)
        G = [[sum(x*y for x, y in zip(Bv[a], Bv[b])) for b in range(n)] for a in range(n)]
        W = [a for a in range(n) if a not in cidx]
        for X in (X1, X2):
            assert all(X[a][b] == ZERO for a in range(n) for b in range(n) if not (a in W and b in cidx))
        assert all(G[a][b] == 0 for a in W for b in cidx)
        assert all(XB[a][b] == cconj(XB[b][a]) for a in range(n) for b in range(n))
        out.append(dict(XB=XB, X1=X1, X2=X2, G=G, W=W, C=cidx, n=n))
    return out, e2


def inertia(M):
    """非奇异厄米高斯有理矩阵的负惯性指数（精确合同对角化）。"""
    A = [row[:] for row in M]
    neg, idx = 0, list(range(len(A)))
    while idx:
        k = next((i for i in idx if A[i][i][0] != 0), None)
        if k is None:
            i0 = idx[0]
            j0 = next((j for j in idx if A[i0][j] != ZERO), None)
            assert j0 is not None, "singular"
            c = A[i0][j0]              # e_i + conj(a_ij) e_j 的对角元为 2|a_ij|^2
            A[i0] = [cadd(x, cmul(c, y)) for x, y in zip(A[i0], A[j0])]
            for r in range(len(A)):
                A[r][i0] = cadd(A[r][i0], cmul(cconj(c), A[r][j0]))
            continue
        p = A[k][k]
        assert p[1] == 0
        neg += p[0] < 0
        idx.remove(k)
        ip = 1/p[0]
        for i in idx:
            if A[i][k] != ZERO:
                f = (A[i][k][0]*ip, A[i][k][1]*ip)
                for j in idx:
                    if A[k][j] != ZERO:
                        A[i][j] = csub(A[i][j], cmul(f, A[k][j]))
    return neg


def solve(F, Bcols):
    n = len(F)
    A = [F[i][:]+Bcols[i][:] for i in range(n)]
    for k in range(n):
        piv = next(i for i in range(k, n) if A[i][k] != ZERO)
        A[k], A[piv] = A[piv], A[k]
        ip = cinv(A[k][k])
        A[k] = [cmul(ip, x) for x in A[k]]
        for i in range(n):
            if i != k and A[i][k] != ZERO:
                f = A[i][k]
                A[i] = [csub(x, cmul(f, y)) for x, y in zip(A[i], A[k])]
    return [row[n:] for row in A]


def schur(blk, tau):
    XB, X1, X2, G, W, C = (blk[k] for k in ("XB", "X1", "X2", "G", "W", "C"))
    n = blk["n"]
    Y = [[(XB[a][b][0]-tau*G[a][b], XB[a][b][1]) for b in range(n)] for a in range(n)]
    F = [[Y[a][b] for b in W] for a in W]
    negF = inertia(F)
    cols = [[[M[a][c] for c in C] for a in W] for M in (Y, X1, X2)]
    sol = solve(F, [sum((col[i] for col in cols), []) for i in range(len(W))])
    m = len(C)

    def H(a, b):
        return [[sum_c(cmul(cconj(cols[a][i][r]), sol[i][b*m+s]) for i in range(len(W)))
                 for s in range(m)] for r in range(m)]
    Hs = {(a, b): H(a, b) for a in range(3) for b in range(3)}
    S = {}
    for r in range(m):
        for s in range(m):
            const = Y[C[r]][C[s]]
            for a in range(3):
                const = csub(const, Hs[a, a][r][s])
            S[r, s] = {"one": const, "U": cneg(Hs[0, 1][r][s]), "Ub": cneg(Hs[1, 0][r][s]),
                       "V": cneg(Hs[0, 2][r][s]), "Vb": cneg(Hs[2, 0][r][s]),
                       "UbV": cneg(Hs[1, 2][r][s]), "UVb": cneg(Hs[2, 1][r][s])}
    return negF, S


def schur_family(blk):
    """S(tau) 对 tau 的精确有理依赖：det F(tau) 次数 <= |W|，B_a^H adj F(tau) B_b 次数 <= |W|-1，
    在 |W|+1 个节点 tau=-1,-2,... 上精确插值（牛顿差商），此后按 tau 求值。"""
    W, C, n = blk["W"], blk["C"], len(blk["W"])
    nodes = [Q(-1-k) for k in range(n+1)]
    vals = []
    for x in nodes:
        negF, S = schur(blk, x)
        assert negF == 0
        F = [[(blk["XB"][a][b][0]-x*blk["G"][a][b], blk["XB"][a][b][1]) for b in W] for a in W]
        dF = det_complex(F)
        assert dF[1] == 0 and dF[0] > 0
        # S = Ycc - sum H；记录 det F * (Ycc - S) 的各单项式系数（多项式，次数 <= n-1 或 n）
        rec = {}
        for key, e in S.items():
            ycc = (blk["XB"][C[key[0]]][C[key[1]]][0]-x*blk["G"][C[key[0]]][C[key[1]]],
                   blk["XB"][C[key[0]]][C[key[1]]][1])
            for mn, c in e.items():
                base = ycc if mn == "one" else ZERO
                rec[key, mn] = cmul((dF[0], Q(0)), csub(base, c))
        rec["det"] = dF
        vals.append(rec)
    keys = list(vals[0])

    def newton(ys):
        coef = list(ys)
        for j in range(1, len(nodes)):
            for i in range(len(nodes)-1, j-1, -1):
                d = nodes[i]-nodes[i-j]
                coef[i] = ((coef[i][0]-coef[i-1][0])/d, (coef[i][1]-coef[i-1][1])/d)
        return coef
    table = {k: newton([v[k] for v in vals]) for k in keys}

    def at(tau):
        out = {}
        for k, coef in table.items():
            acc = coef[-1]
            for i in range(len(coef)-2, -1, -1):
                acc = cadd(coef[i], cmul((tau-nodes[i], Q(0)), acc))
            out[k] = acc
        return out
    return at


def det_complex(M):
    A = [r[:] for r in M]
    n, d = len(A), ONE
    for k in range(n):
        piv = next((i for i in range(k, n) if A[i][k] != ZERO), None)
        if piv is None:
            return ZERO
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            d = cneg(d)
        d = cmul(d, A[k][k])
        ip = cinv(A[k][k])
        for i in range(k+1, n):
            if A[i][k] != ZERO:
                f = cmul(A[i][k], ip)
                A[i] = [csub(a, cmul(f, b)) for a, b in zip(A[i], A[k])]
    return d


def schur_fast(blk, fam, tau):
    """与 schur(blk, tau) 相同的结果：N(F-tau) 由精确合同对角化给出，S 由插值族求值。"""
    W = blk["W"]
    F = [[(blk["XB"][a][b][0]-tau*blk["G"][a][b], blk["XB"][a][b][1]) for b in W] for a in W]
    negF = inertia(F)
    v = fam(tau)
    dF = v["det"]
    assert dF[1] == 0 and dF[0] != 0
    inv = Q(1)/dF[0]
    S = {}
    for k, c in v.items():
        if k == "det":
            continue
        key, mn = k
        S.setdefault(key, {})[mn] = (c[0]*inv, c[1]*inv)
    # det F * (Ycc - S) 的 'one' 项含 Ycc；其余单项式为 det F * (-c)，还原 S：
    for key, e in S.items():
        for mn in list(e):
            ycc = ((blk["XB"][blk["C"][key[0]]][blk["C"][key[1]]][0]-tau*blk["G"][blk["C"][key[0]]][blk["C"][key[1]]],
                    blk["XB"][blk["C"][key[0]]][blk["C"][key[1]]][1]) if mn == "one" else ZERO)
            e[mn] = csub(ycc, e[mn])
    return negF, S


def sum_c(it):
    s = ZERO
    for x in it:
        s = cadd(s, x)
    return s


def padd(p, r):
    o = dict(p)
    for k, v in r.items():
        o[k] = cadd(o.get(k, ZERO), v)
    return {k: v for k, v in o.items() if v != ZERO}


def pmul(p, r):
    o = {}
    for (i, j), a in p.items():
        for (k, l), b in r.items():
            o[i+k, j+l] = cadd(o.get((i+k, j+l), ZERO), cmul(a, b))
    return {k: v for k, v in o.items() if v != ZERO}


def pscale(p, c):
    return {k: w for k, v in p.items() if (w := cmul(c, v)) != ZERO}


def pconj(p):
    return {k: cconj(v) for k, v in p.items()}


def chart_monomials(s1, s2):
    """(1+t^2)(1+s^2) 乘各单项式，U=s1(1-t^2+2it)/(1+t^2)，V=s2(1-s^2+2is)/(1+s^2)。"""
    one_t = {(0, 0): ONE, (2, 0): ONE}
    one_s = {(0, 0): ONE, (0, 2): ONE}
    U = {(0, 0): (Q(s1), Q(0)), (2, 0): (Q(-s1), Q(0)), (1, 0): (Q(0), Q(2*s1))}
    V = {(0, 0): (Q(s2), Q(0)), (0, 2): (Q(-s2), Q(0)), (0, 1): (Q(0), Q(2*s2))}
    Ub, Vb = pconj(U), pconj(V)
    return {"one": pmul(one_t, one_s), "U": pmul(U, one_s), "Ub": pmul(Ub, one_s),
            "V": pmul(one_t, V), "Vb": pmul(one_t, Vb), "UbV": pmul(Ub, V), "UVb": pmul(U, Vb)}


def entry_poly(e, mons):
    p = {}
    for mn, c in e.items():
        if c != ZERO:
            p = padd(p, pscale(mons[mn], c))
    return p


def realpart(p):
    assert all(v[1] == 0 for v in p.values())
    return {k: v[0] for k, v in p.items() if v[0] != 0}


def sign_polys(S, s1, s2):
    mons = chart_monomials(s1, s2)
    if len(S) == 1:
        p = realpart(entry_poly(S[0, 0], mons))
        assert max(max(k) for k in p) <= 2          # 次数 <= (2,2)
        return [p]
    a, b, d = (entry_poly(S[k], mons) for k in ((0, 0), (0, 1), (1, 1)))
    det = realpart(padd(pmul(a, d), pscale(pmul(b, pconj(b)), (Q(-1), Q(0)))))
    assert max(max(k) for k in det) <= 4            # 次数 <= (4,4)
    return [det, realpart(a)]


def to_int(p):
    """乘以正的公分母，得到同号的整系数多项式 (nx, ny, M)。"""
    D = 1
    for v in p.values():
        D = lcm(D, v.denominator)
    nx, ny = max(i for i, j in p), max(j for i, j in p)
    M = [[0]*(ny+1) for _ in range(nx+1)]
    for (i, j), v in p.items():
        M[i][j] = v.numerator*(D//v.denominator)
    return nx, ny, M


def bernstein(P, L, al, ga):
    """二进盒 t in [A/2^L,(A+2)/2^L]、s in [C/2^L,(C+2)/2^L]（A=2al-2^L，C=2ga-2^L）上
    伯恩斯坦系数的正整数倍：代入 t=(A+2x)/2^L 后乘 2^{L(nx+ny)}，再乘 nx! ny!。"""
    nx, ny, M = P
    A, C = 2*al-(1 << L), 2*ga-(1 << L)
    T = [[0]*(ny+1) for _ in range(nx+1)]
    for j in range(ny+1):
        for i in range(nx+1):
            c = M[i][j]
            if c:
                c <<= L*(nx-i)
                for k in range(i+1):
                    T[k][j] += c*comb(i, k)*A**(i-k)*(1 << k)
    S = [[0]*(ny+1) for _ in range(nx+1)]
    for i in range(nx+1):
        for j in range(ny+1):
            c = T[i][j]
            if c:
                c <<= L*(ny-j)
                for k in range(j+1):
                    S[i][k] += c*comb(j, k)*C**(j-k)*(1 << k)
    wx = [factorial(i)*factorial(nx-i) for i in range(nx+1)]
    wy = [factorial(j)*factorial(ny-j) for j in range(ny+1)]
    return [sum(comb(k, i)*wx[i]*comb(l, j)*wy[j]*S[i][j] for i in range(k+1) for j in range(l+1))
            for k in range(nx+1) for l in range(ny+1)]


def certify_sign(P, box, depth=0, maxdepth=12):
    """伯恩斯坦系数同号（必要时四分）即该多项式在闭盒上定号；返回符号。"""
    L, al, ga = box
    bs = bernstein(P, L, al, ga)
    if all(x > 0 for x in bs):
        return 1
    if all(x < 0 for x in bs):
        return -1
    assert depth < maxdepth, ("undecided", box)
    sg = {certify_sign(P, (L+1, 2*al+dx, 2*ga+dy), depth+1, maxdepth) for dy in (0, 1) for dx in (0, 1)}
    assert len(sg) == 1
    return sg.pop()


FULL = (0, 0, 0)                       # 整张图 [-1,1]^2
CHARTS = ((1, 1), (1, -1), (-1, 1), (-1, -1))
BLOCKS, FAMILIES, CACHE = {}, {}, {}


def schur_polys(q, tau, chart):
    key = (q, tau, chart)
    if key not in CACHE:
        CACHE[key] = [(nf, len(S), [to_int(p) for p in sign_polys(S, *chart)])
                      for nf, S in (schur_fast(b, f, tau) for b, f in zip(BLOCKS[q], FAMILIES[q]))]
    return CACHE[key]


def count_on_box(q, tau, chart, box):
    """Lhat_q - tau 在图 chart 的闭盒 box 上每一点的负惯性指数（证明其为常数）。"""
    tot = 0
    for nf, m, ps in schur_polys(q, tau, chart):
        s0 = certify_sign(ps[0], box)
        if m == 1:                      # 1x1：S 本身定号
            neg = int(s0 < 0)
        elif s0 < 0:                    # 2x2：det S < 0，惯性 (1,1)
            neg = 1
        else:                           # det S > 0：两本征值与 S_11 同号
            neg = 0 if certify_sign(ps[1], box) > 0 else 2
        tot += nf+neg
    return tot


def cayley(t, sigma):
    d = 1+t*t
    return (sigma*(1-t*t)/d, sigma*2*t/d)


def point(p):
    return cayley(Q(p[1], p[2]), p[0])


def cpow(z, k):
    r = ONE
    for _ in range(k):
        r = cmul(r, z)
    return r


def point_matrix(blk, U, V, tau):
    n, Ub, Vb = blk["n"], cconj(U), cconj(V)
    Y = []
    for a in range(n):
        row = []
        for b in range(n):
            x = (blk["XB"][a][b][0]-tau*blk["G"][a][b], blk["XB"][a][b][1])
            x = cadd(x, cmul(U, blk["X1"][a][b]))
            x = cadd(x, cmul(V, blk["X2"][a][b]))
            x = cadd(x, cmul(Ub, cconj(blk["X1"][b][a])))
            x = cadd(x, cmul(Vb, cconj(blk["X2"][b][a])))
            row.append(x)
        Y.append(row)
    return Y


def count_at(q, U, V, tau):
    """#{lambda_j(Lhat_q(U,V)) < tau}，要求 tau 不是本征值。"""
    return sum(inertia(point_matrix(b, U, V, tau)) for b in BLOCKS[q])


def qsqrt(x, up, k=60):
    s = isqrt(x.numerator*x.denominator << (2*k))
    r = Q(s+(1 if up else 0), x.denominator << k)
    assert (r*r >= x) if up else (r*r <= x)
    return r


OFF = Q(1, 10**9)


def le_cover(theta):
    """覆盖 {lambda>=0: |lambda(lambda-c)|<=theta}（对一切 c in [CMIN,CMAX]）的有理区间。"""
    a3 = (CMAX+qsqrt(CMAX**2+4*theta, True))/2+OFF
    assert a3 > CMAX and a3*(a3-CMAX) > theta
    if theta >= CMAX**2/4:
        return [(Q(0), a3)]
    assert theta < CMIN**2/4
    r = qsqrt(CMIN**2-4*theta, False)
    a1, a2 = (CMIN-r)/2+OFF, (CMAX+r)/2-OFF
    assert a1 < CMIN/2 and a1*(CMIN-a1) > theta
    assert a2 > CMAX/2 and a2*(CMIN-a2) > theta and a2 < a3
    return [(Q(0), a1), (a2, a3)]


def lt_inside(theta):
    """含于 {lambda>=0: |lambda(lambda-c)|<theta}（对一切 c in [CMIN,CMAX]）的有理区间。"""
    b3 = (CMIN+qsqrt(CMIN**2+4*theta, False))/2-OFF
    assert b3 > CMAX and b3*(b3-CMIN) < theta
    if theta > CMAX**2/4:
        return [(Q(0), b3)]
    assert theta < CMIN**2/4
    r = qsqrt(CMAX**2-4*theta, True)
    b1, b2 = (CMIN-r)/2-OFF, (CMAX+r)/2+OFF
    assert b1 <= CMIN/2 and b1*(CMAX-b1) < theta
    assert b2 >= CMAX/2 and b2*(CMAX-b2) < theta and b2 < b3
    return [(Q(0), b1), (b2, b3)]


def h_le_upper(q, U, V, theta):
    """真实模型 #{j: h_j<=theta} 的上界（一切 delta）。"""
    iv = le_cover(theta)
    tot = count_at(q, U, V, iv[0][1]+EPS)
    for a, b in iv[1:]:
        tot += count_at(q, U, V, b+EPS)-count_at(q, U, V, a-EPS)
    return tot


def h_lt_lower(q, U, V, theta):
    """真实模型 #{j: h_j<theta} 的下界（一切 delta；用到 L_delta 半正定）。"""
    iv = lt_inside(theta)
    tot = count_at(q, U, V, iv[0][1]-EPS)
    for a, b in iv[1:]:
        tot += count_at(q, U, V, b-EPS)-count_at(q, U, V, a+EPS)
    return tot


def dec(x, k):
    """有理数 x 向零截断到 k 位小数的字符串（整数运算）。"""
    sgn = "-" if x < 0 else ""
    n = abs(x.numerator)*10**k//x.denominator
    return sgn+str(n//10**k)+"."+str(n % 10**k).zfill(k)


def frac(x):
    return Q(x[0], x[1])


def covers(paths):
    """四叉树叶子两两不互为前缀，且铺满整个图。"""
    ps = set(paths)
    assert len(ps) == len(paths)
    assert not any(p[:k] in ps for p in ps for k in range(len(p)))

    def full(prefix):
        if prefix in ps:
            return True
        if not any(p.startswith(prefix) for p in ps):
            return False
        return all(full(prefix+str(k)) for k in range(4))
    return full("")


def path_box(path):
    """四叉树路径 -> (L, al, ga)：子盒 0,1,2,3 依次为左下、右下、左上、右上。"""
    L, al, ga = 0, 0, 0
    for ch in path:
        k = int(ch)
        L, al, ga = L+1, 2*al+(k & 1), 2*ga+(k >> 1)
    return L, al, ga


def verify():
    t0 = time.monotonic_ns()
    # ---- A ----
    assert len(EDGES) == 135 and POLE == 14
    deg, soft = [0]*19, [0]*19
    for i, j, *_ in EDGES:
        deg[i] += 1
        deg[j] += 1
        if POLE in (i, j) and (i in TUNED or j in TUNED):
            soft[i] += 1
            soft[j] += 1
    prof = {}
    for key in zip(deg, soft):
        prof[key] = prof.get(key, 0)+1
    assert prof == {(10, 0): 3, (10, 1): 3, (14, 0): 3, (16, 0): 6, (18, 0): 3, (18, 3): 1}
    assert 4*DMAX < Q(1, 10**7)
    # 固定边相位不是有理点：117 条非辐边中 114 条的 N_e=A_e^2+867 D_e^2 不是平方数，
    # 它们按“N_a N_b 为平方数”分成 20 个平方类（z_e 的实部 A_e/sqrt(N_e) 落在 20 个不同的二次域里）。
    Ns = [A*A+R*D*D for i, j, d, A, D, sg in EDGES if POLE not in (i, j)]
    assert len(Ns) == 117
    nonsq = [n for n in Ns if isqrt(n)**2 != n]
    assert len(nonsq) == 114
    classes = []
    for n in nonsq:
        if not any(isqrt(n*m)**2 == n*m for m in classes):
            classes.append(n)
    assert len(classes) == 20
    for q in CHARGES:
        BLOCKS[q], e2 = block_data(q)
        FAMILIES[q] = [schur_family(b) for b in BLOCKS[q]]
        for b, fam in zip(BLOCKS[q], FAMILIES[q]):      # 插值族与直接 Schur 补在一个 tau 上逐项一致
            assert schur_fast(b, fam, Q(1829812, 10**5)) == schur(b, Q(1829812, 10**5))
        assert e2 < ETA**2          # ||E||_2 <= ||E||_F <= sqrt(e2) < ETA
    assert ETA+4*DMAX < EPS
    # 软化项 K 是三条辐边 pole-{0,3,4} 的星图拉普拉斯：det(xI-K) 是首一四次多项式，
    # 在五个整数点上等于 x(x-1)^2(x-4)，故二者相等，K 的本征值为 0,1,1,4。
    spokes = [e for e in EDGES if POLE in e[:2] and ({e[0], e[1]}-{POLE}) <= set(TUNED)]
    assert len(spokes) == 3
    # K 按顶点次序 (极点, 三个软辐边外端) 写出；上一行断言软辐边恰为这 3 条，故 K 就是它们的星图拉普拉斯。
    K = [[3, -1, -1, -1], [-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]]
    for x in range(5):
        M = [[Q((x if i == j else 0)-K[i][j]) for j in range(4)] for i in range(4)]
        assert det_rational(M) == x*(x-1)**2*(x-4)
    print("A 有理化：分母 2^30，max_q ||L_q-Lhat_q|| < 10^-8；三循环对易，[T|N] 可逆，7+12 分块精确")

    # ---- B. 两步秩 15 ----
    d15 = DATA["rank15"]
    tA, tB = frac(d15["tauA"]), frac(d15["tauB"])
    assert tA == Q(1829812, 10**5) and tB == Q(1837162, 10**5)
    for q in CHARGES:
        for tau in (tA, tB):
            counts = {count_on_box(q, tau, ch, FULL) for ch in CHARTS}
            assert counts == {15}, (q, tau, counts)
    LA, LB = tA+ETA, tB-ETA-4*DMAX     # 真实模型：lambda_15 < LA，lambda_16 > LB
    h15, h16 = LA*(LA-CMIN), LB*(LB-CMAX)
    assert LA == Q("18.29812001") and LB > Q("18.37161997")
    assert CMAX**2/4 < h15 < Q("60.3494") and h16 > Q("61.94212") and h15 < h16
    w = d15["witness15"]
    # 见证点一律是有理联合相位 (u,v)，荷 q 读 (u^q, v^q)：既是有理联合点，也是荷 q 单荷环面上的有理点。
    U, V = cpow(point(w["u"]), 3), cpow(point(w["v"]), 3)
    lo15 = frac(w["tau"])
    assert w["u"] == [1, 241439, 10**6] and w["v"] == [1, 4265321, 5*10**6]
    assert w["q"] == 3 and lo15 == Q(1829811, 10**5) and count_at(3, U, V, lo15) == 14
    w = d15["witness16"]
    U, V = cpow(point(w["u"]), 2), cpow(point(w["v"]), 2)
    hi16 = frac(w["tau"])
    assert point(w["u"]) == (Q(0), Q(1)) and U == (Q(-1), Q(0)) and V == ONE and w["q"] == 2   # u=i, v=1
    assert hi16 == Q(1837164, 10**5) and count_at(2, U, V, hi16) == 16
    sup15 = (lo15-ETA-4*DMAX, LA)      # 荷三的 sup lambda_15
    inf16 = (LB, hi16+ETA)             # 荷二的 inf lambda_16
    hs15 = sup15[0]*(sup15[0]-CMAX)
    hi16h = inf16[1]*(inf16[1]-CMIN)
    assert Q("60.34917") < hs15 and hi16h < Q("61.94256")
    assert Q("60.3494")-Q("60.34917") < Q(3, 10**4) and Q("61.94256")-Q("61.94212") < Q(5, 10**4)
    print("B 全环面：每荷每点 lambda_15 <", dec(LA, 8), "，lambda_16 >", dec(LB, 8),
          "；max_q h_(15) < 60.3494 < 61.94212 < min_q h_(16)")
    print("B 端点：sup lambda_{3,15} in (18.29811-3e-8, 18.29812+1e-8)，inf lambda_{2,16} in (18.37162-3e-8,"
          " 18.37164+1e-8)；sup h_(15) in (60.34917, 60.3494)，inf h_(16) in (61.94212, 61.94256)")

    # ---- C. 其余秩与 U=V=1 的秩 3 ----
    TABLE = [
        (1, -1, '-37/250', -1, '37/250', 2, 1, '2.0307', '0.0414'),
        (2, 1, '49/1000', -1, '-49/1000', 6, 1, '7.5118', '0.1598'),
        (3, -1, '49/1000', -1, '-303/1000', 4, 1, '16.0202', '6.9541'),
        (4, -1, '1/4', -1, '1/4', 6, 1, '28.6916', '6.8466'),
        (5, 1, '0', 1, '0', 2, 1, '33.1803', '15.9774'),
        (6, 1, '453/500', 1, '821/1000', 4, 3, '34.2219', '30.7539'),
        (7, 1, '-199/1000', 1, '0', 3, 1, '37.6205', '35.5991'),
        (8, -1, '107/200', -1, '599/1000', 1, 3, '43.9992', '37.112'),
        (9, 1, '303/1000', 1, '107/200', 3, 1, '45.2377', '38.5244'),
        (10, 1, '49/500', -1, '49/500', 1, 2, '48.7046', '44.006'),
        (11, 1, '-107/200', -1, '-107/200', 3, 2, '49.6342', '44.0084'),
        (12, -1, '107/200', -1, '-599/1000', 1, 3, '56.1916', '44.0567'),
        (13, -1, '473/1000', -1, '-473/1000', 1, 6, '56.2493', '54.8667'),
        (14, -1, '107/200', -1, '-599/1000', 1, 3, '56.5552', '56.1759'),
        (16, 1, '107/200', -1, '107/200', 1, 6, '77.8866', '74.022'),
        (17, 1, '-107/200', -1, '107/200', 4, 6, '93.1066', '76.9436'),
        (18, -1, '1', -1, '1', 2, 4, '101.8391', '80.0587'),
    ]
    got = [(it["r"], it["u"][0], str(Q(it["u"][1], it["u"][2])), it["v"][0], str(Q(it["v"][1], it["v"][2])),
            it["q1"], it["q2"], frac(it["hi"]), frac(it["lo"])) for it in DATA["ranks"]]
    assert got == [row[:7]+(Q(row[7]), Q(row[8])) for row in TABLE]      # 与说明中的见证表逐项一致
    for item in DATA["ranks"]:
        r, q1, q2 = item["r"], item["q1"], item["q2"]
        u, v = point(item["u"]), point(item["v"])
        th_hi, th_lo = frac(item["hi"]), frac(item["lo"])
        assert h_le_upper(q1, cpow(u, q1), cpow(v, q1), th_hi) <= r-1   # h_{q1,(r)} > th_hi
        assert h_lt_lower(q2, cpow(u, q2), cpow(v, q2), th_lo) >= r+1   # h_{q2,(r+1)} < th_lo
        item["margin"] = th_hi-th_lo
    ranks = {it["r"]: it["margin"] for it in DATA["ranks"]}
    assert sorted(ranks) == [r for r in range(1, 19) if r != 15]
    assert ranks[14] > Q("0.379") and all(m > Q("1.382") for r, m in ranks.items() if r != 14)
    assert [dec(ranks[r], 4) for r in sorted(ranks)] == [
        "1.9893", "7.3520", "9.0661", "21.8450", "17.2029", "3.4680", "2.0214", "6.8872", "6.7133",
        "4.6986", "5.6258", "12.1349", "1.3826", "0.3793", "3.8646", "16.1630", "21.7804"]
    for q in CHARGES:
        assert h_lt_lower(q, ONE, ONE, Q("7.53559")) >= 3
        assert h_le_upper(q, ONE, ONE, Q("10.6225")) <= 3
    print("C 秩 r in {1..18}\\{15} 各在一个有理联合点上失效，余量", {r: dec(m, 4) for r, m in sorted(ranks.items())})
    print("C U=V=1：max_q h_(3) < 7.53559 < 10.6225 < min_q h_(4)")

    # ---- D. 单步窗端点 ----
    for name in ("a23", "a47"):
        dat = DATA[name]
        theta = frac(dat["theta"])
        per = {}
        for leaf in dat["leaves"]:
            per.setdefault((leaf["q"], tuple(leaf["chart"])), []).append(leaf)
        assert sorted(per) == sorted((q, ch) for q in CHARGES for ch in CHARTS)
        for (q, ch), leaves in per.items():
            assert covers([lf["path"] for lf in leaves])
            for lf in leaves:
                box = path_box(lf["path"])
                taus = [frac(t) for t in lf["taus"]]
                if name == "a23":       # lambda_2<tau_0，lambda_3<tau_1
                    assert len(taus) == 2 and sum(taus) <= 2*theta-Q(1, 10**6)
                    assert count_on_box(q, taus[0], ch, box) >= 2
                    assert count_on_box(q, taus[1], ch, box) >= 3
                else:                   # lambda_j>=tau_{j-4}，j=4..7
                    assert len(taus) == 4 and sum(taus) >= 4*theta+Q(1, 10**6)
                    for k, t in enumerate(taus):
                        assert count_on_box(q, t, ch, box) <= 3+k
        dat["n"] = len(dat["leaves"])
    th23, th47 = frac(DATA["a23"]["theta"]), frac(DATA["a47"]["theta"])
    assert th23 == Q("8.221") and th47 == Q("10.6")
    assert DATA["a23"]["n"] == 506 and DATA["a47"]["n"] == 4874
    # 真实模型：a23 < th23 - 1/(2*10^6) + ETA < th23；a47 > th47 + 1/(4*10^6) - ETA - 4 DMAX > th47。
    assert ETA < Q(1, 2*10**6) and ETA+4*DMAX < Q(1, 4*10**6)
    w = DATA["a23"]["witness"]
    U, V = cpow(point(w["u"]), 2), cpow(point(w["v"]), 2)
    lo2, lo3 = frac(w["lo"][0]), frac(w["lo"][1])
    assert w["u"] == [1, 532343, 10**7] and w["v"] == [-1, -66433, 1250000]
    assert lo2 == Q("7.9983099") and lo3 == Q("8.4429877")
    assert w["q"] == 2 and count_at(2, U, V, lo2) <= 1 and count_at(2, U, V, lo3) <= 2
    a23w = (lo2+lo3)/2-ETA-4*DMAX
    w = DATA["a47"]["witness"]
    U, V = cpow(point(w["u"]), 6), cpow(point(w["v"]), 6)
    his = [frac(x) for x in w["hi"]]
    assert w["u"] == [1, 3390803, 10**7] and w["v"] == [1, 2022581, 10**7]
    assert his == [Q("10.0716827"), Q("10.599381"), Q("10.8499148"), Q("11.021486")]
    assert w["q"] == 6 and all(count_at(6, U, V, h) >= 4+k for k, h in enumerate(his))
    a47w = sum(his)/4+ETA
    assert a23w > Q("8.220648") and a47w < Q("10.635617")
    assert Q("8.221")-Q("8.220648") < Q(36, 10**5) and Q("10.635617")-Q("10.6") < Q(357, 10**4)
    print("D 单步窗：sup a23 in (8.220648, 8.221)，inf a47 in (10.6, 10.635617)；叶子数",
          DATA["a23"]["n"], DATA["a47"]["n"])

    # ---- E. 雷吉向量扇区不变量 ----
    # c=(4(X^2+X^-2)+22+(X-X^-1)(Y-Y^-1))/50，X=x+iy，Y=u+iv，在圆上 X^-1=x-iy，Y^-1=u-iv。
    # 以 (x,y,u,v) 的整系数多项式逐项展开，验证恒等式（不是取样检验）。
    def P(**kw):
        return {tuple(kw.get(k, 0) for k in "xyuv"): Q(1)}

    def ppadd(*ps):
        o = {}
        for p in ps:
            for k, c in p.items():
                o[k] = o.get(k, 0)+c
        return {k: c for k, c in o.items() if c}

    def ppmul(p, r):
        o = {}
        for k1, c1 in p.items():
            for k2, c2 in r.items():
                k = tuple(a+b for a, b in zip(k1, k2))
                o[k] = o.get(k, 0)+c1*c2
        return {k: c for k, c in o.items() if c}

    def sc(c, p):
        return {k: c*v for k, v in p.items()}
    one = {(0, 0, 0, 0): Q(1)}
    x, y, u, v = P(x=1), P(y=1), P(u=1), P(v=1)
    # 复数用 (实部, 虚部) 两个实多项式表示
    def zmul(a, b):
        return (ppadd(ppmul(a[0], b[0]), sc(-1, ppmul(a[1], b[1]))), ppadd(ppmul(a[0], b[1]), ppmul(a[1], b[0])))
    X, Xi, Y, Yi = (x, y), (x, sc(-1, y)), (u, v), (u, sc(-1, v))
    for sgn in (1, -1):            # (+,-) 扇区 c(t,s)，(-,+) 扇区 c(t,-s)
        YY = (Y[0], sc(sgn, Y[1]))
        YYi = (Yi[0], sc(sgn, Yi[1]))
        sq = zmul(X, X)
        sqi = zmul(Xi, Xi)
        dx = (ppadd(X[0], sc(-1, Xi[0])), ppadd(X[1], sc(-1, Xi[1])))
        dy = (ppadd(YY[0], sc(-1, YYi[0])), ppadd(YY[1], sc(-1, YYi[1])))
        cross = zmul(dx, dy)
        re = ppadd(sc(4, sq[0]), sc(4, sqi[0]), sc(22, one), cross[0])
        im = ppadd(sc(4, sq[1]), sc(4, sqi[1]), cross[1])
        assert im == {}
        lhs = ppadd(sc(Q(25, 50), re), sc(-5, one))
        ymv = ppadd(y, sc(-sgn, v))
        rhs = ppadd(sc(9, ppmul(x, x)), ppmul(u, u), ppmul(ymv, ymv),
                    sc(5, ppadd(one, sc(-1, ppmul(x, x)), sc(-1, ppmul(y, y)))),
                    ppadd(one, sc(-1, ppmul(u, u)), sc(-1, ppmul(v, v))))
        assert lhs == rhs, sgn
    print("E 25c-5 = 9x^2+u^2+(y-v)^2 + 5(1-x^2-y^2) + (1-u^2-v^2)；c(t,-s) 时 (y-v) 换为 (y+v)")
    print("用时", (time.monotonic_ns()-t0)//10**9, "秒")


if __name__ == "__main__":
    verify()

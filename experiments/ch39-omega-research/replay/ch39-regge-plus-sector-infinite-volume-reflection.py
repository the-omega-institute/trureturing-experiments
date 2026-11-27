"""正偏振扇区两种格点共形围道旋转在对角时间方向无限体积中的站点反射正性。

说明见同目录《正偏振扇区两种格点共形围道旋转在对角时间方向取无限体积后对所检横向动量仍不反射正且迹旋转下的负期望有与对角时间方向盒长无关的上界》。

几何与约定全部沿用 ch39-regge-plus-sector-conformal-rotation.py（投影动量族上的 (+,+) 扇区，
对角时间 n、横向坐标 m、X=e^{it}、Y=e^{is}，K=-H，站点反射 π∈{(13)(24),(14)(23)}）。
本脚本把 n 方向取为整条 Z（m 方向仍是盒，或固定横向动量 s），并对 n 方向周期 N 的盒给出与 N 无关的界。

固定 s 时，局域实观测量 f 的反射期望密度
    E(s) = (1/2π) ∫_0^{2π} I(t,s) dt,    I = a(X,Y) Ĉ(X,Y) b(X,Y),
a 为 Θf 的扇区符号，b 为 f 的共轭符号；n 方向周期 N 的盒上同一量是梯形和
    E_N(s) = (1/N) Σ_{j<N} I(2πj/N, s)。
规范不变观测量上
    迹旋转  Ĉ = K^+ - 2 w w^#/(w^# K w)，
    谱旋转  Ĉ = K^+ - 2 u u^#/(λ u^# D u)，K u = λ D u，λ<0。
a K^+ b 化为精确有理函数 P/Q'（Q' 与 Q 只差常数与单项式，Q 为规范不变行列式的三角多项式）。
"""
from fractions import Fraction as Fr
from math import lcm
from pathlib import Path
from runpy import run_path
import time

import numpy as np
import flint
from flint import acb, acb_mat, acb_poly, arb, ctx, fmpq

ROOT = Path(__file__).parent
PS = run_path(str(ROOT/"ch39-regge-plus-sector-conformal-rotation.py"))
MASKS, BASIS, DD = PS["MASKS"], PS["BASIS"], PS["DD"]
step, perm_mask, SITE_REFL, QCOEF = PS["step"], PS["perm_mask"], PS["SITE_REFL"], PS["QCOEF"]
ctx.prec = 128
BM = [[0]*8 for _ in range(15)]
for _j, (_, _o) in enumerate(BASIS):
    for _k, _v in _o.items():
        BM[MASKS.index(_k)][_j] = int(_v)


# ------------------------------------------------------------ Laurent 代数 ----
def mul(a, b):
    o = {}
    for (p, q), c in a.items():
        for (r, s), e in b.items():
            o[(p+r, q+s)] = o.get((p+r, q+s), 0)+c*e
    return {k: v for k, v in o.items() if v}


def add(a, b, sg=1):
    o = dict(a)
    for k, v in b.items():
        o[k] = o.get(k, 0)+sg*v
    return {k: v for k, v in o.items() if v}


def conj(a):
    return {(-p, -q): c for (p, q), c in a.items()}


def sect(terms, sg=1):
    """边型坐标上的项 ((n,m),S,c) → 扇区坐标符号；sg=-1 给出共轭符号（X→1/X, Y→1/Y）。"""
    v = [{} for _ in range(15)]
    for (n, m), S, c in terms:
        i = MASKS.index(S)
        v[i] = add(v[i], {(sg*n, sg*m): Fr(c)})
    out = []
    for j in range(8):
        z = {}
        for i in range(15):
            if BM[i][j]:
                z = add(z, v[i])
        out.append(z)
    return out


def theta(terms, pi):
    out = []
    for (n, m), S, c in terms:
        dn, dm = step(S)
        out.append(((-n-dn, m+dm), perm_mask(pi, S), c))
    return out


def lin(*pairs):
    return [(x, S, Fr(c)*k) for terms, k in pairs for x, S, c in terms]


# ------------------------------------------------------- 精确有理被积函数 ----
MC = flint.fmpz_mpoly_ctx.get(('X', 'Y'), 'lex')
XG, YG = MC.gens()
SH = 4                                    # 统一单项式平移，使指数非负


def bareiss(M, zero):
    M = [r[:] for r in M]; n = len(M); prev = zero+1; sign = 1
    for k in range(n-1):
        if M[k][k] == 0:
            for r in range(k+1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]; sign = -sign; break
            else:
                return zero
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k]-M[i][k]*M[k][j])/prev
        prev = M[k][k]
    return sign*M[n-1][n-1]


def todict(p):
    return {tuple(int(x) for x in k): int(v) for k, v in p.to_dict().items()}


class Exact:
    def __init__(self, K, Gd, W):
        self.K, self.Gd, self.W = K, Gd, W
        N0 = [[dict(K[i][j]) for j in range(8)] for i in range(8)]
        for i in range(8):
            for j in range(8):
                for l in range(2):
                    N0[i][j] = add(N0[i][j], mul(Gd[i][l], conj(Gd[j][l])))
        self.L = 1
        for r in N0:
            for e in r:
                for v in e.values():
                    self.L = lcm(self.L, v.denominator)
        self.M = [[self.poly(N0[i][j]) for j in range(8)] for i in range(8)]
        self.D = bareiss(self.M, 0*XG)                       # det(K+GG^#)，乘 (L X^SH Y^SH)^8
        q = {}
        for i in range(8):
            for j in range(8):
                q = add(q, mul(mul(conj(W[i]), K[i][j]), W[j]))
        self.q = q                                            # w^# K w

    def poly(self, dc):
        p = 0*XG
        for (a, b), c in dc.items():
            v = c*self.L
            assert v.denominator == 1
            p += int(v)*XG**(a+SH)*YG**(b+SH)
        return p

    def witness(self, f, thf):
        a, b = sect(thf), sect(f, -1)
        for l in range(2):                                    # 规范不变：a·G ≡ 0，G^#·b ≡ 0（Laurent 恒等式）
            z1, z2 = {}, {}
            for i in range(8):
                z1 = add(z1, mul(a[i], self.Gd[i][l])); z2 = add(z2, mul(conj(self.Gd[i][l]), b[i]))
            assert z1 == {} and z2 == {}
        Mb = [r[:]+[self.poly(b[i])] for i, r in enumerate(self.M)]+[[self.poly(a[j]) for j in range(8)]+[0*XG]]
        B9 = bareiss(Mb, 0*XG)                                # det[[N0,b],[a,0]] = -(a adj N0 b)
        g = self.D.gcd(B9)
        P, Dr = B9/g, self.D/g
        aw, wb = {}, {}
        for i in range(8):
            aw = add(aw, mul(a[i], self.W[i])); wb = add(wb, mul(conj(self.W[i]), b[i]))
        # a K^+ b = -B9/(D·L·X^SH Y^SH) = P'/Dr，P' = -P/(L X^SH Y^SH)
        Pl = {(x-SH, y-SH): Fr(-v, self.L) for (x, y), v in todict(P).items()}
        return dict(a=a, b=b, P=Pl, Dr={k: Fr(v) for k, v in todict(Dr).items()}, aw=aw, wb=wb)


def check_Q(Dr):
    """约化分母 Dr 等于常数乘单项式乘 Q（规范不变行列式 c=-3Q/1280000 的三角多项式）。"""
    assert len(Dr) == len(QCOEF)
    sx = max(x for x, _ in Dr)-max(p for p, _ in QCOEF)
    sy = max(y for (x, y) in Dr if x == max(x2 for x2, _ in Dr))-max(q for (p, q) in QCOEF if p == max(p2 for p2, _ in QCOEF))
    cs = {Fr(Dr.get((p+sx, q+sy), 0))/a for (p, q), a in QCOEF.items()}
    assert len(cs) == 1
    return cs.pop(), (sx, sy)


# -------------------------------------------------------------- 迹旋转数值层 ----
def lp(dc, Yv):
    lo = min(a for a, b in dc); out = {}
    for (a, b), c in dc.items():
        yy = Yv**b if b >= 0 else (1/Yv)**(-b)
        out[a-lo] = out.get(a-lo, acb(0))+acb(fmpq(Fr(c).numerator, Fr(c).denominator))*yy
    return lo, acb_poly([out.get(i, acb(0)) for i in range(max(out)+1)])


def cen(p, x):
    m = acb(x.real.mid(), x.imag.mid())
    return p(acb_poly([m, 1]))(x-m)


def lev(lpv, X, centered=False):
    lo, p = lpv
    v = cen(p, X) if centered else p(X)
    return v*(X**lo if lo >= 0 else (1/X)**(-lo))


class Trace:
    """固定 s（Y 为点）的迹旋转被积函数 I = P/Dr - 2 (a·w)(w^#·b)/q 及其基础部分 P/Dr = a K^+ b。"""
    def __init__(self, ex, wit, Yv):
        self.P, self.Dr = lp(wit["P"], Yv), lp(wit["Dr"], Yv)
        self.q, self.aw, self.wb = lp(ex.q, Yv), lp(wit["aw"], Yv), lp(wit["wb"], Yv)

    def parts(self, X, centered=False):
        den, qq = lev(self.Dr, X, centered), lev(self.q, X, centered)
        if den.contains(0) or qq.contains(0):
            return None
        base = lev(self.P, X, centered)/den
        return base, base-2*lev(self.aw, X, centered)*lev(self.wb, X, centered)/qq

    def integral(self, which):
        k = {"base": 0, "trace": 1}[which]

        def f(t, _):
            r = self.parts((acb(0, 1)*t).exp())
            return acb("nan") if r is None else r[k]
        return acb.integral(f, 0, 2*arb.pi())/(2*arb.pi())

    def EN(self, N, nodes=None):
        tot = acb(0)
        for j in range(N):
            r = self.parts(acb.exp_pi_i(acb(fmpq(2*j, N))))
            assert r is not None
            tot += r[1]
        return tot/N

    def pole_gap(self):
        """Dr、q 作为 X 的多项式的全部根（Arb 隔离），返回 min |log|X|| 的下界。"""
        m = None
        for lo, p in (self.Dr, self.q):
            for r in p.roots():
                lg = abs(abs(r).log())
                m = lg if m is None else m.min(lg)
        return arb(m.lower())

    def sup_circle(self, rho, ncell=1500):
        Mb, Mt = arb(0), arb(0)
        stack = [(fmpq(2*j, ncell), fmpq(2*j+2, ncell)) for j in range(ncell)]
        while stack:
            lo, hi = stack.pop()
            ang = arb.pi()*(arb(lo)+arb(hi))/2+arb(0, (arb.pi()*(arb(hi)-arb(lo))/2).upper())
            r = self.parts(acb(rho)*(acb(0, 1)*acb(ang)).exp(), centered=True)
            if r is None or any(abs(x).rad() > abs(x).mid()*2+1 for x in r):
                assert hi-lo > fmpq(1, 10**7)
                mid = (lo+hi)/2; stack += [(lo, mid), (mid, hi)]
                continue
            Mb = Mb.max(abs(r[0]).upper()); Mt = Mt.max(abs(r[1]).upper())
        return Mb, Mt


def strip_data(T):
    """a < 极点距离；在 |log|X||=±a 两个圆上取上确界（最大模原理给出整个环上的界）。"""
    gap = T.pole_gap()
    a = arb(fmpq(3, 4))*gap
    a = arb(a.lower())
    b1, t1 = T.sup_circle(a.exp()); b2, t2 = T.sup_circle((-a).exp())
    return gap, a, b1.max(b2), t1.max(t2)


def tail(M, a, N):
    """梯形误差 |E_N - E_∞| ≤ 2M/(e^{aN}-1)（被积函数在 |Im t|≤a 的闭带上解析、|I|≤M）。"""
    return 2*M/((a*N).exp()-1)


# -------------------------------------------------------------- 谱旋转数值层 ----
SQ = [arb(x).sqrt() for x in DD]


def pw(X, k):
    t = {0: acb(1)}; Xi = X.conjugate()               # |X|=1
    for p in range(1, k+1):
        t[p] = t[p-1]*X; t[-p] = t[-p+1]*Xi
    return t


def lacb(dc, PX, PY):
    z = acb(0)
    for (p, q), c in dc.items():
        z += acb(fmpq(Fr(c).numerator, Fr(c).denominator))*PX[p]*PY[q]
    return z


class Spectral:
    """单位圆节点上的谱旋转附加项 Δ = -2 (a D^{-1/2} v)(v† D^{-1/2} b)/(λ v†v)，A v = λ v，A = D^{-1/2} K D^{-1/2}。"""
    def __init__(self, ex, wits, sq_):
        self.K = ex.K; self.Y = acb.exp_pi_i(acb(sq_)); self.PY = pw(self.Y, 40)
        self.wits = wits

    def node(self, X):
        PX = pw(X, 40)
        A = acb_mat([[lacb(self.K[i][j], PX, self.PY)/(SQ[i]*SQ[j]) for j in range(8)] for i in range(8)])
        An = np.array([[complex(float(A[i, j].real.mid()), float(A[i, j].imag.mid())) for j in range(8)] for i in range(8)])
        e, U = np.linalg.eigh(An)
        assert (e < -1e-9).sum() == 1
        I8 = acb_mat([[int(i == j) for j in range(8)] for i in range(8)])
        lo, hi = arb(e[0])-arb(1e-12), arb(e[0])+arb(1e-12)
        # μ<0 时 det(A-μ) 的符号等于 λ-μ 的符号（惯性 (5,1,2)：其余本征值 ≥ 0）
        assert hi < 0 and (A-acb(lo)*I8).det().real > 0 and (A-acb(hi)*I8).det().real < 0
        lam = arb(lo.mid()).union(arb(hi.mid()))
        Ml = A-acb(lam)*I8; col = int(np.argmax(np.abs(U[:, 0])))
        rows = [r for r in range(8) if r != col]
        Mr = acb_mat([[Ml[r, c] for c in rows] for r in rows])
        rhs = acb_mat([[-Ml[r, col]] for r in rows])
        xs = Mr.solve(rhs)
        v = [acb(1) if i == col else xs[rows.index(i), 0] for i in range(8)]
        vv = sum((x*x.conjugate()).real for x in v)
        out = []
        for w in self.wits:
            av = sum((lacb(w["a"][i], PX, self.PY)*v[i]/SQ[i] for i in range(8) if w["a"][i]), acb(0))
            vb = sum((v[i].conjugate()*lacb(w["b"][i], PX, self.PY)/SQ[i] for i in range(8) if w["b"][i]), acb(0))
            out.append(-2*av*vb/(acb(lam)*vv))
        return lam, out


def lip_A(K):
    """A 系数优势矩阵：W1_ij = Σ|c||p|（Weyl 用），W(η)_ij = Σ|c|(e^{|p|η}-1)（Bauer–Fike 用）。"""
    rows = []
    for i in range(8):
        r = []
        for j in range(8):
            byp = {}
            for (p, q), c in K[i][j].items():
                byp[p] = byp.get(p, 0)+abs(c)
            r.append({p: v/Fr(1) for p, v in byp.items()})
        rows.append(r)
    return rows


def frob(rows, fn):
    tot = arb(0)
    for i in range(8):
        for j in range(8):
            e = sum((arb(fmpq(v.numerator, v.denominator))*fn(p) for p, v in rows[i][j].items()), arb(0))/(SQ[i]*SQ[j])
            tot += e*e
    return tot.sqrt()


def vec_bound(dcs, eta):
    """strip |Im t|≤η 上 ‖D^{-1/2} a(X)‖ 的上界（|Y|=1）。"""
    tot = arb(0)
    for i, dc in enumerate(dcs):
        s = sum((arb(fmpq(abs(Fr(c)).numerator, abs(Fr(c)).denominator))*(arb(abs(p))*eta).exp()
                 for (p, q), c in dc.items()), arb(0))
        tot += s*s/arb(DD[i])
    return tot.sqrt()


def spectral_block(ex, wits, sq_, eta, Ninf, small_N, Mbase, abase, base_inf):
    """返回各见证的 E_∞(谱) 区间、小 N 的 E_N(谱)、尾部常数。"""
    S = Spectral(ex, wits, sq_)
    Tb = [Trace(ex, w, S.Y) for w in wits]
    sums = [acb(0) for _ in wits]; lmin = None
    for j in range(Ninf):
        lam, o = S.node(acb.exp_pi_i(acb(fmpq(2*j, Ninf))))
        lmin = abs(lam) if lmin is None else lmin.min(abs(lam))
        for k, x in enumerate(o):
            sums[k] += x
    rows = lip_A(ex.K)
    L1 = frob(rows, lambda p: arb(abs(p)))
    lam_lo = arb((lmin-L1*arb.pi()/Ninf).lower())          # 全圆上 |λ(t)| 的下界（Weyl）
    assert lam_lo > 0
    delta = frob(rows, lambda p: (arb(abs(p))*eta).exp()-1)
    assert delta.upper() < (lam_lo/4).lower()               # δ ≤ ρ/2，ρ = λ_min/2
    res = []
    for k, w in enumerate(wits):
        Ma, Mb_ = vec_bound(w["a"], eta), vec_bound(w["b"], eta)
        MD = 2*Ma*Mb_*2/(arb(3)/4*lam_lo)                    # ‖P‖≤2，|λ(X)|≥3λ_min/4
        err = tail(MD, eta, Ninf)
        Delta = sums[k]/Ninf
        Einf = base_inf[k].real+Delta.real+arb(0, err.upper())
        res.append(dict(Einf=Einf, MD=MD, err=err))
    EN = {}
    for N in small_N:
        tot = [acb(0) for _ in wits]
        for j in range(N):
            X = acb.exp_pi_i(acb(fmpq(2*j, N)))
            lam, o = S.node(X)
            for k in range(len(wits)):
                tot[k] += Tb[k].parts(X)[0]+o[k]
        EN[N] = [x/N for x in tot]
    return res, EN, dict(lam_lo=lam_lo, delta=delta, L1=L1)


# ------------------------------------------------------------------ 主流程 ----
def verify():
    t0 = time.time()
    K, Gd, W, _, _ = PS["symbols"]()
    ex = Exact(K, Gd, W)
    F1 = [(x, S, Fr(int(c.p), int(c.q))) for x, S, c in PS["F1"]]
    G = [(x, S, Fr(int(c.p), int(c.q))) for x, S, c in PS["G_PI2"]]
    H = lin((G, 1), (F1, -16))
    WIT = {"f1": F1, "g": G, "h": H}
    pi0 = SITE_REFL[0]
    span = lambda terms: (min(n for (n, m), S, c in terms), max(n+step(S)[0] for (n, m), S, c in terms))
    supp = {"f1": (0, 2), "g": (0, 4), "h": (0, 4)}
    for name, f in WIT.items():                       # 见证与其反射像的支集（按边的两端点）
        lo, hi = span(f); assert lo >= supp[name][0] and hi <= supp[name][1]
        for pi in SITE_REFL:
            lo, hi = span(theta(f, pi)); assert lo >= -4 and hi <= 0
    for g_ in [(1, 0, 2, 3), (0, 1, 3, 2)]:           # 三个见证在 (12)、(34) 下不变
        for f in WIT.values():
            dct = {}
            for x, S, c in f:
                dct[(x, S)] = dct.get((x, S), 0)+c
            dct = {k: v for k, v in dct.items() if v}
            assert {(x, perm_mask(g_, S)): v for (x, S), v in dct.items()} == dct
    wd = {}
    for name, f in WIT.items():
        for pi in SITE_REFL:
            wd[(name, pi)] = ex.witness(f, theta(f, pi))
    cQ, shQ = check_Q(wd[("f1", pi0)]["Dr"])
    assert cQ == -3 and shQ == (6, 3)
    for k, w in wd.items():
        assert check_Q(w["Dr"]) == (cQ, shQ)
    t1 = time.time()
    print(f"[exact] a K^+ b = P/Dr，Dr = {cQ}·X^{shQ[0]}Y^{shQ[1]}·Q；规范不变恒等式成立（{t1-t0:.1f} 秒）")

    # ---------- 迹旋转：四个盒块，n 方向无限体积与一切偶 N≥8 ----
    blocks = {  # 名称: (s/π, 见证, 块内 s 值/π, 块在盒中所占动量比例, 已有 N=8 精确值)
        "pi/4": (fmpq(1, 4), "f1", [fmpq(1, 4), fmpq(3, 4)], fmpq(1), fmpq(-270721840955, 3088612034704)),
        "pi/6": (fmpq(1, 6), "f1", [fmpq(1, 6), fmpq(5, 6)], fmpq(2, 3), fmpq(-30700157525965, 155183768488104)),
        "pi/3": (fmpq(1, 3), "f1", [fmpq(1, 3), fmpq(2, 3)], fmpq(2, 3), fmpq(-182051589335, 15708820029048)),
        "pi/2": (fmpq(1, 2), "g", [fmpq(1, 2)], fmpq(1, 3), fmpq(-8016671, 7114338)),
    }
    expect_tr = {"pi/4": (-0.0567, -0.0567), "pi/6": (-0.1657, -0.1657), "pi/3": (-0.00878, -0.00878), "pi/2": (-3.997, -3.997)}
    TR = {}
    for name, (sq_, wn, svals, frac, box8) in blocks.items():
        Yv = acb.exp_pi_i(acb(sq_))
        Ts = {pi: Trace(ex, wd[(wn, pi)], Yv) for pi in SITE_REFL}
        Ei = {pi: Ts[pi].integral("trace") for pi in SITE_REFL}
        assert Ei[SITE_REFL[0]].real.overlaps(Ei[SITE_REFL[1]].real)
        for pi in SITE_REFL:
            assert abs(Ei[pi].imag) < arb(1e-20) and Ei[pi].real.rad() < arb(1e-15)
        for s2 in svals[1:]:                                   # E(π-s)=E(s)
            E2 = Trace(ex, wd[(wn, pi0)], acb.exp_pi_i(acb(s2))).integral("trace")
            assert E2.real.overlaps(Ei[pi0].real)
        T = Ts[pi0]
        gap, a, Mb, Mt = strip_data(T)
        E8 = T.EN(8)
        assert (E8.real*arb(frac)).overlaps(arb(box8.p)/arb(box8.q))          # 复现已有 8 盒精确值
        Einf = Ei[pi0].real
        assert Einf < 0
        N1 = 8
        while not (tail(Mt, a, N1) < -Einf/1000):
            N1 += 2
        EN = {N: T.EN(N).real for N in range(8, N1, 2)}
        U = Einf+tail(Mt, a, N1)
        for v in EN.values():
            U = U.max(v)
        assert U < 0
        TR[name] = dict(EN=EN, Einf=Einf, gap=gap, a=a, Mb=Mb, Mt=Mt, N1=N1, U=U, E8=E8.real, T=T,
                        base=[Trace(ex, wd[(wn, pi0)], Yv).integral("base")])
    # 严格半空间见证（f1 平移一层，支撑在 1≤n≤3）在 s=π/4 的无限体积值
    F1s = [((x[0]+1, x[1]+1), S, c) for x, S, c in F1]
    lo, hi = span(F1s); assert lo >= 1 and hi <= 3
    ws = ex.witness(F1s, theta(F1s, pi0))
    Estrict = Trace(ex, ws, acb.exp_pi_i(acb(fmpq(1, 4)))).integral("trace").real
    assert Estrict < 0 and check_Q(ws["Dr"]) == (cQ, shQ)
    t2 = time.time()
    print("[迹旋转, Arb] n 方向无限体积 E_∞ 与一切偶 N≥8 的一致上界 U：")
    for name, r in TR.items():
        print(f"  s={name}: E_∞={r['Einf'].str(10, radius=True)}  E_8={r['E8'].str(8)}  极点距离>{r['gap'].str(5)}  "
              f"a={r['a'].str(4)}  M≤{r['Mt'].str(4)}  N1={r['N1']}  sup_N E_N ≤ {r['U'].str(6)}")
    print(f"  严格半空间见证 s=π/4：E_∞={Estrict.str(10, radius=True)}（{t2-t1:.1f} 秒）")
    assert TR["pi/4"]["Einf"] < arb(-0.0567) and TR["pi/6"]["Einf"] < arb(-0.165)
    assert TR["pi/3"]["Einf"] < arb(-0.0087) and TR["pi/2"]["Einf"] < arb(-3.99)
    for name, ub in {"pi/4": -0.0566, "pi/6": -0.1655, "pi/3": -0.00878, "pi/2": -3.38}.items():
        assert TR[name]["U"] < arb(ub)
    assert [TR[k]["N1"] for k in ["pi/4", "pi/6", "pi/3", "pi/2"]] == [48, 68, 40, 26]
    for k in ["pi/4", "pi/6", "pi/3"]:                     # 所算范围内 E_N 随 N 单调上升趋于 E_∞
        Ns = sorted(TR[k]["EN"])
        assert all(TR[k]["EN"][a] < TR[k]["EN"][b] for a, b in zip(Ns, Ns[1:]))
    assert all(TR["pi/2"]["EN"][N] < TR["pi/2"]["EN"][8] for N in TR["pi/2"]["EN"] if N > 8)
    r2 = TR["pi/2"]                                        # N≥N1 时 E_N ≤ E+误差 < E_8，故 sup_N E_N = E_8
    assert r2["Einf"]+tail(r2["Mt"], r2["a"], r2["N1"]) < r2["E8"]
    near = lambda x, v, r: abs(x-arb(v)) < arb(r)
    doc = {"pi/4": ("-0.05674521474", 3e-12, 0.3438, "0.2579", 5.583), "pi/6": ("-0.1657630017", 3e-11, 0.2321, "0.1741", 11.38),
           "pi/3": ("-0.008789596127", 7e-13, 0.4489, "0.3367", 2.832), "pi/2": ("-3.997512698", 7e-10, 0.6263, "0.4698", 217.4)}
    for k, (v, r, gp, al, Mx) in doc.items():
        assert near(TR[k]["Einf"], v, r) and TR[k]["gap"] > arb(gp) and near(TR[k]["a"], al, 1e-4) and TR[k]["Mt"] < arb(Mx)
    assert near(TR["pi/4"]["E8"], "-0.087651618", 1e-9) and near(TR["pi/2"]["E8"], "-3.3804991", 8e-8)
    assert near(Estrict, "-0.04704275856", 8e-12)

    # ---------- 迹旋转：m 方向长度 M≤12 的全部非零横向动量类 ----
    SLIST = [fmpq(1, 12), fmpq(1, 10), fmpq(1, 8), fmpq(1, 6), fmpq(1, 5), fmpq(1, 4), fmpq(3, 10), fmpq(1, 3),
             fmpq(3, 8), fmpq(2, 5), fmpq(5, 12), fmpq(1, 2)]
    EXT = {}
    for sq_ in SLIST:
        wn = "h" if sq_ <= fmpq(1, 3) else "g"
        T = Trace(ex, wd[(wn, pi0)], acb.exp_pi_i(acb(sq_)))
        v = T.integral("trace").real
        assert v < 0 and v.rad() < arb(1e-8)
        gap, a, _, Mt = strip_data(T)                        # 与定理 2 同一流程：一切偶 N≥8 的一致上界
        N1 = 8
        while not (tail(Mt, a, N1) < -v/1000):
            N1 += 2
        EN = {N: T.EN(N).real for N in range(8, N1, 2)}
        U = v+tail(Mt, a, N1)
        for x in EN.values():
            U = U.max(x)
        assert U < 0
        EXT[sq_] = (wn, v, U, N1, EN)
    EXTF1 = {sq_: Trace(ex, wd[("f1", pi0)], acb.exp_pi_i(acb(sq_))).integral("trace").real for sq_ in SLIST}
    assert all((v < 0) if sq_ <= fmpq(1, 3) else (v > 0) for sq_, v in EXTF1.items())
    ubx = {fmpq(1, 12): -471.6, fmpq(1, 10): -366.4, fmpq(1, 8): -262.1, fmpq(1, 6): -159.7, fmpq(1, 5): -109.9, fmpq(1, 4): -62.09,
           fmpq(3, 10): -32.38, fmpq(1, 3): -18.76, fmpq(3, 8): -1.618, fmpq(2, 5): -2.515, fmpq(5, 12): -2.982, fmpq(1, 2): -3.997}
    assert all(EXT[k][1] < arb(v) for k, v in ubx.items())
    ubU = {fmpq(1, 12): (-471.2, 138), fmpq(1, 10): (-366.1, 114), fmpq(1, 8): (-261.9, 90), fmpq(1, 6): (-159.5, 68),
           fmpq(1, 5): (-109.8, 58), fmpq(1, 4): (-62.03, 46), fmpq(3, 10): (-32.35, 40), fmpq(1, 3): (-18.75, 38),
           fmpq(3, 8): (-0.2967, 36), fmpq(2, 5): (-1.487, 32), fmpq(5, 12): (-2.094, 32), fmpq(1, 2): (-3.38, 26)}
    for k, (ub, n1) in ubU.items():
        assert EXT[k][2] < arb(ub) and EXT[k][3] == n1
    print("  一切偶 N≥8 的一致上界 (N1, U)：", {str(k): (e[3], e[2].str(6)) for k, e in EXT.items()})
    assert abs(EXTF1[fmpq(1, 2)]-arb("0.0195751")) < arb(7e-8)
    t3 = time.time()
    print("[迹旋转, Arb] m 方向长度 ≤12 的横向动量类 s/π 上的 E_∞（见证 h=g-16f1 或 g）：")
    print("  ", {str(k): (e[0], e[1].str(6)) for k, e in EXT.items()})
    print("  f1 单独的 E_∞（s/π≤1/3 为负，≥3/8 为正）：", {str(k): v.str(6) for k, v in EXTF1.items()}, f"（{t3-t2:.1f} 秒）")

    # ---------- 谱旋转：四个块 ----
    SPEC = {
        "pi/4": (fmpq(1, 4), ["f1"], fmpq(1, 250), 3200),
        "pi/6": (fmpq(1, 6), ["f1"], fmpq(1, 560), 6800),
        "pi/3": (fmpq(1, 3), ["f1", "g"], fmpq(1, 144), 2400),
        "pi/2": (fmpq(1, 2), ["g"], fmpq(1, 72), 900),
    }
    small_N = list(range(8, 66, 2))
    SR = {}
    for name, (sq_, wns, eta, Ninf) in SPEC.items():
        wits = [wd[(wn, pi0)] for wn in wns]
        Yv = acb.exp_pi_i(acb(sq_))
        base_inf = [Trace(ex, w, Yv).integral("base") for w in wits]
        res, EN, info = spectral_block(ex, wits, sq_, arb(eta), Ninf, small_N, None, None, base_inf)
        for k, w in enumerate(wits):                       # 大 N 尾部：|E_N-E_∞| ≤ 2M_base/(e^{aN}-1)+2M_Δ/(e^{ηN}-1)
            _, a_, Mb_, _ = strip_data(Trace(ex, w, Yv))
            E = res[k]["Einf"]; N1 = 66
            while not (tail(Mb_, a_, N1)+tail(res[k]["MD"], arb(eta), N1) < abs(E).lower()/2):
                N1 += 2
            res[k]["N1"] = N1
        SR[name] = (wns, res, EN, info, eta, Ninf)
    t4 = time.time()
    print("[谱旋转, Arb] E_∞（单位圆节点和 + 带状解析误差界）与 8≤N≤64 的 E_N：")
    for name, (wns, res, EN, info, eta, Ninf) in SR.items():
        for k, wn in enumerate(wns):
            print(f"  s={name} 见证 {wn}: E_∞ ∈ {res[k]['Einf'].str(8, radius=True)}  (N={Ninf}, η={eta}, 误差≤{res[k]['err'].str(3)})"
                  f"  E_8={EN[8][k].real.str(6)}  E_64={EN[64][k].real.str(6)}  N1={res[k]['N1']}")
        print(f"    |λ|≥{info['lam_lo'].str(4)}  δ(η)≤{info['delta'].str(4)}")
    f1i = SR["pi/3"][0].index("f1"); gi = SR["pi/3"][0].index("g")
    E_f1_pi3 = SR["pi/3"][2]
    assert SR["pi/4"][1][0]["Einf"] < arb(-0.053) and SR["pi/6"][1][0]["Einf"] < arb(-0.164)
    assert SR["pi/3"][1][gi]["Einf"] < arb(-1.613) and SR["pi/2"][1][0]["Einf"] < arb(-6.67)
    assert [SR[nm][1][0]["N1"] for nm in ["pi/4", "pi/6", "pi/2"]] == [2368, 5118, 474]
    assert SR["pi/3"][1][f1i]["N1"] == 1742 and SR["pi/3"][1][gi]["N1"] == 1246
    assert SR["pi/3"][1][f1i]["Einf"] > arb(0.0022) and SR["pi/3"][1][f1i]["Einf"] < arb(0.0024)   # f1 不再是无限体积见证
    tab = {"pi/4": (0.0692, 0.01606, [0.000965], [("-0.137890", 8e-7)]), "pi/6": (0.0314, 0.007148, [0.00409], [("-0.430456", 7e-7)]),
           "pi/3": (0.1203, 0.02795, [1.18e-5, 0.000267], [("-0.0275845", 1e-7), ("-0.0013995", 1e-7)]),
           "pi/2": (0.2367, 0.05629, [0.00893], [("-5.70542", 1e-5)])}
    for nm, (lw, dl, errs, e8) in tab.items():
        wns, res, EN, info, _, _ = SR[nm]
        assert info["lam_lo"] > arb(lw) and info["delta"] < arb(dl)
        for k in range(len(wns)):
            assert res[k]["err"] < arb(errs[k]) and abs(EN[8][k].real-arb(e8[k][0])) < arb(e8[k][1])
    # 与已有 8 盒 Arb 值（块因子 κ）一致：<-0.137、<-0.286、<-0.0183、<-1.90
    assert SR["pi/4"][2][8][0].real < arb(-0.137) and arb(fmpq(2, 3))*SR["pi/6"][2][8][0].real < arb(-0.286)
    assert arb(fmpq(2, 3))*SR["pi/3"][2][8][f1i].real < arb(-0.0183) and arb(fmpq(1, 3))*SR["pi/2"][2][8][0].real < arb(-1.90)
    for N in small_N:
        v = E_f1_pi3[N][f1i].real
        assert (v < 0) if N <= 12 else (v > 0), (N, v)
        assert E_f1_pi3[N][gi].real < 0
        for nm in ["pi/4", "pi/6", "pi/2"]:
            assert SR[nm][2][N][0].real < 0
    assert (SR["pi/4"][2][8][0].real).overlaps(arb(-0.13789)+arb(0, 3e-5))
    print("  s=π/3 的 f1：N≤12 时 E_N<0，14≤N≤64 时 E_N>0；g 在 8≤N≤64 全为负")
    print(f"  用时 {t4-t3:.1f} 秒；总用时 {time.time()-t0:.1f} 秒")


if __name__ == "__main__":
    verify()

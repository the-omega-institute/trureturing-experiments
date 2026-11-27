"""三个归一无关读数比在两点格以外的局部可识别性：十九点壳只保留最短弦的子图（人为子图）的精确约化，
三角形图与三点开链在静态极限附近的有效重叠结构与雅可比行列式的首个非零阶 t^16。
配套说明：三比值在三角形图与三点链上于参数取三分之一与二分之一处到t16阶才局部确定重叠模长.md

探针与 ch39-proton-current-overlap-hopping.py 相同（人为定义的读数）：每格点一代 Q,L 及其共轭，单点接触项 kappa H_x，最近邻标量跳跃 t，
各格点参数相同。读数在格点 0：X_AB, X_AA, X_BB，作为 (u,v,rho,t) 的函数，u=kappa K11, v=kappa K22。

方法（全部精确有理运算）：
(1) 原始 eps(x)eps 基、物理 n_f=4。单点 exp(kappa H) 按缺失生成元集合 (M_u, M_b) 展开为振幅
    A(M_u,M_b)（|M_u|,|M_b|<=8，单侧元素由 D^d (QL)^l 的直接格拉斯曼乘积得到）。跳跃把点 x 的共轭生成元
    与邻点同分量的非共轭生成元配对，Z 的 t^n 系数是所有配置上振幅之积乘格拉斯曼符号之和。
(2) 校验：两点格 t^0, t^4, t^8 系数（含点 0 的一、二阶导数）与 ch39-proton-current-overlap-hopping.py 的 225 项多项式逐项相等。
(3) 三角形 t^6 项、三点开链（中心失去 8 个生成元）的 t^8 联合项；log Z 的团簇展开给出读数到 t^8。
(4) 齐次性 R_{2k}(lam u, lam v) = lam^{-k} R_{2k}(u,v) 在 u=2 处核对；u 方向导数由欧拉关系取得。
(5) 雅可比行列式按 t 展开到 t^16：按多重线性逐项求和。未算出的列（rho 方向的 t^6 及以上、u,v 方向的
    t^10, t^12）在 t^16 以内只出现在另两列平行的项中，脚本逐项断言这一平行性；断言 t^8..t^14 系数恒为零，
    t^16 系数等于给定有理函数。
(6) 两点格的完整有理函数在静态射线 U=V=L 上的精确值与 t^16 系数一致。
"""
import itertools
import runpy
import time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from math import factorial as fa
from pathlib import Path

import flint
import sympy as sp

HERE = Path(__file__).resolve().parent
_static = runpy.run_path(str(HERE / "ch39-proton-current-overlap-static.py"))
_hop = runpy.run_path(str(HERE / "ch39-proton-current-overlap-hopping.py"))
mul = _static["mul"]
pc = lambda m: bin(m).count("1")
FULL = (1 << 16)-1
DERS = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
KEYS = ["AB", "AA", "BB"]


# ---------------------------------------------------------------- 十九点壳只保留最短弦的子图（人为子图）与模型壳图 G19
def shell_graph():
    P = [(a, b, c) for a in range(30) for b in range(30) for c in range(30) if a*a+b*b+c*c == 867]
    d2 = lambda p, q: sum((x-y)**2 for x, y in zip(p, q))
    m = min(d2(p, q) for p, q in itertools.combinations(P, 2))
    E = [(p, q) for p, q in itertools.combinations(P, 2) if d2(p, q) == m]
    deg = Counter(x for e in E for x in e)
    assert len(P) == 19 and m == 12 and len(E) == 3
    assert set(deg.values()) == {1} and len(deg) == 6  # 三条互不相交的边，十三个孤立点
    # 模型的壳图 G19 取 d^2 <= 804：135 条边、529 个三角形、度分布 10x6, 14x3, 16x6, 18x4。
    # 上面只取其中 d^2 = 12 的 3 条边，是人为子图。
    G = {frozenset((p, q)) for p, q in itertools.combinations(P, 2) if d2(p, q) <= 804}
    gdeg = Counter(x for e in G for x in e)
    ntri = sum(1 for a, b, c in itertools.combinations(P, 3)
               if frozenset((a, b)) in G and frozenset((b, c)) in G and frozenset((a, c)) in G)
    assert len(G) == 135 and ntri == 529 and Counter(gdeg.values()) == {10: 6, 14: 3, 16: 6, 18: 4}
    assert all(frozenset(e) in G for e in E)
    adj = {x: {y for y in P if frozenset((x, y)) in G} for x in P}
    a0, c0 = next((a, c) for a, c in itertools.combinations(P, 2) if len(adj[a] & adj[c]) >= 2)
    assert len(adj[a0] & adj[c0]) >= 2  # G19 含四边形 a0-b-c0-b'-a0
    return P, E


# ---------------------------------------------------------------- 单点振幅
def site_amplitudes():
    D, QL = _static["bilinears"](False)
    Db, T = _static["bilinears"](True)

    def elements(F, S):
        pf, ps = [], []
        for b in F:
            o = [{0: 1}]
            for _ in range(4):
                o.append(mul(o[-1], b))
            pf.append(o)
        for b in S:
            o = [{0: 1}]
            for _ in range(4):
                o.append(mul(o[-1], b))
            ps.append(o)
        out = {}

        def rec(a, acc, key):
            if a == 3:
                if sum(key) >= 4:
                    out[key] = acc
                return
            for d in range(5):
                for l in range(5):
                    if sum(key)+d+l > 8:
                        continue
                    r = mul(mul(acc, pf[a][d]), ps[a][l])
                    if r:
                        rec(a+1, r, key+(d, l))
        rec(0, {0: 1}, ())
        return out

    def bymiss(E):
        out = defaultdict(list)
        for k, poly in E.items():
            for m, c in poly.items():
                out[FULL ^ m].append((k, c))
        return out
    mu_, mb_ = bymiss(elements(D, QL)), bymiss(elements(Db, T))
    wcache = {}

    def Wdict(a, b):
        if (a, b) in wcache:
            return wcache[a, b]
        res = {(0, 0, 0, 0): Fr(1)}
        for c in range(3):
            d, l, db, lb = a[2*c], a[2*c+1], b[2*c], b[2*c+1]
            part = {}
            for n11 in range(5):
                n21, n12, n22 = d-n11, db-n11, l-d+n11
                if min(n21, n12, n22) < 0 or n22+n12 != lb:
                    continue
                part[(n11, n22, n21, n12)] = Fr(1, fa(n11)*fa(n22)*fa(n21)*fa(n12))
            if not part:
                res = {}
                break
            new = defaultdict(Fr)
            for e, cv in res.items():
                for e2, c2 in part.items():
                    new[tuple(x+y for x, y in zip(e, e2))] += cv*c2
            res = new
        wcache[a, b] = dict(res)
        return wcache[a, b]
    bysize_u, bysize_b = defaultdict(list), defaultdict(list)
    for Mu in mu_:
        bysize_u[pc(Mu)].append(Mu)
    for Mb in mb_:
        bysize_b[pc(Mb)].append(Mb)
    AMP = {}
    for su in bysize_u:
        for sb in bysize_b:
            if su+sb not in (0, 4, 8):
                continue
            for Mu in bysize_u[su]:
                for Mb in bysize_b[sb]:
                    acc = defaultdict(Fr)
                    for ka, ca in mu_[Mu]:
                        for kb, cb in mb_[Mb]:
                            for e, cv in Wdict(ka, kb).items():
                                acc[e] += cv*ca*cb
                    acc = {e: v for e, v in acc.items() if v}
                    if acc:
                        AMP[Mu, Mb] = acc
    cnt = Counter(pc(a)+pc(b) for a, b in AMP)
    assert cnt == {0: 1, 4: 168, 8: 6592}
    return AMP


def ff(n, i):
    r = 1
    for k in range(i):
        r *= n-k
    return r


def amp_poly(acc, rho, der=(0, 0), uv=1):
    """u=uv, v=uv*w, x=rho*uv^2*w, y=1；der=(i,j) 为对 u、v 在固定 x,y 下的导数。"""
    i, j = der
    co = defaultdict(Fr)
    for (n11, n22, n21, n12), c in acc.items():
        f = ff(n11, i)*ff(n22, j)
        if f:
            co[n22-j+n21] += c*f*rho**n21*Fr(uv)**(n11-i+n22-j+2*n21)
    if not co:
        return flint.fmpq_poly([])
    return flint.fmpq_poly([flint.fmpq(co[k].numerator, co[k].denominator) for k in range(max(co)+1)])


def sign_of(sites, hops):
    """全局次序 [u0][b0][u1][b1]...；单点单项按此次序已排好，依次乘跳跃双线性 bbar_{x,k} psi_{y,k}。"""
    M = 0
    for x, (pu, pb) in enumerate(sites):
        M |= (pu << (32*x)) | (pb << (32*x+16))
    s = 1
    for (bx, k, uy) in hops:
        for g in (32*bx+16+k, 32*uy+k):
            if M >> g & 1:
                return 0
            if pc(M >> (g+1)) & 1:
                s = -s
            M |= 1 << g
    return s


comp = lambda key: (FULL ^ key[0], FULL ^ key[1])
bits = lambda m: [k for k in range(16) if m >> k & 1]


def zero_polys():
    return {d: flint.fmpq_poly([]) for d in DERS}


def two_point(AMP, rho, loss, uv):
    out = zero_polys()
    for (Mu0, Mb0), acc0 in AMP.items():
        if pc(Mu0)+pc(Mb0) != loss or (Mb0, Mu0) not in AMP:
            continue
        hops = [(0, k, 1) for k in bits(Mb0)]+[(1, k, 0) for k in bits(Mu0)]
        s = sign_of([comp((Mu0, Mb0)), comp((Mb0, Mu0))], hops)
        a1 = amp_poly(AMP[Mb0, Mu0], rho, uv=uv)
        for d in DERS:
            out[d] += s*amp_poly(acc0, rho, d, uv)*a1
    return out


def triangle_t6(AMP, rho, uv):
    """三角形 t^6：每条边两次跳跃，每点失去 4 个生成元（其余分配由宇称排除）。"""
    hl = {}
    for (x, y) in [(0, 1), (0, 2), (1, 2)]:
        hs = [(x, k, y) for k in range(16)]+[(y, k, x) for k in range(16)]
        hl[x, y] = list(itertools.combinations(hs, 2))

    def miss(site, hops):
        Mu = Mb = 0
        for (bx, k, uy) in hops:
            if bx == site:
                if Mb >> k & 1:
                    return None
                Mb |= 1 << k
            if uy == site:
                if Mu >> k & 1:
                    return None
                Mu |= 1 << k
        return (Mu, Mb)
    A4 = {k: v for k, v in AMP.items() if pc(k[0])+pc(k[1]) == 4}
    A1 = {k: amp_poly(v, rho, uv=uv) for k, v in A4.items()}
    A0 = {k: {d: amp_poly(v, rho, d, uv) for d in DERS} for k, v in A4.items()}
    out, n = zero_polys(), 0
    for h01 in hl[0, 1]:
        for h02 in hl[0, 2]:
            m0 = miss(0, h01+h02)
            if m0 is None or m0 not in A0:
                continue
            for h12 in hl[1, 2]:
                m1 = miss(1, h01+h12)
                if m1 is None or m1 not in A1:
                    continue
                m2 = miss(2, h02+h12)
                if m2 is None or m2 not in A1:
                    continue
                s = sign_of([comp(m0), comp(m1), comp(m2)], h01+h02+h12)
                if s:
                    n += 1
                    pr = A1[m1]*A1[m2]
                    for d in DERS:
                        out[d] += s*A0[m0][d]*pr
    return out, n


def chain_joint(AMP, rho, uv):
    """三点开链 t^8 联合项：中心 0 在两条键上各四次跳跃（失去 8 个），端点各失去 4 个。
    返回（对中心求导，对端点 1 求导）。"""
    A4 = {k: v for k, v in AMP.items() if pc(k[0])+pc(k[1]) == 4}
    P4 = {k: amp_poly(v, rho, uv=uv) for k, v in A4.items()}
    D4 = {k: {d: amp_poly(v, rho, d, uv) for d in DERS} for k, v in A4.items()}
    outc, oute, n = zero_polys(), zero_polys(), 0
    for (Mu, Mb), acc in AMP.items():
        if pc(Mu)+pc(Mb) != 8:
            continue
        items = [("u", k) for k in bits(Mu)]+[("b", k) for k in bits(Mb)]
        cp = None
        for sub in itertools.combinations(range(8), 4):
            s1 = [items[i] for i in sub]
            s2 = [items[i] for i in range(8) if i not in sub]
            S1u, S1b = sum(1 << k for t, k in s1 if t == "u"), sum(1 << k for t, k in s1 if t == "b")
            S2u, S2b = sum(1 << k for t, k in s2 if t == "u"), sum(1 << k for t, k in s2 if t == "b")
            e1, e2 = (S1b, S1u), (S2b, S2u)
            if e1 not in A4 or e2 not in A4:
                continue
            hops = ([(0, k, 1) for k in bits(S1b)]+[(1, k, 0) for k in bits(S1u)]
                    + [(0, k, 2) for k in bits(S2b)]+[(2, k, 0) for k in bits(S2u)])
            s = sign_of([comp((Mu, Mb)), comp(e1), comp(e2)], hops)
            if not s:
                continue
            n += 1
            if cp is None:
                cp = {d: amp_poly(acc, rho, d, uv) for d in DERS}
            for d in DERS:
                outc[d] += s*cp[d]*P4[e1]*P4[e2]
                oute[d] += s*cp[(0, 0)]*D4[e1][d]*P4[e2]
    return outc, oute, n


# ---------------------------------------------------------------- 读数的团簇展开
QW, WF = sp.field("w", sp.QQ)  # 有理函数域 Q(w)，w = v/u
w = sp.Symbol("w")
ORD = (0, 4, 6, 8)


def to_field(p):
    return sum((QW(sp.Rational(int(c.p), int(c.q)))*WF**i for i, c in enumerate(p.coeffs())), QW(0))


class Jet:
    """点 0 两个耦合 (u0, v0) 上截断到二阶的泰勒式，分量依次为 f, f_u, f_v, f_uv, f_uu, f_vv。"""
    def __init__(self, c):
        self.c = list(c)

    @staticmethod
    def const(x):
        return Jet([x]+[QW(0)]*5)

    def __add__(self, o):
        o = o if isinstance(o, Jet) else Jet.const(o)
        return Jet([x+y for x, y in zip(self.c, o.c)])
    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, Jet) else Jet.const(o)
        return Jet([x-y for x, y in zip(self.c, o.c)])

    def __rmul__(self, k):
        return Jet([k*x for x in self.c])

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return Jet([o*x for x in self.c])
        f, fu, fv, fuv, fuu, fvv = self.c
        g, gu, gv, guv, guu, gvv = o.c
        return Jet([f*g, fu*g+f*gu, fv*g+f*gv, fuv*g+fu*gv+fv*gu+f*guv,
                    fuu*g+2*fu*gu+f*guu, fvv*g+2*fv*gv+f*gvv])

    def inv(self):
        f, fu, fv, fuv, fuu, fvv = self.c
        return Jet([1/f, -fu/f**2, -fv/f**2, 2*fu*fv/f**3-fuv/f**2,
                    2*fu**2/f**3-fuu/f**2, 2*fv**2/f**3-fvv/f**2])

    def __truediv__(self, o):
        return self*(o.inv() if isinstance(o, Jet) else Jet.const(1/o))


def jet(V):
    return Jet([V[d] for d in [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]])


def log_jet(J):
    f, fu, fv, fuv, fuu, fvv = J.c
    return Jet([None, fu/f, fv/f, fuv/f-fu*fv/f**2, fuu/f-fu**2/f**2, fvv/f-fv**2/f**2])


def cluster_data(AMP, rho, uv=1):
    raw = {"s0": {d: amp_poly(AMP[(0, 0)], rho, d, uv) for d in DERS},
           "b4": two_point(AMP, rho, 4, uv), "b8": two_point(AMP, rho, 8, uv)}
    raw["tri"], ntri = triangle_t6(AMP, rho, uv)
    raw["pc"], raw["pe"], nch = chain_joint(AMP, rho, uv)
    return {k: {d: to_field(v[d]) for d in DERS} for k, v in raw.items()}, ntri, nch


def series_readouts(Dd):
    s, s0 = Dd["s0"][(0, 0)], jet(Dd["s0"])
    b4, b8 = jet(Dd["b4"])/(s0*s), jet(Dd["b8"])/(s0*s)
    b4u, b8u = Dd["b4"][(0, 0)]/s**2, Dd["b8"][(0, 0)]/s**2
    g6 = jet(Dd["tri"])/(s0*s**2)
    pic, pie = jet(Dd["pc"])/(s0*s**2), jet(Dd["pe"])/(s0*s**2)
    L0 = log_jet(s0)
    c8 = b8-(b4*b4)*QW(sp.Rational(1, 2))
    graphs = {
        # log Z 减去与点 0 无关的常数，到 t^8。单键 t^8 的连通部分为 b8 - b4^2/2，
        # 共点两键的连通部分为（联合项）- b4 b4'；三角形的 t^6 项只来自三角形本身。
        "two-point": (1, {4: b4, 6: None, 8: c8}),
        "triangle": (2, {4: 2*b4+b4u, 6: g6,
                         8: 2*c8+(b8u-b4u**2*QW(sp.Rational(1, 2)))+(pic-b4*b4)+2*(pie-b4*b4u)}),
        "chain-center": (2, {4: 2*b4, 6: None, 8: 2*c8+(pic-b4*b4)}),
        "chain-end": (1, {4: b4+b4u, 6: None, 8: c8+(b8u-b4u**2*QW(sp.Rational(1, 2)))+(pie-b4*b4u)}),
    }
    out = {}
    for name, (z, comps) in graphs.items():
        dd = {0: L0.c}
        for k in (4, 6, 8):
            dd[k] = comps[k].c if comps[k] is not None else [QW(0)]*6
        # 分量下标：1=F_u, 2=F_v, 3=F_uv, 4=F_uu, 5=F_vv
        def ratio(nk, d1k, d2k):
            n = {k: dd[k][nk] for k in ORD}
            d1 = {k: dd[k][d1k] for k in ORD}
            d2 = {k: dd[k][d2k] for k in ORD}
            e1 = {k: d1[k]/d1[0] for k in (4, 6, 8)}
            e2 = {k: d2[k]/d2[0] for k in (4, 6, 8)}
            inv = lambda e: {0: QW(1), 4: -e[4], 6: -e[6], 8: -e[8]+e[4]**2}
            i1, i2 = inv(e1), inv(e2)
            res = {}
            for k in range(9):
                acc = QW(0)
                for a_ in ORD:
                    for b_ in ORD:
                        c_ = k-a_-b_
                        if c_ in ORD:
                            acc += n[a_]*i1[b_]*i2[c_]
                res[k] = acc/(d1[0]*d2[0])
            res[0] = res[0]+1
            return res
        out[name] = (z, [ratio(3, 1, 2), ratio(4, 1, 1), ratio(5, 2, 2)])
    return out


# ---------------------------------------------------------------- 两点格的符号首阶数据
def symbolic_t4():
    """各点同参的两点格：静态读数 S_i(rho) 与 t^4 系数 g_i(u,v,rho)。"""
    Z4 = _hop["Z_of"](_hop["generating_series"](), 4)
    u0, u1, v0, v1, x, y, t = (_hop[k] for k in ("u0", "u1", "v0", "v1", "x", "y", "t"))
    u, v, R = sp.symbols("u v rho", positive=True)
    z = Z4.subs(y, 1)
    uni = {u0: u, u1: u, v0: v, v1: v, x: R*u*v}
    co = lambda e, k: sp.expand(sp.Poly(e, t).coeff_monomial(t**k).subs(uni))
    P = {}
    for n, ds in {"Z": [], "a": [u0], "b": [v0], "ab": [u0, v0], "aa": [u0, u0], "bb": [v0, v0]}.items():
        e = z
        for dd in ds:
            e = sp.diff(e, dd)
        P[n] = (co(e, 0), co(e, 4))
    eps = sp.Symbol("eps")
    ser = lambda n: P[n][0]+eps*P[n][1]
    X = {"AB": ser("Z")*ser("ab")/(ser("a")*ser("b")), "AA": ser("Z")*ser("aa")/ser("a")**2,
         "BB": ser("Z")*ser("bb")/ser("b")**2}
    S = {k: sp.cancel(e.subs(eps, 0)) for k, e in X.items()}
    g = {k: sp.cancel(sp.diff(e, eps).subs(eps, 0)) for k, e in X.items()}
    return Z4, (u, v, R), S, g


def main():
    t0 = time.time()
    P, E = shell_graph()
    Z4, (u, v, R), S, g = symbolic_t4()
    p = 1+24*R+90*R**2+80*R**3+15*R**4
    cub = 1+18*R+45*R**2+20*R**3
    Wp = 1+12*R+225*R**2+2030*R**3+7455*R**4+17640*R**5+22925*R**6+16050*R**7+6750*R**8+1000*R**9
    Qf = ((16+288*R+720*R**2+320*R**3)*u**2-(32+240*R+320*R**2+80*R**3)*u*v+(1+18*R+45*R**2+20*R**3)*v**2)
    Sp = {k: sp.diff(S[k], R) for k in KEYS}
    lamq = -R*Qf/(192*u**2*v**2*p)
    # 定理 A：静态 X_AA = X_BB；t^4 修正等于有效重叠的平移 rho -> rho + z t^4 lamq（对一切 rho）。
    assert sp.simplify(S["AA"]-S["BB"]) == 0
    assert S["AB"].subs(R, 1) == sp.Rational(5, 4) and S["AB"].subs(R, 0) == 1
    assert sp.simplify(S["AA"]-3*(1+12*R+15*R**2)*p/(4*cub**2)) == 0
    assert sp.simplify(Sp["AB"]-3*Wp/(2*cub**3)) == 0
    for k in KEYS:
        assert sp.simplify(g[k]-lamq*Sp[k]) == 0
    print(f"[{time.time()-t0:.0f}s] 符号 t^4 数据：g_i = lamq * S_i'(rho)")

    AMP = site_amplitudes()
    print(f"[{time.time()-t0:.0f}s] 单点振幅 {len(AMP)} 个缺失集合")
    # 校验：两点格 t^0,t^4,t^8 系数与 225 项多项式（rho=1/2，各导数版本）。
    u0, u1, v0, v1, x, y, t = (_hop[k] for k in ("u0", "u1", "v0", "v1", "x", "y", "t"))
    half = Fr(1, 2)
    for loss in (0, 4, 8):
        mine = two_point(AMP, half, loss, 1)
        ck = sp.Poly(Z4, t).coeff_monomial(t**loss)
        for d in DERS:
            e = ck
            if d[0]:
                e = sp.diff(e, u0, d[0])
            if d[1]:
                e = sp.diff(e, v0, d[1])
            e = sp.expand(e.subs({u0: 1, u1: 1, v0: w, v1: w, x: sp.Rational(1, 2)*w, y: 1}))
            assert QW.from_expr(e) == to_field(mine[d])
    print(f"[{time.time()-t0:.0f}s] 两点格 t^0,t^4,t^8 与 225 项多项式逐项相等")

    results = {}
    for rho in (Fr(1, 2), Fr(1, 3), Fr(0)):
        Dd, ntri, nch = cluster_data(AMP, rho)
        assert (ntri, nch) == (3216, 17880)
        results[rho] = series_readouts(Dd)
    # 齐次性：rho=1/2，u=2 处 R_{2k}(2,2w) = 2^{-k} R_{2k}(1,w)。
    Dd2, _, _ = cluster_data(AMP, half, uv=2)
    r2 = series_readouts(Dd2)
    for name in r2:
        for i in range(3):
            for k, deg in [(4, -2), (6, -3), (8, -4)]:
                assert r2[name][1][i][k] == results[half][name][1][i][k]*QW(sp.Rational(2)**deg)
    print(f"[{time.time()-t0:.0f}s] 三角形 t^6、开链 t^8 联合项与读数级数完成；齐次性核对通过")

    F = lambda e: QW.from_expr(sp.sympify(e))
    cross = lambda a_, b_: [a_[1]*b_[2]-a_[2]*b_[1], a_[2]*b_[0]-a_[0]*b_[2], a_[0]*b_[1]-a_[1]*b_[0]]
    det3 = lambda a_, b_, c_: sum(x*y for x, y in zip(cross(b_, c_), a_))
    J16 = {}
    for rho in (Fr(1, 2), Fr(1, 3)):
        rs = sp.Rational(rho.numerator, rho.denominator)
        Spv = [F(Sp[k].subs(R, rs)) for k in KEYS]
        for name, (z, ser) in results[rho].items():
            f = {k: [ser[i][k] for i in range(3)] for k in range(9)}
            assert all(x == 0 for k in (1, 2, 3, 5, 7) for x in f[k])
            assert f[0] == [F(S[k].subs(R, rs)) for k in KEYS]
            R4 = [z*g[k] for k in KEYS]
            sub = {u: 1, v: w, R: rs}
            assert f[4] == [F(x.subs(sub)) for x in R4]
            # t^6：三角形上非零、与 S' 平行（仍是有效重叠的平移）且 X_AA = X_BB；其余图上为零。
            if name == "triangle":
                assert any(x != 0 for x in f[6]) and all(x == 0 for x in cross(f[6], Spv)) and f[6][1] == f[6][2]
            else:
                assert all(x == 0 for x in f[6])
            # X_BB - X_AA 的首项在 t^8，且不是零。
            assert f[8][2] != f[8][1]
            # 雅可比的列：rho 列 = S' + t^4 dR4/drho + t^6 (?) + ...，u、v 列从 t^4 起。
            # 已知：rho 列的 0、2（为零）、4 阶；u、v 列的 4、6、8 阶（u 方向由欧拉关系 u d_u = deg - w d_w）。
            rcol = {0: Spv, 2: [QW(0)]*3, 4: [F(sp.diff(x, R).subs(sub)) for x in R4]}
            ucol = {4: [F(sp.diff(x, u).subs(sub)) for x in R4]}
            vcol = {4: [F(sp.diff(x, v).subs(sub)) for x in R4]}
            for k, deg in [(6, -3), (8, -4)]:
                vcol[k] = [x.diff(WF) for x in f[k]]
                ucol[k] = [QW(deg)*x-WF*y for x, y in zip(f[k], vcol[k])]
            par = lambda a_, b_: all(x == 0 for x in cross(a_, b_))
            # u、v 两列在 4、6 阶都与 S' 平行（有效重叠结构）。
            for k in (4, 6):
                assert par(ucol[k], Spv) and par(vcol[k], Spv)
            coef = {}
            for n in range(8, 17, 2):
                acc = QW(0)
                for a_ in range(0, n+1, 2):
                    for b_ in range(4, n+1, 2):
                        c_ = n-a_-b_
                        if c_ < 4:
                            continue
                        known = (a_ in rcol, b_ in ucol, c_ in vcol)
                        if all(known):
                            acc += det3(rcol[a_], ucol[b_], vcol[c_])
                        else:
                            # 恰有一列未算出（n<=16 时不会有两列）：另两列平行，该项为零。
                            assert sum(known) == 2
                            others = [cl for cl, kn in zip((rcol.get(a_), ucol.get(b_), vcol.get(c_)), known) if kn]
                            assert par(*others), (name, rho, n, a_, b_, c_)
                coef[n] = acc
            for n in range(8, 16, 2):
                assert coef[n] == 0, (name, rho, n)
            assert coef[16] != 0
            J16[rho, name] = (z, coef[16])
    for rho in (Fr(1, 2), Fr(1, 3)):
        base = J16[rho, "two-point"][1]
        for name in ("triangle", "chain-center", "chain-end"):
            z, c = J16[rho, name]
            assert c == QW(z**2)*base
    j2 = J16[half, "two-point"][1]
    poly6 = (8555769273768192*w**6-188484971655861335*w**5-5571824064272581248*w**4-20828682148361300688*w**3
             - 89149185028361299968*w**2-48252152743900501760*w+35044430945354514432)
    assert j2 == F(-poly6/(46083730603479903436800000000*w**8))
    j2_1 = sp.sympify(j2.as_expr()).subs(w, 1)
    assert j2_1 == sp.Rational(1031498737935386099, 368669844827839227494400000)
    assert sp.Rational(2797, 10**12) < j2_1 < sp.Rational(2798, 10**12)
    j3_1 = sp.sympify(J16[Fr(1, 3), "two-point"][1].as_expr()).subs(w, 1)
    assert j3_1 == sp.Rational(10593838074357783855, 1380913646891877604595335168)
    print(f"[{time.time()-t0:.0f}s] 雅可比：t^8..t^14 恒为零；t^16 系数 = z^2 * J16(two-point)；"
          f"J16(1,1,1/2) = {float(j2_1):.6e}，J16(1,1,1/3) = {float(j3_1):.6e}")

    # 正交支 rho=0：到 t^8 为止 X_AB = 1 + z t^8/(4608 u^2 v^2)，X_AA = X_BB = 3/4（四种图）。
    for name in ("two-point", "triangle", "chain-center", "chain-end"):
        z, ser = results[Fr(0)][name]
        assert ser[0][0] == 1 and ser[1][0] == QW(sp.Rational(3, 4)) and ser[2][0] == QW(sp.Rational(3, 4))
        for k in range(1, 8):
            assert all(ser[i][k] == 0 for i in range(3))
        assert ser[0][8] == QW(z)/(4608*WF**2) and ser[1][8] == 0 and ser[2][8] == 0

    # 两点格完整有理函数：静态射线 U=V=L 上的精确雅可比与 J16 L^{-10} 一致。
    Pr = _hop["readouts"](Z4)
    U, V, RR = _hop["U"], _hop["V"], _hop["R"]
    num = {"AB": Pr["Z"]*Pr["Zab"], "AA": Pr["Z"]*Pr["Zaa"], "BB": Pr["Z"]*Pr["Zbb"]}
    den = {"AB": Pr["Za"]*Pr["Zb"], "AA": Pr["Za"]**2, "BB": Pr["Zb"]**2}
    dn = {k: [num[k].diff(q) for q in (U, V, RR)] for k in num}
    dd = {k: [den[k].diff(q) for q in (U, V, RR)] for k in den}

    def Jex(pt):
        M = []
        for k in KEYS:
            n0, d0 = num[k](*pt), den[k](*pt)
            M.append([(dn[k][i](*pt)*d0-n0*dd[k][i](*pt))/d0**2 for i in range(3)])
        return sp.Matrix(M).det()
    j111 = Jex((1, 1, sp.Rational(1, 2)))
    assert sp.Rational(21, 10**10) < j111 < sp.Rational(22, 10**10)  # 跳跃读数命题 7 的值，重新断言
    for L, tol in [(10**3, sp.Rational(1, 10**6)), (10**4, sp.Rational(1, 10**8))]:
        jv = Jex((L, L, sp.Rational(1, 2)))
        assert jv > 0 and abs(jv*L**10-j2_1) < tol*j2_1
    print(f"[{time.time()-t0:.0f}s] 两点格 L^10 J(L,L,1/2) 在 L=10^3,10^4 与 J16 的相对差分别 <1e-6, <1e-8")
    print(f"十九点壳：{len(P)} 点；G19 有 135 条边、529 个三角形；只保留最短弦（平方 12）的子图为三条互不相交的边 {E}")
    print(f"用时 {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()

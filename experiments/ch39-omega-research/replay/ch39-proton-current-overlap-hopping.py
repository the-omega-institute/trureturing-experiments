"""辅助电流重叠 rho=|c|^2 的 Q-L 接触探针加入最近邻跳跃 t 之后：两点格的精确配分多项式、
读数比 X 的归一依赖、端值窗口与反例，不含三角形的图上 O(t^4) 首阶修正的成键分解，以及三角形上的 t^6 项。
配套说明：最近邻跳跃下接触读数比随Mt漂移与三比值局部确定重叠模长.md

探针（人为定义，不是模型配分函数的极限）：每格点一代 Q,L 及其共轭，单点接触项
kappa H_x，H = sum_a [K11 Dbar_a D_a + K22 (QL)_a T_a + K21 (QL)_a D_a + K12 Dbar_a T_a]，
最近邻跳跃 t sum_<xy> sum_k (psibar_{x,k} psi_{y,k} + psibar_{y,k} psi_{x,k})，k 取遍 Q,L 全部 16 个分量
（色、弱同位旋、旋量指标逐分量对角，Q 与 L 同一个 t）。Z = Berezin 顶阶系数。

方法：(1) 把 (i,alpha) 指标上的对称形 g=eps(x)eps 用实正交变换与相位化为单位阵，模型等价于 n_f=4 个
"味" 拷贝，每个拷贝每格点 4 个非共轭生成元 (Q^1,Q^2,Q^3,L) 与 4 个共轭生成元；跳跃核不变。
(2) 以形式高斯（Wick 矩）把接触项线性化，每个味的格拉斯曼积分是一个多项式 f(J)，Z=E[f^{n_f}]。
(3) 两点格上 f=P01 P10，每色的矩母函数为 det(I-SB)^{-1/2}=1/Delta（det 为完全平方），Z 由
1/Delta^3 的截断级数精确抽取。n_f=1,2 与不经线性化的直接格拉斯曼积分逐点精确比对；
n_f=4 的静态极限与上一份脚本的多项式 2211840 (K11K22)^4 p(rho) 逐项一致。
(4) 另在原始 eps(x)eps 基下（不经味约化、不经线性化）直接计算物理 n_f=4 两点格：单点元素
E_u(d,l)=prod_a D_a^{d_a}(QL)_a^{l_a} 与共轭一侧 E_b，跳跃把点 0 的 (U,B) 唯一配到点 1 的 (B,U)，
符号为 (-1)^{m(m-1)/2+n(n-1)/2}（m,n 为两方向跳跃数），逐阶 t^0..t^32 与 (2)(3) 的多项式精确相等。
全部为精确整数/有理运算。
"""
import runpy
import time
from collections import defaultdict
from fractions import Fraction as Fr
from math import factorial as fa
from pathlib import Path

import flint
import sympy as sp

HERE = Path(__file__).resolve().parent
_prev = runpy.run_path(str(HERE / "ch39-proton-current-overlap-static.py"))
word, mul = _prev["word"], _prev["mul"]

EPS3 = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


# ---------------------------------------------------------------- 格拉斯曼直接积分（味约化基）
def add(A, B):
    out = defaultdict(int, A)
    for k, v in B.items():
        out[k] += v
    return {k: v for k, v in out.items() if v != 0}


def exp_even(S):
    r = {0: 1}
    for m, c in S.items():
        r = mul(r, {0: 1, m: c})
    return r


class Lat:
    def __init__(self, L, nf):
        self.L, self.nf, self.N = L, nf, 8*L*nf

    def u(self, x, A, mu):
        return ((x*self.nf+A)*4+mu)*2

    def b(self, x, A, mu):
        return ((x*self.nf+A)*4+mu)*2+1

    def bil(self, x, kind):
        """非共轭 D_a=sum eps_abc psi_b psi_c, (QL)_a=psi_a psi_L；共轭一侧取反序共轭。"""
        out = []
        for a in range(3):
            d = defaultdict(int)
            for A in range(self.nf):
                if kind in ("D", "Db"):
                    for (aa, b, c), e in EPS3.items():
                        if aa != a:
                            continue
                        pair = [self.u(x, A, b), self.u(x, A, c)] if kind == "D" else [self.b(x, A, c), self.b(x, A, b)]
                        m, s = word(pair)
                        d[m] += e*s
                else:
                    pair = [self.u(x, A, a), self.u(x, A, 3)] if kind == "QL" else [self.b(x, A, 3), self.b(x, A, a)]
                    m, s = word(pair)
                    d[m] += s
            out.append({k: v for k, v in d.items() if v})
        return out

    def hop(self, x, y, t):
        d = {}
        for A in range(self.nf):
            for mu in range(4):
                m, s = word([self.b(x, A, mu), self.u(y, A, mu)])
                d[m] = s*t
        return d

    def contact(self, x, par):
        D, QL, Db, T = [self.bil(x, k) for k in ("D", "QL", "Db", "T")]
        u, v, xx, yy = par
        S = {}
        for a in range(3):
            for (A, B, c) in [(Db, D, u), (QL, T, v), (QL, D, xx), (Db, T, yy)]:
                S = add(S, {k: c*w for k, w in mul(A[a], B[a]).items()})
        return S


def direct_chain(L, nf, pars, bonds, t):
    lat = Lat(L, nf)
    S = {}
    for x in range(L):
        S = add(S, lat.contact(x, pars[x]))
    for i, j in bonds:
        S = add(S, lat.hop(i, j, t))
        S = add(S, lat.hop(j, i, t))
    return exp_even(S).get((1 << lat.N)-1, 0)


def direct_two_site(nf, par0, par1, t):
    """两点格：先展开单点 exp(kappa H_x)，跳跃把点 0 的单项唯一地配到点 1 的单项。"""
    lat = Lat(2, nf)
    P0, P1 = exp_even(lat.contact(0, par0)), exp_even(lat.contact(1, par1))
    idx = [(A, mu) for A in range(nf) for mu in range(4)]
    tot = 0
    for m0, c0 in P0.items():
        Us = {k for k in idx if m0 >> lat.u(0, *k) & 1}
        Bs = {k for k in idx if m0 >> lat.b(0, *k) & 1}
        m1 = sum(1 << lat.u(1, *k) for k in Bs) | sum(1 << lat.b(1, *k) for k in Us)
        c1 = P1.get(m1, 0)
        if c1 == 0:
            continue
        prod = mul({m0: c0}, {m1: c1})
        for k in idx:
            if k not in Bs:
                m, s = word([lat.b(0, *k), lat.u(1, *k)])
                prod = mul(prod, {m: s*t})
            if k not in Us:
                m, s = word([lat.b(1, *k), lat.u(0, *k)])
                prod = mul(prod, {m: s*t})
        tot += prod.get((1 << lat.N)-1, 0)
    return tot


def z_original_basis(p0, p1):
    """原始基两点格的 Z，按 t 的幂返回 {幂: 系数}。

    推导：生成元按 [u0][b0][u1][b1] 排序。点 0 的单项 [u0^U][b0^B] 迫使点 1 取 [u1^B][b1^U]，
    跳跃补齐余集：prod_{k notin B} b0_k u1_k = eps(n)[b0^Bc][u1^Bc]，prod_{k notin U} b1_k u0_k = eps(m)[b1^Uc][u0^Uc]，
    eps(n)=(-1)^{n(n-1)/2}。各块都是偶次，可整块交换；每个生成元族的 S 与余集合并给出同一符号 sigma(S)，
    它在 U、B 上各出现两次而相消，故总符号只剩 eps(m) eps(n)。"""
    D, QL = _prev["bilinears"](False)
    Db, T = _prev["bilinears"](True)

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
    Eu, Eb = elements(D, QL), elements(Db, T)
    size = lambda k: 2*sum(k)
    M = {}
    for a, ea in Eu.items():
        for b, fb in Eb.items():
            if size(a) == size(b):
                v = sum(c*fb.get(m, 0) for m, c in ea.items())
                if v:
                    M[a, b] = v

    def W(par, a, b):
        u, v, xx, yy = par
        tot = Fr(1)
        for c in range(3):
            d, l, db, lb = a[2*c], a[2*c+1], b[2*c], b[2*c+1]
            sc = Fr(0)
            for n11 in range(5):
                n21, n12, n22 = d-n11, db-n11, l-d+n11
                if min(n21, n12, n22) < 0 or n22+n12 != lb:
                    continue
                sc += u**n11*v**n22*xx**n21*yy**n12/(fa(n11)*fa(n22)*fa(n21)*fa(n12))
            if sc == 0:
                return 0
            tot *= sc
        return tot
    W0 = {(a, b): w for a in Eu for b in Eb if (w := W(p0, a, b))}
    W1col = defaultdict(list)
    for a in Eu:
        for b in Eb:
            w = W(p1, a, b)
            if w:
                W1col[b].append((a, w))
    MW = defaultdict(Fr)
    for (a, bp), mv in M.items():
        for ap, wv in W1col[bp]:
            MW[a, ap] += mv*wv
    MWrow = defaultdict(list)
    for (a, ap), v in MW.items():
        MWrow[a].append((ap, v))
    out = defaultdict(Fr)
    for (a, b), w in W0.items():
        m, n = 16-size(a), 16-size(b)
        sg = (-1)**(m*(m-1)//2+n*(n-1)//2)
        acc = sum((v*M.get((ap, b), 0) for ap, v in MWrow[a]), Fr(0))
        if acc:
            out[m+n] += sg*w*acc
    return dict(out)


# ---------------------------------------------------------------- 线性化后的两点格精确多项式
u0, u1, v0, v1, x, y, t = sp.symbols("u0 u1 v0 v1 x y t")
PARS = sp.symbols("s0 sb1 r s1 sb0 rp")


def per_flavor_factor():
    """n_f=1 两点格（符号源 J）的直接积分 = P01*P10，并读出 P01 的形式。"""
    lat = Lat(2, 1)
    S, J = {}, {}
    for xs in range(2):
        for kind in ("D", "QL", "Db", "T"):
            bl = lat.bil(xs, kind)
            for a in range(3):
                J[kind, xs, a] = sp.Symbol(f"{kind}{xs}{a}")
                S = add(S, {k: c*J[kind, xs, a] for k, c in bl[a].items()})
    S = add(S, lat.hop(0, 1, t))
    S = add(S, lat.hop(1, 0, t))
    f = sp.expand(exp_even(S).get((1 << lat.N)-1, 0))
    dot = lambda p, q, x0, x1: sum(J[p, x0, a]*J[q, x1, a] for a in range(3))
    P01 = 4*dot("D", "QL", 0, 0)*dot("Db", "T", 1, 1)+t**2*(4*dot("D", "Db", 0, 1)+dot("QL", "T", 0, 1))+t**4
    P10 = 4*dot("D", "QL", 1, 1)*dot("Db", "T", 0, 0)+t**2*(4*dot("D", "Db", 1, 0)+dot("QL", "T", 1, 0))+t**4
    assert sp.expand(f-P01*P10) == 0
    return f


def generating_series():
    """每色矩母函数 det(I-SB)^{-1/2} = 1/Delta，三色得 1/Delta^3；截断展开。"""
    s0, sb1, r, s1, sb0, rp = PARS
    S = sp.zeros(8)
    for xs, (uu, vv) in enumerate([(u0, v0), (u1, v1)]):
        J, Jp, Jb, Jbp = [4*xs+k for k in range(4)]
        for i, j, val in [(J, Jb, uu), (Jp, Jbp, vv), (J, Jp, x), (Jb, Jbp, y)]:
            S[i, j] = S[j, i] = val
    B = sp.zeros(8)

    def put(i, j, c):
        B[i, j] += c
        B[j, i] += c
    put(0, 1, s0); put(6, 7, sb1); put(0, 6, 4*r); put(1, 7, r)
    put(4, 5, s1); put(2, 3, sb0); put(4, 2, 4*rp); put(5, 3, rp)
    d = sp.expand((sp.eye(8)-S*B).det())
    fl = sp.factor_list(d)
    assert len(fl[1]) == 1 and fl[1][0][1] == 2
    Delta = fl[1][0][0]
    if Delta.subs({p: 0 for p in PARS}) == -1:
        Delta = -Delta
    assert Delta.subs({p: 0 for p in PARS}) == 1 and sp.expand(Delta**2-d) == 0
    delta = sp.expand(1-Delta)
    KV = ["u0", "u1", "v0", "v1", "x", "y"]
    ctx = flint.fmpz_mpoly_ctx.get(KV, "lex")
    dterms = [(m[:6], ctx.from_dict({m[6:]: int(c)}))
              for m, c in sp.Poly(delta, *PARS, u0, u1, v0, v1, x, y).terms()]
    ok = lambda e: e[0]+e[2] <= 4 and e[1]+e[2] <= 4 and e[3]+e[5] <= 4 and e[4]+e[5] <= 4
    one = ctx.from_dict({(0,)*6: 1})
    G, cur = {(0,)*6: one}, {(0,)*6: one}
    for n in range(1, 17):
        new = {}
        for e, c in cur.items():
            for pe, m in dterms:
                e2 = tuple(a+b for a, b in zip(e, pe))
                if ok(e2):
                    new[e2] = new[e2]+c*m if e2 in new else c*m
        cur = {k: v for k, v in new.items() if v != 0}
        for k, v in cur.items():
            G[k] = G.get(k, 0*one)+((n+1)*(n+2)//2)*v
    syms = dict(zip(KV, [u0, u1, v0, v1, x, y]))
    return {k: sp.sympify(str(v), locals=syms) for k, v in G.items()}


def Z_of(G, nf):
    tot = 0
    for i in range(nf+1):
        for j in range(nf+1-i):
            for i2 in range(nf+1):
                for j2 in range(nf+1-i2):
                    key = (i, i, j, i2, i2, j2)
                    if key not in G:
                        continue
                    k, k2 = nf-i-j, nf-i2-j2
                    m = Fr(fa(nf), fa(i)*fa(j)*fa(k))*Fr(fa(nf), fa(i2)*fa(j2)*fa(k2))
                    tot += (sp.Rational(m.numerator, m.denominator)*4**(i+i2)*t**(2*j+4*k+2*j2+4*k2)
                            * fa(i)**2*fa(j)*fa(i2)**2*fa(j2)*G[key])
    return sp.expand(tot)


def flavor_reduction():
    """g_{(i a),(j b)}=eps_ij eps_ab：存在 S=O diag(1,1,i,i) 使 S^T g S = I 且 conj(S)^T S = I。"""
    g = sp.Matrix(4, 4, lambda A, B: ({(0, 1): 1, (1, 0): -1}.get((A//2, B//2), 0)
                                      * {(0, 1): 1, (1, 0): -1}.get((A % 2, B % 2), 0)))
    assert g == g.T and g*g == sp.eye(4)
    O = sp.Matrix([[1, 0, 0, 1], [0, 1, -1, 0], [0, 1, 1, 0], [1, 0, 0, -1]]).T/sp.sqrt(2)  # 列为本征向量
    assert sp.simplify(O.T*O-sp.eye(4)) == sp.zeros(4)
    assert sp.simplify(O.T*g*O-sp.diag(1, 1, -1, -1)) == sp.zeros(4)
    Sm = O*sp.diag(1, 1, sp.I, sp.I)
    assert sp.simplify(Sm.T*g*Sm-sp.eye(4)) == sp.zeros(4)
    assert sp.simplify(Sm.conjugate().T*Sm-sp.eye(4)) == sp.zeros(4)


# ---------------------------------------------------------------- 读数
U, V, R = sp.symbols("U V rho")


def readouts(z4):
    """t=1、y=1、x=rho*U*V（即 K21K12=rho K11K22），各点同参后的 Z 与点 0 的导数。"""
    zz = z4.subs({y: 1, t: 1})
    uni = {u0: U, u1: U, v0: V, v1: V, x: R*U*V}
    out = {}
    for name, ders in [("Z", []), ("Za", [u0]), ("Zb", [v0]), ("Zab", [u0, v0]),
                       ("Zaa", [u0, u0]), ("Zbb", [v0, v0])]:
        e = zz
        for dv in ders:
            e = sp.diff(e, dv)
        out[name] = sp.Poly(sp.expand(e.subs(uni)), U, V, R)
    return out


def lead(p, a, b):
    w = max(a*m[0]+b*m[1] for m in p.monoms())
    return sum(c*U**m[0]*V**m[1]*R**m[2] for m, c in p.terms() if a*m[0]+b*m[1] == w)


def main():
    t0 = time.time()
    flavor_reduction()
    per_flavor_factor()
    G = generating_series()
    Z1, Z2, Z4 = Z_of(G, 1), Z_of(G, 2), Z_of(G, 4)
    # n_f=1：符号核对；n_f=1,2：与不经线性化的直接积分逐点精确比对（点 0、1 参数不同，t 不为 1）。
    assert Z1 == sp.expand(t**8+48*t**4*u0*u1+3*t**4*v0*v1+72*t**4*x*y+(12*u0*v0+36*x*y)*(12*u1*v1+36*x*y))
    pts = [((Fr(2), Fr(5), Fr(11), Fr(13)), (Fr(3), Fr(7), Fr(11), Fr(13)), Fr(1)),
           ((Fr(-1, 2), Fr(5, 3), Fr(1, 3), Fr(-4)), (Fr(3, 5), Fr(2), Fr(1, 3), Fr(-4)), Fr(3, 2))]
    for nf, Zs in [(1, Z1), (2, Z2)]:
        for p0, p1, tt in pts:
            sub = {u0: p0[0], v0: p0[1], u1: p1[0], v1: p1[1], x: p0[2], y: p0[3], t: tt}
            hz = Zs.subs({k: sp.Rational(v.numerator, v.denominator) for k, v in sub.items()})
            dz = direct_two_site(nf, p0, p1, tt)
            assert sp.Rational(dz.numerator, dz.denominator) == hz
    # 静态极限：t^0 系数 = 两个单点静态不变量之积，常数与上一份脚本相同。
    rho = sp.Symbol("rho")
    p = 1+24*rho+90*rho**2+80*rho**3+15*rho**4
    site = lambda uu, vv: 2211840*uu**4*vv**4*p.subs(rho, x*y/(uu*vv))
    assert sp.expand(Z4.subs(t, 0)-site(u0, v0)*site(u1, v1)) == 0
    # 原始基、物理 n_f=4、不经味约化与线性化的直接计算：t^0..t^32 逐阶精确相等。
    for p0, p1 in [((Fr(2), Fr(5), Fr(11), Fr(13)), (Fr(3), Fr(7), Fr(11), Fr(13))),
                   ((Fr(-1, 2), Fr(5, 3), Fr(1, 3), Fr(-4)), (Fr(3, 5), Fr(2), Fr(1, 3), Fr(-4)))]:
        zo = z_original_basis(p0, p1)
        sub = {u0: p0[0], v0: p0[1], u1: p1[0], v1: p1[1], x: p0[2], y: p0[3]}
        Pz = sp.Poly(Z4.subs({k: sp.Rational(v.numerator, v.denominator) for k, v in sub.items()}), t)
        assert set(zo) == set(range(0, 33, 4))
        for k in range(33):
            c = zo.get(k, 0)
            assert (sp.Rational(c.numerator, c.denominator) if c else 0) == Pz.coeff_monomial(t**k)
    # t 只以 t^{4k} 出现；Z 只经 x*y 依赖 K21,K12（c 的相位不进入）。
    tp = sorted({m[0] for m in sp.Poly(Z4, t).monoms()})
    assert tp == [0, 4, 8, 12, 16, 20, 24, 28, 32]
    lam = sp.Symbol("lam")
    assert sp.expand(Z4.subs({x: lam*x, y: y/lam}, simultaneous=True)-Z4) == 0
    assert len(Z4.args) == 225

    # ---- 首阶修正（任意 u,v,rho；t 保留）
    uu, vv = sp.symbols("u v", positive=True)
    uniT = {u0: uu, u1: uu, v0: vv, v1: vv, x: R*uu*vv}
    z = Z4.subs(y, 1)
    co = lambda e, k: sp.expand(sp.Poly(e, t).coeff_monomial(t**k).subs(uniT))
    parts = {n: [co(e, 0), co(e, 4)] for n, e in
             [("Z", z), ("Za", sp.diff(z, u0)), ("Zb", sp.diff(z, v0)), ("Zab", sp.diff(z, u0, v0))]}
    Xs = sp.cancel(parts["Z"][0]*parts["Zab"][0]/(parts["Za"][0]*parts["Zb"][0]))
    Xstatic = ((10*R**3+45*R**2+27*R+2)*p.subs(rho, R)/(2*(20*R**3+45*R**2+18*R+1)**2))
    assert sp.simplify(Xs-Xstatic) == 0
    dX = sp.cancel(Xs*(parts["Zab"][1]/parts["Zab"][0]+parts["Z"][1]/parts["Z"][0]
                       - parts["Za"][1]/parts["Za"][0]-parts["Zb"][1]/parts["Zb"][0]))
    Wpol = (1000*R**9+6750*R**8+16050*R**7+22925*R**6+17640*R**5+7455*R**4+2030*R**3+225*R**2+12*R+1)
    Qf = ((16+288*R+720*R**2+320*R**3)*uu**2-(32+240*R+320*R**2+80*R**3)*uu*vv
          + (1+18*R+45*R**2+20*R**3)*vv**2)
    dX_closed = -R*Wpol*Qf/(128*uu**2*vv**2*(20*R**3+45*R**2+18*R+1)**3*p.subs(rho, R))
    assert sp.simplify(dX-dX_closed) == 0
    # Q 在 rho=0 时为 16u^2-32uv+v^2（f=v/u=|F|^2 在 16-4sqrt15<f<16+4sqrt15 时为负）；rho=1 时为 84(v-4u)^2。
    assert sp.expand(Qf.subs(R, 0)-(16*uu**2-32*uu*vv+vv**2)) == 0
    assert sp.expand(Qf.subs(R, 1)-84*(vv-4*uu)**2) == 0
    assert set(sp.solve(Qf.subs({R: 0, uu: 1}), vv)) == {16-4*sp.sqrt(15), 16+4*sp.sqrt(15)}
    assert Qf.subs({R: 0, uu: 1, vv: 1}) == -15 and Qf.subs({R: 1, uu: 1, vv: 1}) == 756
    assert all(c > 0 for c in sp.Poly(Wpol, R).all_coeffs())
    # 成键分解：三点开链（n_f=1，各点参数不同）t^4 系数 = 两条键的两点修正之和。
    T = flint.fmpq_poly([0, 1])
    q = flint.fmpq
    pars = [(q(2), q(3), q(5), q(7)), (q(1, 2), q(-3), q(2, 3), q(1)), (q(4), q(1, 5), q(-1), q(3))]
    wrap = lambda P: [tuple(flint.fmpq_poly([c]) for c in pp) for pp in P]
    z3 = direct_chain(3, 1, wrap(pars), [(0, 1), (1, 2)], T)
    z01 = direct_chain(2, 1, wrap(pars[:2]), [(0, 1)], T)
    z12 = direct_chain(2, 1, wrap(pars[1:]), [(0, 1)], T)
    s = [direct_chain(1, 1, wrap([pp]), [], T)[0] for pp in pars]
    assert z3[2] == 0 and z3[0] == s[0]*s[1]*s[2]
    assert z3[4] == s[2]*z01[4]+s[0]*z12[4]
    # 三角形（三条键两两相邻）：t^6 项非零，开链上为零。成键分解只对不含三角形的图成立。
    z3t = direct_chain(3, 1, wrap(pars), [(0, 1), (1, 2), (0, 2)], T)
    assert z3[6] == 0 and z3t[6] == flint.fmpq(110034, 5) and z3t[6] != 0
    assert z3t[2] == 0 and z3t[0] == z3[0]

    # ---- t=1 的读数：X_AB, X_AA, X_BB 作为 (U,V,rho) 的有理函数
    P = readouts(Z4)
    X = lambda a, b, c: (P["Z"](a, b, c)*P["Zab"](a, b, c))/(P["Za"](a, b, c)*P["Zb"](a, b, c))
    # 极限（沿 U=lam^a, V=lam^b, lam->oo 的首项面）
    lim = {}
    for a, b in [(1, 1), (-1, -1), (1, -1), (-1, 1), (2, -1), (1, -2)]:
        L = {k: lead(P[k], a, b) for k in ("Z", "Za", "Zb", "Zab")}
        lim[a, b] = sp.factor(sp.cancel(L["Z"]*L["Zab"]/(L["Za"]*L["Zb"])))
    assert sp.simplify(lim[1, 1]-Xstatic) == 0
    assert sp.simplify(lim[-1, -1]-sp.Rational(3, 4)*(1+R)) == 0
    for k in [(1, -1), (-1, 1), (2, -1), (1, -2)]:
        assert lim[k] == 1
    # 非单调：U=1, V=1/10 时 X(0) > X(1/40)
    a0, a1 = X(1, sp.Rational(1, 10), 0), X(1, sp.Rational(1, 10), sp.Rational(1, 40))
    assert a0 > a1 and a0 > 1 > a1
    a2 = X(1, sp.Rational(1, 10), 1)
    assert a2 > a1 and sp.Rational(11194, 10000) < a2 < sp.Rational(11195, 10000)  # 先降后升：既非递增也非递减
    assert sp.Rational(10002, 10000) < a0 < sp.Rational(10003, 10000) and sp.Rational(9981, 10000) < a1 < sp.Rational(9982, 10000)
    # 两支窗口重叠：正交支在 (U,V)=(1/10,1) 的读数高于平行支在 (10^4,10^-4) 的读数
    xo, xp = X(sp.Rational(1, 10), 1, 0), X(10**4, sp.Rational(1, 10**4), 1)
    assert xo > sp.Rational(10069, 10000) and 1 < xp < sp.Rational(1000001, 1000000) and xo > xp
    # Z, Za, Zb 的系数全正 => 在 U,V>0、rho>=0 上分母不为零，X 连续；(0,oo)^2 连通，配合自由极限 3/4
    # 与静态极限 5/4，由介值定理两支都取遍 [1.000001, 1.0069]。
    for k in ("Z", "Za", "Zb"):
        assert all(c > 0 for c in P[k].coeffs())
    # 局部可识别性：(X_AB, X_AA, X_BB) 对 (U,V,rho) 的雅可比行列式非零；任两个读数不能确定 rho。
    num = {"AB": P["Z"]*P["Zab"], "AA": P["Z"]*P["Zaa"], "BB": P["Z"]*P["Zbb"]}
    den = {"AB": P["Za"]*P["Zb"], "AA": P["Za"]*P["Za"], "BB": P["Zb"]*P["Zb"]}
    jac = {}
    pt = (1, 1, sp.Rational(1, 2))
    for k in num:
        row = []
        for var in (U, V, R):
            n_, d_ = num[k], den[k]
            row.append((n_.diff(var)(*pt)*d_(*pt)-n_(*pt)*d_.diff(var)(*pt))/d_(*pt)**2)
        jac[k] = row
    Jm = sp.Matrix([jac["AB"], jac["AA"], jac["BB"]])
    detJ = Jm.det()
    assert sp.Rational(21, 10**10) < detJ < sp.Rational(22, 10**10)
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        assert Jm[[i, j], [0, 1]].det() != 0  # 纤维方向的 rho 分量非零
    print(f"Z(n_f=4) 共 {len(Z4.args)} 项；t 的幂 {tp}；静态极限 = 两点单点不变量之积")
    print("首阶修正 dX/t^4 =", sp.factor(dX_closed))
    print("极限：静态 -> X_static(rho)；自由 -> 3(1+rho)/4；U->oo,V->0 或反向 -> 1")
    print(f"非单调：X(1,1/10,0)={float(a0):.6f} > X(1,1/10,1/40)={float(a1):.6f} < X(1,1/10,1)={float(a2):.6f}")
    print(f"三角形 n_f=1 的 t^6 系数 {z3t[6]}（开链为 0）")
    print(f"窗口重叠：X_orth(1/10,1)={float(xo):.6f} > X_par(1e4,1e-4)={float(xp):.8f}")
    print(f"雅可比 det = {float(detJ):.6e}（精确有理，非零）")
    print(f"用时 {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()

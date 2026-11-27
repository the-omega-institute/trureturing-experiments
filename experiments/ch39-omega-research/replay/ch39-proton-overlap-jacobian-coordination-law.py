"""三个归一无关读数比的雅可比行列式在任意有限图上的首阶结构：t^16 系数为 z^2 J16^(2)(u,v,rho)
（z 为读数格点的配位数，u,v>0），一般参数处首个非零阶为 t^16；J16^(2) 作为 (u,v,rho) 的有理函数精确求出，对一切 rho 成立。应用于模型壳图 G19 的四类格点。
配套说明：三比值雅可比行列式在任意有限图上的t16系数为配位数平方乘两点格系数且一般参数处首个非零阶为t16由此判定G19四类格点.md

探针与 ch39-proton-current-overlap-hopping.py、ch39-proton-overlap-identifiability-graphs.py 相同（人为定义的读数）：
每格点一代 Q,L 及其共轭，单点接触项 kappa H_x，最近邻标量跳跃 t，各格点参数相同、各边跳跃相同。

全部为精确有理运算（python-flint 多元有理多项式），对符号 (u,v,rho,z) 成立，不在个别 rho 上取值：
(1) 引理：缺失 4 个生成元的 168 个单点振幅，每一个作为点 0 的扰动，读数的一阶响应都与 S'(rho) 平行
    （对 u,v,rho 恒等）。
(2) 星图 K_{1,z}（z 为符号）：Z 在点 0 的 6 分量导数组（未取对数）为 s0 + t^4 z b4 + t^8 (z b8 + C(z,2) pc)；
    读数 X = Z Z_xy/(Z_x Z_y) 的 t^4, t^8 系数是显式分母的有理函数。断言：t^8, t^12 系数为零，t^16 系数等于
    z^2 乘 z=1 的值；R8 的 z^2 系数与 l^2 S''/2 之差与 S' 平行（有效重叠的二阶平移）。
(3) J16^(2)(1,w,rho) 的因式分解、与已发表的 rho=1/2（作为 w 的有理函数）、rho=1/3（w=1）精确值一致，
    P(w,rho) 对 w 的系数无公共根，对角线 w=1 上的零点计数与隔离区间。
(4) 模型壳图 G19 的度分布 10,14,16,18 与对应的 t^16 系数倍数。
"""
import itertools
import json
import runpy
import time
from collections import Counter, defaultdict
from fractions import Fraction as Fr
from pathlib import Path

import flint
import sympy as sp

HERE = Path(__file__).resolve().parent
_g = runpy.run_path(str(HERE / "ch39-proton-overlap-identifiability-graphs.py"))
pc, sign_of, comp, bits, ff = _g["pc"], _g["sign_of"], _g["comp"], _g["bits"], _g["ff"]
DERS = [(0, 0), (1, 0), (0, 1), (1, 1), (2, 0), (0, 2)]
KEYS = ["AB", "AA", "BB"]
PAIRS = {"AB": ((0, 0), (1, 1), (1, 0), (0, 1)), "AA": ((0, 0), (2, 0), (1, 0), (1, 0)),
         "BB": ((0, 0), (0, 2), (0, 1), (0, 1))}
CTX = flint.fmpq_mpoly_ctx.get(("u", "v", "r", "z"), "lex")
U, V, RR, Z = CTX.gens()
ZERO = CTX.from_dict({})
Q = lambda a, b=1: flint.fmpq(a, b)


def amp3(acc, der=(0, 0)):
    """单点振幅作为 (u,v,r) 的多项式（y=1，x=r u v，导数在固定 x 下取）。"""
    i, j = der
    d = defaultdict(Fr)
    for (n11, n22, n21, n12), c in acc.items():
        f = ff(n11, i)*ff(n22, j)
        if f:
            d[(n11-i+n21, n22-j+n21, n21, 0)] += c*f
    return CTX.from_dict({k: Q(c.numerator, c.denominator) for k, c in d.items() if c})


def two_point(AMP, loss):
    out = {d: ZERO for d in DERS}
    for (Mu0, Mb0), acc0 in AMP.items():
        if pc(Mu0)+pc(Mb0) != loss or (Mb0, Mu0) not in AMP:
            continue
        hops = [(0, k, 1) for k in bits(Mb0)]+[(1, k, 0) for k in bits(Mu0)]
        s = sign_of([comp((Mu0, Mb0)), comp((Mb0, Mu0))], hops)
        a1 = amp3(AMP[Mb0, Mu0])
        for d in DERS:
            out[d] += s*amp3(acc0, d)*a1
    return out


def center_pair(AMP):
    """中心 0 在两条键上各失去 4 个生成元（端点各失去 4 个）的 t^8 联合项，对中心求导。"""
    A4 = {k: v for k, v in AMP.items() if pc(k[0])+pc(k[1]) == 4}
    P4 = {k: amp3(v) for k, v in A4.items()}
    out, n = {d: ZERO for d in DERS}, 0
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
                cp = {d: amp3(acc, d) for d in DERS}
            pr = P4[e1]*P4[e2]
            for d in DERS:
                out[d] += s*cp[d]*pr
    return out, n


# ---------------------------------------------------------------- 分母显式的有理函数
class Fq:
    """n / base^k；base 是该读数的静态分母 D0（每个读数一行，行内同一 base）。"""
    def __init__(self, n, base, k):
        self.n, self.base, self.k = n, base, k

    def der(self, x):
        return Fq(self.n.derivative(x)*self.base-self.k*self.n*self.base.derivative(x), self.base, self.k+1)

    def lift(self, K):
        return self.n*self.base**(K-self.k)

    def iszero(self):
        return self.n == 0


def series_readouts(Zj):
    """Zj[d] = {0: ., 4: ., 8: .}；返回每个读数的 (X0, X4, X8)。"""
    def smul(a, b):
        return {k: sum((a[i]*b[k-i] for i in (0, 4, 8) if k-i in (0, 4, 8)), ZERO) for k in (0, 4, 8)}
    out = {}
    for key, (a, b, c, d) in PAIRS.items():
        N, D = smul(Zj[a], Zj[b]), smul(Zj[c], Zj[d])
        D0 = D[0]
        out[key] = (Fq(N[0], D0, 1), Fq(N[4]*D0-N[0]*D[4], D0, 2),
                    Fq(N[8]*D0**2-N[4]*D[4]*D0-N[0]*D[8]*D0+N[0]*D[4]**2, D0, 3))
    return out


def cross_zero(a, b):
    """a, b 为三行列向量（第 i 分量属于读数 i）；判定叉积恒为零。"""
    L = []
    for i in range(3):
        K = max(a[i].k, b[i].k)
        L.append((a[i].lift(K), b[i].lift(K)))
    return all(L[i][0]*L[j][1] == L[j][0]*L[i][1] for i, j in ((0, 1), (1, 2), (0, 2)))


def det3(c1, c2, c3, K=4):
    """返回 det 的分子；分母为 prod_i D0_i^K。"""
    m = [[c.lift(K) for c in (c1[i], c2[i], c3[i])] for i in range(3)]
    return (m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
            + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]))


def main():
    t0 = time.time()
    AMP = _g["site_amplitudes"]()   # 内部断言缺失 0,4,8 个生成元的振幅数为 1,168,6592
    s0 = {d: amp3(AMP[(0, 0)], d) for d in DERS}
    s = s0[(0, 0)]
    one = {d: {0: s0[d], 4: ZERO, 8: ZERO} for d in DERS}
    X0 = {k: series_readouts(one)[k][0] for k in KEYS}
    for k in KEYS:
        assert X0[k].der("u").iszero() and X0[k].der("v").iszero()   # 静态读数只依赖 rho
    Sp = [X0[k].der("r") for k in KEYS]
    Spp = [x.der("r") for x in Sp]
    print(f"[{time.time()-t0:.0f}s] 单点振幅与静态读数")

    # ------------------------------------------------ (1) 引理：失去 4 个生成元的每个振幅，一阶响应平行于 S'
    nlem = 0
    for key, acc in AMP.items():
        if pc(key[0])+pc(key[1]) != 4:
            continue
        a = {d: amp3(acc, d) for d in DERS}
        ser = series_readouts({d: {0: s0[d], 4: a[d], 8: ZERO} for d in DERS})
        assert cross_zero([ser[k][1] for k in KEYS], Sp)
        nlem += 1
    assert nlem == 168
    print(f"[{time.time()-t0:.0f}s] 引理：168 个失去四个生成元的振幅，一阶响应都与 S'(rho) 平行（对 u,v,rho 恒等）")

    # ------------------------------------------------ (2) 星图，z 为符号
    B4, B8 = two_point(AMP, 4), two_point(AMP, 8)
    PCt, npc = center_pair(AMP)
    assert npc == 17880
    half = Q(1, 2)
    Zj = {d: {0: s*s*s0[d], 4: Z*s*B4[d], 8: Z*s*B8[d]+(Z*Z-Z)*half*PCt[d]} for d in DERS}
    ser = series_readouts(Zj)
    Sp = [ser[k][0].der("r") for k in KEYS]     # 与星图读数同一行分母的 S'
    Spp = [x.der("r") for x in Sp]
    X4 = [ser[k][1] for k in KEYS]
    X8 = [ser[k][2] for k in KEYS]
    # t^4：X4 = z g，g 平行于 S'（对一切 rho；与成键分解一致）
    assert cross_zero(X4, Sp)
    for x in X4:
        assert x.n.subs({"z": 0}) == 0
    print(f"[{time.time()-t0:.0f}s] 星图读数到 t^8（z 为符号）")

    col = lambda xs, v: [x.der(v) for x in xs]
    r0, r4, u4, v4, u8, v8 = Sp, col(X4, "r"), col(X4, "u"), col(X4, "v"), col(X8, "u"), col(X8, "v")
    assert cross_zero(u4, Sp) and cross_zero(v4, Sp)
    J8 = det3(r0, u4, v4)
    J12 = det3(r0, u4, v8)+det3(r0, u8, v4)+det3(r4, u4, v4)
    J16n = det3(r0, u8, v8)+det3(r4, u4, v8)+det3(r4, u8, v4)
    J16d = r0[0].base**4*r0[1].base**4*r0[2].base**4
    assert J8 == 0 and J12 == 0
    J1n = J16n.subs({"z": 1})
    assert J16n-Z*Z*J1n == 0     # 配位数平方律，对一切 (u,v,rho)
    assert J16d.degrees()[3] == 0 and J1n != 0  # 分母与 z 无关；系数不恒为零
    print(f"[{time.time()-t0:.0f}s] 星图：t^8, t^12 系数为零；t^16 系数 = z^2 * (z=1 的值)，对 (u,v,rho,z) 恒等")

    # 机理：R8 的 z^2 系数 = l^2 S''/2 (mod S')，l = g_AB/S'_AB。乘以 S'_AB^2 后逐分量化为多项式恒等式。
    def z2part(p):
        return CTX.from_dict({(a, b, c, 0): co for (a, b, c, e), co in p.to_dict().items() if e == 2})
    A8 = [(z2part(x.n), x.base**3) for x in X8]                       # R8 的 z^2 系数
    G = [(x.n.subs({"z": 1}), x.base**2) for x in X4]                 # g
    P1 = [(x.n, x.base**2) for x in Sp]                               # S'
    P2 = [(x.n, x.base**3) for x in Spp]                              # S''
    sab = P1[0]
    # Et_i = S'_AB^2 A8_i - g_AB^2 S''_i / 2
    Et = []
    for i in range(3):
        n1, d1 = sab[0]**2*A8[i][0], sab[1]**2*A8[i][1]
        n2, d2 = G[0][0]**2*P2[i][0]*half, G[0][1]**2*P2[i][1]
        Et.append((n1*d2-n2*d1, d1*d2))
    for i, j in ((0, 1), (1, 2), (0, 2)):
        assert Et[i][0]*P1[j][0]*Et[j][1]*P1[i][1] == Et[j][0]*P1[i][0]*Et[i][1]*P1[j][1]
    print(f"[{time.time()-t0:.0f}s] 机理：R8 的 z^2 系数与 l^2 S''/2 之差平行于 S'")

    # ------------------------------------------------ (3) J16^(2)(1,w,rho) 的结构
    c2 = flint.fmpq_mpoly_ctx.get(("w", "r"), "lex")

    def to_wr(p):
        p = p.subs({"u": 1, "z": 1})
        return c2.from_dict({(b, c): co for (a, b, c, e), co in p.to_dict().items()})
    num, den = to_wr(J16n), to_wr(J16d)
    cn, fn = num.factor()
    cd, fd = den.factor()
    W, Rv = c2.gens()
    cub = 20*Rv**3+45*Rv**2+18*Rv+1
    pq = 15*Rv**4+80*Rv**3+90*Rv**2+24*Rv+1
    fnd, fdd = dict((str(f), m) for f, m in fn), dict((str(f), m) for f, m in fd)
    F = [(f, m) for f, m in fn if f.degrees()[0] >= 1 and len(f) > 1]
    assert len(F) == 1 and F[0][1] == 1
    F = F[0][0]
    assert F.degrees() == (6, 26) and len(F) == 181
    expn = {"r": 2, str(cub): fnd[str(cub)]-fdd[str(cub)], str(pq): fnd[str(pq)]-fdd[str(pq)], "w": fnd["w"]-fdd["w"]}
    assert fnd["r"] == 2 and "r" not in fdd
    assert set(fnd) == {"r", "w", str(cub), str(pq), str(F)} and set(fdd) == {"w", str(cub), str(pq)}
    assert expn[str(cub)] == -8 and expn[str(pq)] == -2 and expn["w"] == -8
    # 规范化多项式 P（代码中变量名 F/Fint）：整系数、本原、在 (w,rho)=(1,1/2) 处取正号；常数 C 使 J16^(2)(1,w,rho) = C rho^2 P/(w^8 q^8 p^2)
    from math import gcd, lcm
    L = 1
    for v in F.to_dict().values():
        L = lcm(L, int(v.q))
    dz = {k: int((v*L).p) for k, v in F.to_dict().items()}
    gg = 0
    for v in dz.values():
        gg = gcd(gg, v)
    dz = {k: v//gg for k, v in dz.items()}
    ev = lambda p, sub: p.subs(sub).to_dict().get((0,)*len(p.context().names()), Q(0))
    Fint = c2.from_dict({k: Q(v) for k, v in dz.items()})
    if ev(Fint, {"w": 1, "r": Q(1, 2)}) < 0:
        Fint = -Fint
        dz = {k: -v for k, v in dz.items()}
    # 由 rho=1/2, w=1 的已发表值反推 C，再在 rho=1/3 与 rho=1/2 的整条 w 曲线上独立核对
    ev = lambda p, sub: p.subs(sub).to_dict().get((0,)*len(p.context().names()), Q(0))
    shape = lambda wv, rv: (Q(rv)**2*ev(Fint, {"w": wv, "r": rv})
                            / (Q(wv)**8*ev(cub, {"w": wv, "r": rv})**8*ev(pq, {"w": wv, "r": rv})**2))
    J2 = lambda wv, rv: ev(J1n, {"u": 1, "v": wv, "r": rv, "z": 1})/ev(J16d, {"u": 1, "v": wv, "r": rv, "z": 1})
    known12 = Q(1031498737935386099, 368669844827839227494400000)
    known13 = Q(10593838074357783855, 1380913646891877604595335168)
    assert J2(1, Q(1, 2)) == known12 and J2(1, Q(1, 3)) == known13
    C = known12/shape(1, Q(1, 2))
    for wv in (Q(1), Q(2), Q(1, 3), Q(7, 5)):
        for rv in (Q(1, 2), Q(1, 3), Q(1), Q(1, 10)):
            assert J2(wv, rv) == C*shape(wv, rv)
    assert C > 0
    # 齐次性 J16^(2)(2u,2v,rho) = 2^{-10} J16^(2)(u,v,rho)，作为多项式恒等式：
    # J1n(2u,2v) * J16d(u,v) = 2^{-10} J1n(u,v) * J16d(2u,2v)
    sc = lambda p: p.compose(2*U, 2*V, RR, Z)
    assert sc(J1n)*J16d == Q(1, 1024)*J1n*sc(J16d)
    # rho=1/2 时与已发表的 P6(w) 闭式一致（作为 w 的有理函数）
    c1 = flint.fmpq_poly
    P6 = c1([35044430945354514432, -48252152743900501760, -89149185028361299968, -20828682148361300688,
             -5571824064272581248, -188484971655861335, 8555769273768192])
    Fhalf = c1([0]*7)
    coeffs = defaultdict(lambda: Q(0))
    for (a, b), v in Fint.to_dict().items():
        coeffs[a] += v*Q(1, 2)**b
    Fhalf = c1([coeffs[i] for i in range(7)])
    lhs = C*Q(1, 4)/(ev(cub, {"w": 1, "r": half})**8*ev(pq, {"w": 1, "r": half})**2)
    assert Fhalf*lhs == P6*Q(-1, 46083730603479903436800000000)
    # P 作为 w 的多项式，各系数（rho 的多项式）的最大公因式为 1：没有 rho 使 P(.,rho) 恒为零
    rho = sp.Symbol("rho")
    wcoef = defaultdict(lambda: sp.Integer(0))
    for (a, b), v in dz.items():
        wcoef[a] += v*rho**b
    gw = wcoef[0]
    for a in range(1, 7):
        gw = sp.gcd(gw, wcoef[a])
    assert sp.Poly(gw, rho).degree() == 0
    assert sp.Poly(wcoef[6], rho).degree() >= 1
    # rho=0：P(w,0) 非零，故 J16 在 rho -> 0 时恰按 rho^2 趋零
    assert any(v != 0 for (a, b), v in dz.items() if b == 0)
    # 对角线 w=1：P(1,rho) 在 (0,1] 中的零点
    F1 = sp.Poly(sum(v*rho**b for (a, b), v in dz.items()), rho)
    nroot = F1.count_roots(0, 1)
    iv = [(a, b) for (a, b), m in F1.intervals(inf=0, sup=1, eps=sp.Rational(1, 10**6))]
    assert F1.eval(0) != 0 and F1.eval(1) != 0
    print(f"[{time.time()-t0:.0f}s] J16^(2)(1,w,rho) = C rho^2 P(w,rho)/(w^8 q^8 p^2)，P 不可约、181 项、"
          f"次数 (6,26)；C = {C}；对角线 (0,1] 上 P(1,rho) 有 {nroot} 个零点 {[(str(a), str(b)) for a, b in iv]}")
    # 个别 rho 上 w>0 的零点个数（精确 Sturm 计数）
    wpos = {}
    ws = sp.Symbol("w")
    for rv in (sp.Rational(1, 3), sp.Rational(1, 2), sp.Integer(1)):
        Pw = sp.Poly(sum(v*ws**a*rv**b for (a, b), v in dz.items()), ws)
        wpos[str(rv)] = (Pw.count_roots(0, None) - (1 if Pw.eval(0) == 0 else 0),
                         [(str(a), str(b)) for (a, b), m in Pw.intervals(inf=0, eps=sp.Rational(1, 10**4))])
    print("w>0 的零点（rho=1/3,1/2,1）：", wpos)
    # 断言：常数 C；P(w,0)；对角线 P(1,rho) 在 (0,1] 内恰一个零点及其括号；rho=1/3,1/2,1 时 w>0 上恰两个零点及其括号
    assert C == Q(1, 12582912)
    assert sp.expand(sum(v*ws**a for (a, b), v in dz.items() if b == 0) + 96*ws**2*(5*ws**2-96*ws+80)) == 0
    sgn = lambda P, x: sp.sign(P.eval(sp.Rational(x)))
    assert nroot == 1 and sgn(F1, "0.92958312")*sgn(F1, "0.92958313") < 0
    brackets = {"1/3": [("0.1187164", "0.1187165"), ("134.77488", "134.77489")],
                "1/2": [("0.4000162", "0.4000163"), ("39.998373", "39.998374")],
                "1": [("1.0856553", "1.0856554"), ("14.737641", "14.737642")]}
    for rv, br in brackets.items():
        Pw = sp.Poly(sum(v*ws**a*sp.Rational(rv)**b for (a, b), v in dz.items()), ws)
        assert wpos[rv][0] == 2
        for lo, hi in br:
            assert sgn(Pw, lo)*sgn(Pw, hi) < 0

    # ------------------------------------------------ (4) 模型壳图 G19
    P, _ = _g["shell_graph"]()   # 内部断言 G19 有 135 条边、529 个三角形、度分布 10x6,14x3,16x6,18x4
    d2 = lambda p, q: sum((x-y)**2 for x, y in zip(p, q))
    deg = Counter(sum(1 for q in P if q != p and d2(p, q) <= 804) for p in P)
    assert deg == {10: 6, 14: 3, 16: 6, 18: 4}
    mult = {z: z*z for z in deg}
    assert mult == {10: 100, 14: 196, 16: 256, 18: 324}
    print(f"[{time.time()-t0:.0f}s] G19 四类格点的 t^16 系数倍数：{mult}")

    data = {"description": "P(w,rho): J16^(2)(1,w,rho) = C rho^2 P(w,rho)/(w^8 q(rho)^8 p(rho)^2), q = 1+18rho+45rho^2+20rho^3, "
                           "p = 1+24rho+90rho^2+80rho^3+15rho^4; w = |F|^2 with F the Wilson factor; keys 'a,b' = coefficient of w^a rho^b",
            "C": str(C), "P": {f"{a},{b}": str(v) for (a, b), v in sorted(dz.items())},
            "diagonal_zeros_in_(0,1]": [[str(a), str(b)] for a, b in iv]}
    out = HERE / "ch39-proton-overlap-jacobian-coordination-law.json"
    if out.exists():
        assert json.loads(out.read_text()) == data
    else:
        out.write_text(json.dumps(data, indent=1, ensure_ascii=False))
    print(f"用时 {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()

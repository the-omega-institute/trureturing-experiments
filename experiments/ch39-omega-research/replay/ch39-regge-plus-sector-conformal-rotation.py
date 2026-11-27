"""正偏振扇区的两种格点共形围道旋转良定，旋转后的协方差在所检有限盒上对对角时间切片的两类站点反射都不反射正。

裸欧氏雷吉二次型在投影动量族 (t+s,t+s,t-s,t-s) 上的 (+,+) 扇区（点式对称 (12)、(34) 下双偶）
（说明见同目录《正偏振扇区的共形围道旋转在投影动量族上良定而两种格点旋转下的协方差在所检有限盒上对对角切片两类站点反射均不反射正》）。

约定沿用 ch39-regge-vector-sector-reflection.py：场是平方边长扰动，按基点记在 (n,m) 上，
n=Σx_i，m=x1+x2-x3-x4，X=e^{it}、Y=e^{is}，平移 (Δn,Δm) 的相位为 X^Δn Y^Δm，K=-H。
扇区坐标取 8 个轨道代表 S=1,3,4,5,7,12,13,15（边型 e_S 的掩码），坐标值即轨道上每条边的取值。

两种围道旋转（模型选择）：
  迹旋转   ：w(k) 为顶点外尔模 q_S(x)=L_S(σ(x)+σ(x+e_S))，L_S=|S|(5-|S|)/5，
             K_w = K - 2 Kw w†K/(w†Kw)；
  谱旋转   ：在平方边长的平直度量（扇区坐标下 D=diag(轨道长度)）中取 |K|。
"""
from fractions import Fraction as Fr
from itertools import product
from math import comb
from pathlib import Path
from runpy import run_path
import time

import numpy as np
import sympy as sp
from flint import acb, acb_mat, arb, ctx, fmpq, fmpq_mat

ROOT = Path(__file__).parent
VEC = run_path(str(ROOT/"ch39-regge-vector-sector-reflection.py"))
Box, step, MASKS = VEC["Box"], VEC["step"], VEC["MASKS"]
sector_basis, perm_mask, inertia = VEC["sector_basis"], VEC["perm_mask"], VEC["inertia"]
colbasis, nullspace, ST = VEC["colbasis"], VEC["nullspace"], VEC["ST"]
X, Y = VEC["X"], VEC["Y"]
CHARS = (1, 1)
BASIS = sector_basis(CHARS)
REPS = [S for S, _ in BASIS]
DD = [len(o) for _, o in BASIS]                       # 轨道长度 = 平方边长度量
SITE_REFL = [(2, 3, 0, 1), (3, 2, 1, 0)]
ctx.prec = 200


# ---------------------------------------------------------------- 符号 ----
def laurent(e):
    e = sp.expand(sp.cancel(e)*X**12*Y**12)
    if e == 0:
        return {}
    return {(p-12, q-12): Fr(int(c.p), int(c.q)) for (p, q), c in sp.Poly(e, X, Y).terms()}


def symbols():
    H, G = VEC["full_symbol"]()
    Bm = sp.zeros(15, 8)
    for j, (_, o) in enumerate(BASIS):
        for k, v in o.items():
            Bm[MASKS.index(k), j] = v
    assert [S for S, _ in BASIS] == [1, 3, 4, 5, 7, 12, 13, 15]
    assert Bm.T*Bm == sp.diag(*DD) and DD == [2, 1, 2, 4, 2, 1, 2, 1]
    Ks = (-(Bm.T*H*Bm)).applyfunc(sp.cancel)
    Gs = ((Bm.T*Bm).inv()*Bm.T*G).applyfunc(sp.cancel)
    G2 = (Gs*sp.Matrix([[1, 0], [1, 0], [0, 1], [0, 1]])).applyfunc(sp.cancel)  # ξ∝(1,1,0,0),(0,0,1,1)
    assert (Gs.subs({X: sp.Rational(2), Y: sp.Rational(3)})).rank() == 2
    assert (Ks*G2).applyfunc(sp.cancel) == sp.zeros(8, 2)
    K = [[laurent(Ks[i, j]) for j in range(8)] for i in range(8)]
    Gd = [[laurent(G2[i, j]) for j in range(2)] for i in range(8)]
    W = []
    for S in REPS:
        dn, dm = step(S)
        L = Fr(bin(S).count("1")*(5-bin(S).count("1")), 5)
        W.append({(0, 0): L, (dn, dm): L} if (dn, dm) != (0, 0) else {(0, 0): 2*L})
    return K, Gd, W, Ks, G2


def conj_d(d):
    return {(-p, -q): c for (p, q), c in d.items()}


def mul_d(a, b):
    out = {}
    for (p, q), c in a.items():
        for (r, s), e in b.items():
            out[(p+r, q+s)] = out.get((p+r, q+s), 0)+c*e
    return {k: v for k, v in out.items() if v}


def add_d(a, b, sb=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0)+sb*v
    return {k: v for k, v in out.items() if v}


# 规范不变行列式 c = -3 Q / 1280000，Q = Σ a_{pq} X^p Y^q（在环面上为实）。
QCOEF = {(6, 0): 24, (5, 1): 4, (5, -1): 4, (4, 2): 4, (4, 0): -164, (4, -2): 4, (3, 3): -1,
         (3, 1): -29, (3, -1): -29, (3, -3): -1, (2, 2): -62, (2, 0): -716, (2, -2): -62,
         (1, 3): -1, (1, 1): -333, (1, -1): -333, (1, -3): -1, (0, 2): 4, (0, 0): 3376, (0, -2): 4}
QCOEF.update({(-p, -q): v for (p, q), v in list(QCOEF.items())})
CINV = {k: Fr(-3*v, 1280000) for k, v in QCOEF.items()}


def check_compound(K, Gd):
    """对三组行 {a,b}：删去 a,b 的主子式 det K_E = c·u_ab·ū_ab（u_ab 为 G 的 2×2 Plücker 子式）。
    以 Laurent 多项式恒等式验证：乘单项式化为多项式后在超过次数的整点网格上逐点为零。"""
    plucker = {}
    for a, b in [(0, 2), (0, 3), (2, 3)]:
        u = add_d(mul_d(Gd[a][0], Gd[b][1]), mul_d(Gd[a][1], Gd[b][0]), -1)
        plucker[(a, b)] = u
        idx = [i for i in range(8) if i not in (a, b)]
        rhs = mul_d(mul_d(CINV, u), conj_d(u))
        # 次数界
        lo = [sum(min(p for d in (K[i][j] for j in idx) if d for (p, _) in d) for i in idx),
              sum(min(q for d in (K[i][j] for j in idx) if d for (_, q) in d) for i in idx)]
        hi = [sum(max(p for d in (K[i][j] for j in idx) if d for (p, _) in d) for i in idx),
              sum(max(q for d in (K[i][j] for j in idx) if d for (_, q) in d) for i in idx)]
        lo = [min(lo[0], min(p for p, _ in rhs)), min(lo[1], min(q for _, q in rhs))]
        hi = [max(hi[0], max(p for p, _ in rhs)), max(hi[1], max(q for _, q in rhs))]
        Kq = {(i, j): [(p, q, fmpq(c.numerator, c.denominator)) for (p, q), c in K[i][j].items()] for i in idx for j in idx}
        Rq = [(p, q, fmpq(c.numerator, c.denominator)) for (p, q), c in rhs.items()]
        npts = 0
        for x0, y0 in product(range(1, hi[0]-lo[0]+2), range(1, hi[1]-lo[1]+2)):
            xp = {p: fmpq(x0)**p if p >= 0 else 1/fmpq(x0)**(-p) for p in range(-40, 41)}
            yp = {q: fmpq(y0)**q if q >= 0 else 1/fmpq(y0)**(-q) for q in range(-40, 41)}
            val = lambda L: sum((c*xp[p]*yp[q] for p, q, c in L), fmpq(0))
            M = fmpq_mat([[val(Kq[(i, j)]) for j in idx] for i in idx])
            assert M.det() == val(Rq)
            npts += 1
        assert npts == (hi[0]-lo[0]+1)*(hi[1]-lo[1]+1)
    # Plücker 子式的显式因式；在实环面上三者同时为零当且仅当 k≡0（X=Y=±1）。
    Xs, Ys = X, Y
    tosym = lambda d: sum(sp.Rational(c.numerator, c.denominator)*Xs**p*Ys**q for (p, q), c in d.items())
    f02, f03, f23 = (sp.factor(tosym(plucker[k])) for k in [(0, 2), (0, 3), (2, 3)])
    assert sp.simplify(f02-4*(X-Y)*(X*Y-1)/(5*Y)) == 0
    assert sp.simplify(f03-4*(X-1)*(X+1)*(X*Y-1)/5) == 0
    assert sp.simplify(f23+4*(X-1)*(X+1)*(X-Y)/(5*Y)) == 0
    return plucker


def q_in_uv():
    """Q 换到 4 维动量 u=t+s=k1=k2、v=t-s=k3=k4：Q = A + sin u sin v · B，a=1-cos u，b=1-cos v。"""
    u, v, a, b, cu, cv, su, sv = sp.symbols("u v a b cu cv su sv")
    expr = sum(c*sp.cos(sp.Rational(p+q, 2)*u+sp.Rational(p-q, 2)*v) for (p, q), c in QCOEF.items())
    assert all((p+q) % 2 == 0 for p, q in QCOEF)
    expr = sp.expand(sp.expand_trig(sp.expand(expr))).subs(
        {sp.cos(u): cu, sp.cos(v): cv, sp.sin(u): su, sp.sin(v): sv})
    pe = sp.Poly(sp.expand(expr), su, sv)
    A = Bc = 0
    for (i, j), co in pe.terms():
        assert i % 2 == j % 2
        if i % 2 == 0:
            A += co*(1-cu**2)**(i//2)*(1-cv**2)**(j//2)
        else:
            Bc += co*(1-cu**2)**((i-1)//2)*(1-cv**2)**((j-1)//2)
    A = sp.expand(sp.expand(A).subs({cu: 1-a, cv: 1-b}))
    Bc = sp.expand(sp.expand(Bc).subs({cu: 1-a, cv: 1-b}))
    F = sp.expand(A**2-(2*a-a**2)*(2*b-b**2)*Bc**2)
    assert A.subs({a: 0, b: 0}) == 0 and Bc.subs({a: 0, b: 0}) == 2400
    assert sp.Poly(A, a, b).coeff_monomial(a) == 3600 and sp.Poly(A, a, b).coeff_monomial(b) == 3600
    return A, Bc, F, (a, b)


def positivity_certificate(A, F, ab):
    """在 [0,2]^2\\{0} 上 A>0 且 F>0（有理区间分支定界；含原点的盒用二次型分解）。
    由此 Q=A+zB 在 z^2=(2a-a^2)(2b-b^2) 时恒正，即 c<0 于去心环面。"""
    a, b = ab
    terms = lambda P: {(i, j): Fr(int(c.p), int(c.q)) for (i, j), c in sp.Poly(P, a, b).terms()}

    def shift(T, ac, bc):
        out = {}
        for (i, j), c in T.items():
            for k in range(i+1):
                for l in range(j+1):
                    out[(k, l)] = out.get((k, l), 0)+c*comb(i, k)*comb(j, l)*ac**(i-k)*bc**(j-l)
        return out

    def bounds(T, a0, a1, b0, b1):
        ac, bc, ha, hb = (a0+a1)/2, (b0+b1)/2, (a1-a0)/2, (b1-b0)/2
        S = shift(T, ac, bc)
        e = sum(abs(c)*ha**k*hb**l for (k, l), c in S.items() if (k, l) != (0, 0))
        return S.get((0, 0), 0)-e, S.get((0, 0), 0)+e

    TA, TF = terms(A), terms(F)
    assert (0, 0) not in TA and all(i+j >= 2 for i, j in TF)
    alpha = {(i-1, j): c for (i, j), c in TA.items() if i >= 1}
    beta = {(0, j-1): c for (i, j), c in TA.items() if i == 0}
    pq = {(i-2, j): c for (i, j), c in TF.items() if i >= 2}
    qq = {(0, j-1): c for (i, j), c in TF.items() if i == 1}
    rq = {(0, j-2): c for (i, j), c in TF.items() if i == 0}
    boxes = [0]

    def near(a0, a1, b0, b1):   # A = aα+bβ，F = a²p+ab q+b² r，α,β,p,r>0 且 4pr>q²
        if bounds(alpha, a0, a1, b0, b1)[0] <= 0 or bounds(beta, a0, a1, b0, b1)[0] <= 0:
            return False
        pl, rl = bounds(pq, a0, a1, b0, b1)[0], bounds(rq, a0, a1, b0, b1)[0]
        ql, qu = bounds(qq, a0, a1, b0, b1)
        return pl > 0 and rl > 0 and 4*pl*rl > max(abs(ql), abs(qu))**2

    def cert(a0, a1, b0, b1, depth=0):
        boxes[0] += 1
        if a0 == 0 and b0 == 0:
            if near(a0, a1, b0, b1):
                return True
        elif bounds(TA, a0, a1, b0, b1)[0] > 0 and bounds(TF, a0, a1, b0, b1)[0] > 0:
            return True
        assert depth < 30
        am, bm = (a0+a1)/2, (b0+b1)/2
        return all(cert(*bx, depth+1) for bx in
                   [(a0, am, b0, bm), (am, a1, b0, bm), (a0, am, bm, b1), (am, a1, bm, b1)])

    assert cert(Fr(0), Fr(2), Fr(0), Fr(2))
    return boxes[0]


def leading_orders():
    """c ≈ -(9/640)(t²+s²/5)，w†Kw ≈ -24(t²+s²/5)：Laurent 系数的二阶矩。"""
    t, s = sp.symbols("t s")
    assert sum(QCOEF.values()) == 0
    q2 = sp.expand(-sum(c*(p*t+q*s)**2 for (p, q), c in QCOEF.items())/2)
    assert q2 == 6000*t**2+1200*s**2
    c2 = sp.expand(-sp.Rational(3, 1280000)*q2)
    assert c2 == -sp.Rational(9, 640)*(t**2+s**2/5)
    return c2, t, s


WKW = {(4, 0): 16, (3, 1): 6, (3, -1): 6, (2, 2): 3, (2, 0): 12, (2, -2): 3, (1, 1): 32, (1, -1): 32}
WKW.update({(-p, -q): v for (p, q), v in list(WKW.items())})
WKW[(0, 0)] = -220


def check_conformal(K, W):
    """w†Kw = (6/125)(Σ_{(p,q)≠0} c_pq X^pY^q - 220)，c_pq ≥ 0 且 Σc_pq = 220。"""
    Kw = [{} for _ in range(8)]
    for i in range(8):
        for j in range(8):
            Kw[i] = add_d(Kw[i], mul_d(K[i][j], W[j]))
    q = {}
    for i in range(8):
        q = add_d(q, mul_d(conj_d(W[i]), Kw[i]))
    assert q == {k: Fr(6*v, 125) for k, v in WKW.items()}
    assert all(v > 0 for k, v in WKW.items() if k != (0, 0)) and sum(v for k, v in WKW.items() if k != (0, 0)) == 220
    # Kw = O(k^2)：每个分量的零阶与一阶矩为零；于是旋转修正 -2Kw w†K/(w†Kw) = O(k^2)，K_w(0)=K(0)。
    for d in Kw:
        assert sum(d.values()) == 0
        assert sum(c*p for (p, qq), c in d.items()) == 0 and sum(c*qq for (p, qq), c in d.items()) == 0
    assert any(sum(c*p*p for (p, qq), c in d.items()) != 0 for d in Kw)
    # 等号要求 cos(t+s)=cos(t-s)=1（系数 32 的项），即 k≡0。
    t, s = sp.symbols("t s")
    q2 = sp.expand(-sp.Rational(3, 125)*sum(v*(p*t+qq*s)**2 for (p, qq), v in WKW.items() if (p, qq) != (0, 0)))
    assert q2 == -24*t**2-sp.Rational(24, 5)*s**2
    # 迹旋转不是有限程：|(Kw)_0|^2 的分子不能被 w†Kw 的（不可约）分子整除。
    num = sp.expand(sum(sp.Rational(v)*X**(p+4)*Y**(qq+2) for (p, qq), v in WKW.items()))
    assert len(sp.factor_list(num)[1]) == 1 and sp.factor_list(num)[1][0][1] == 1
    e0 = mul_d(Kw[0], conj_d(Kw[0]))
    sh = 20
    e0p = sp.expand(sum(sp.Rational(c.numerator, c.denominator)*X**(p+sh)*Y**(qq+sh) for (p, qq), c in e0.items()))
    rem = sp.rem(sp.Poly(e0p, X, domain=sp.QQ.frac_field(Y)), sp.Poly(num, X, domain=sp.QQ.frac_field(Y)))
    assert not rem.is_zero
    return Kw, q


def zero_momentum(Ks):
    K0 = Ks.subs({X: 1, Y: 1})
    ev = K0.eigenvals()
    pos = sum(m for v, m in ev.items() if v > 0)
    neg = sum(m for v, m in ev.items() if v < 0)
    assert (pos, neg, 8-pos-neg) == (3, 0, 5)
    return K0


def limit_lines(K, Gd):
    """沿射线 k=ε(α,β)：负本征线（D 度量）的极限由约化束 (K_eff, D_N0) 的负本征向量给出；
    三个方向 (1,0)、(0,1)、(1,1) 的极限线两两不同，所以负本征线在 k=0 处没有连续延拓。"""
    al, be = sp.symbols("alpha beta")
    K0 = sp.zeros(8); K1 = sp.zeros(8); K2 = sp.zeros(8); G1 = sp.zeros(8, 2)
    for i in range(8):
        for j in range(8):
            for (p, q), c in K[i][j].items():
                c = sp.Rational(c.numerator, c.denominator); l = p*al+q*be
                K0[i, j] += c; K1[i, j] += sp.I*c*l; K2[i, j] += -c*l**2/2
        for j in range(2):
            for (p, q), c in Gd[i][j].items():
                G1[i, j] += sp.I*sp.Rational(c.numerator, c.denominator)*(p*al+q*be)
    N0 = sp.Matrix.hstack(*K0.nullspace())
    assert N0.shape == (8, 5)
    assert (N0.T*K1*N0).applyfunc(sp.expand) == sp.zeros(5, 5)
    Keff = (N0.T*(K2-K1*K0.pinv()*K1)*N0).applyfunc(sp.expand)
    D = sp.diag(*DD); Dn = N0.T*D*N0
    lines = {}
    for dirn in [(1, 0), (0, 1), (1, 1)]:
        Ke = Keff.subs({al: dirn[0], be: dirn[1]}).applyfunc(sp.nsimplify)
        assert Ke == Ke.T.conjugate() and all(e.is_real for e in Ke)
        mu = sp.symbols("mu")
        cp = sp.Poly(sp.expand((Ke-mu*Dn).det()), mu)
        assert cp.all_coeffs()[-1] == 0 and cp.all_coeffs()[-2] == 0       # 两个规范零根
        cub = sp.Poly(sp.cancel(cp.as_expr()/mu**2), mu)
        assert cub.degree() == 3
        # 精确计数：一个负根、两个正根（实根，因 Dn 正定）。
        assert cub.eval(0) != 0 and cub.count_roots(-sp.oo, 0) == 1 and cub.count_roots(0, sp.oo) == 2
        negs = cub.intervals(sup=0)
        assert len(negs) == 1
        lo, hi = cub.refine_root(*negs[0][0], eps=sp.Rational(1, 10**30))
        assert hi <= 0
        lam = arb(sp.Rational(lo).p)/arb(sp.Rational(lo).q)
        lam = lam.union(arb(sp.Rational(hi).p)/arb(sp.Rational(hi).q))
        Ka = acb_mat([[acb(arb(sp.Rational(Ke[i, j]).p)/arb(sp.Rational(Ke[i, j]).q)) for j in range(5)] for i in range(5)])
        Da = acb_mat([[acb(int(Dn[i, j])) for j in range(5)] for i in range(5)])
        Ml = Ka-acb(lam)*Da
        col = None
        for cidx in range(5):
            v = [((-1)**(i+cidx))*acb_mat([[Ml[r, c] for c in range(5) if c != i] for r in range(5) if r != cidx]).det()
                 for i in range(5)]
            if any(abs(x) > 0 for x in v):
                col = v
                break
        assert col is not None
        vec = [sum(acb(int(N0[r, i]))*col[i] for i in range(5)) for r in range(8)]
        lines[dirn] = vec
    # 两两不同：存在 2×2 子式 v_i w_j - v_j w_i 的球不含 0。
    keys = list(lines)
    for i1 in range(3):
        for i2 in range(i1+1, 3):
            v, w = lines[keys[i1]], lines[keys[i2]]
            assert any(abs(v[i]*w[j]-v[j]*w[i]) > 0 for i in range(8) for j in range(i+1, 8))
    out = {}
    for k, v in lines.items():
        f = [complex(float(x.real.mid()), float(x.imag.mid())) for x in v]
        top = max(f, key=abs)
        out[k] = [(z/top).real for z in f]
        assert all(abs((z/top).imag) < 1e-12 for z in f)
    return out


# -------------------------------------------------------------- 有限盒 ----
class PBox(Box):
    def red(self, n, m):
        return (n % self.N, m % self.M), 1


def build(N, M, periodic):
    box = (PBox if periodic else Box)(N, M)
    cols = []
    for x in box.sites:
        for _, o in BASIS:
            cols.append({box.idx[(x, k)]: fmpq(int(v.p), int(v.q)) for k, v in o.items()})
    Vm = fmpq_mat(box.D, len(cols))
    for j, c in enumerate(cols):
        for k, v in c.items():
            Vm[k, j] = v
    Vp = (Vm.transpose()*Vm).inv()*Vm.transpose()
    K = Vm.transpose()*(-box.Q)*Vm
    Gs = Vp*box.G
    Wc = fmpq_mat(K.nrows(), len(box.sites))
    for a, x in enumerate(box.sites):
        for j, S in enumerate(REPS):
            n = bin(S).count("1"); L = fmpq(n*(5-n), 5); dn, dm = step(S)
            for base in [x, (x[0]-dn, x[1]-dm)]:
                y, sg = box.red(*base)
                Wc[box.sites.index(y)*8+j, a] += sg*L
    T = Vp*box.mshift(2)*Vm+Vp*box.mshift(-2)*Vm         # τ_2+τ_2^{-1}，本征值 2cos2s
    return dict(box=box, Vm=Vm, Vp=Vp, K=K, Gs=Gs, W=Wc, T=T)


def projector(T, c, others):
    n = T.nrows(); I = fmpq_mat(n, n)
    for i in range(n):
        I[i, i] = 1
    P, sc = fmpq_mat(I), fmpq(1)
    for c2 in others:
        P = P*(T-fmpq(c2)*I); sc *= fmpq(c)-fmpq(c2)
    return P*(1/sc)


def half(B, lo, hi):
    box, Vm = B["box"], B["Vm"]
    hcols = []
    for j in range(Vm.ncols()):
        k = next(k for k in range(box.D) if Vm[k, j] != 0); (x, S) = box.var[k]
        if x[0] >= lo and x[0]+step(S)[0] <= hi:
            hcols.append(j)
    E = fmpq_mat(Vm.ncols(), len(hcols))
    for c, j in enumerate(hcols):
        E[j, c] = 1
    Ga = colbasis(B["Gs"])
    return E*nullspace(Ga.transpose()*E)


def raw_functional(box, terms):
    r = fmpq_mat(box.D, 1)
    for x, S, cf in terms:
        y, sg = box.red(*x)
        r[box.idx[(y, S)], 0] += sg*cf
    return r


# 局域见证：f1 = ½(q_{e1}+q_{e2}-q_{e12})(0,0) + ½(q_{e1}+q_{e2})(1,1)。
F1 = [((0, 0), 1, fmpq(1, 2)), ((0, 0), 2, fmpq(1, 2)), ((0, 0), 3, fmpq(-1, 2)),
      ((1, 1), 1, fmpq(1, 2)), ((1, 1), 2, fmpq(1, 2))]
# s=π/2 块的局域见证 g（支撑在 0≤n≤4，(12)、(34) 不变）。
G_PI2 = [((0, 0), 1, 1), ((0, 0), 2, 1), ((0, 0), 5, Fr(1, 2)), ((0, 0), 6, Fr(1, 2)), ((0, 0), 7, -1),
         ((0, 0), 9, Fr(1, 2)), ((0, 0), 10, Fr(1, 2)), ((0, 0), 11, -1), ((0, 2), 4, Fr(5, 2)),
         ((0, 2), 7, Fr(1, 2)), ((0, 2), 8, Fr(5, 2)), ((0, 2), 11, Fr(1, 2)), ((0, 2), 12, -1),
         ((0, 2), 13, -1), ((0, 2), 14, -1), ((1, 1), 3, -1), ((1, 1), 5, Fr(3, 2)), ((1, 1), 6, Fr(3, 2)),
         ((1, 1), 7, Fr(-1, 2)), ((1, 1), 9, Fr(3, 2)), ((1, 1), 10, Fr(3, 2)), ((1, 1), 11, Fr(-1, 2)),
         ((2, 0), 3, 1), ((3, 3), 4, Fr(1, 2)), ((3, 3), 8, Fr(1, 2))]
G_PI2 = [(x, S, fmpq(Fr(c).numerator, Fr(c).denominator)) for x, S, c in G_PI2]


def box_block(N, M, periodic, c, others, witness):
    """横向动量块（τ_2+τ_2^{-1}=c）上的迹旋转：精确有理的反射型式惯性与见证值。"""
    B = build(N, M, periodic); box, Vm, Vp, K = B["box"], B["Vm"], B["Vp"], B["K"]
    P = projector(B["T"], c, others)
    assert P*P == P and B["T"]*P == fmpq(c)*P and P == P.transpose()
    Bc = colbasis(P); BtB = Bc.transpose()*Bc
    coord = lambda Mv: BtB.inv()*Bc.transpose()*(P*Mv)
    Kc = Bc.transpose()*K*Bc
    Gc = colbasis(coord(B["Gs"])); Wc = colbasis(coord(B["W"]))
    assert Kc*Gc == fmpq_mat(Kc.nrows(), Gc.ncols())
    nk = Kc.nrows()//8
    assert inertia(Kc) == (5*nk, nk, 2*nk) and Gc.ncols() == 2*nk and Wc.ncols() == nk
    WKW_ = Wc.transpose()*Kc*Wc
    assert inertia(WKW_) == (0, nk, 0)
    KW = Kc*Wc; Kr = Kc-2*KW*WKW_.inv()*KW.transpose()
    assert inertia(Kr) == (6*nk, 0, 2*nk)
    Cr = Bc*(Kr+Gc*Gc.transpose()).inv()*Bc.transpose()
    C0 = Bc*(Kc+Gc*Gc.transpose()).inv()*Bc.transpose()
    Mw = -(WKW_.inv())
    Wfull = Bc*Wc
    out = {"momenta": nk}
    for pi in SITE_REFL:
        Pt = box.theta(pi, (0, 0, 0, 0))
        assert Pt.transpose()*box.Q*Pt == box.Q and Pt*Pt == box.pointwise((0, 1, 2, 3))
        Ths = Vp*Pt*Vm                                      # Θ 在扇区坐标上的作用
        assert Ths*Vp == Vp*Pt and Ths.transpose()*K*Ths == K and Ths*B["T"] == B["T"]*Ths
        # Θ 保持顶点外尔模子空间，所以 K_w 与 Θ 交换。
        TW, W0 = Ths*B["W"], B["W"]
        cat = fmpq_mat([[(TW[i, j] if j < TW.ncols() else W0[i, j-TW.ncols()]) for j in range(2*W0.ncols())]
                        for i in range(W0.nrows())])
        assert cat.rank() == W0.rank() == TW.rank()
        th = lambda Fx: Vm.transpose()*(Pt.transpose()*(Vp.transpose()*Fx))
        for lo, hi in [(0, 4), (1, 3)]:
            Fb = colbasis(P.transpose()*half(B, lo, hi))
            R = th(Fb).transpose()*Cr*Fb
            assert R == R.transpose()
            out[(pi, lo, hi)] = (Fb.ncols(), inertia(R))
            if (lo, hi) == (0, 4):
                # 单方向旋转的秩一公式：C_w - C_0 = 2 W M W^T（规范不变观测量上）。
                assert th(Fb).transpose()*(Cr-C0)*Fb == 2*th(Fb).transpose()*Wfull*Mw*Wfull.transpose()*Fb
                Nn = nullspace(Wfull.transpose()*Fb)
                Rn = th(Fb*Nn).transpose()*Cr*(Fb*Nn)
                assert Rn == th(Fb*Nn).transpose()*C0*(Fb*Nn)
                out[(pi, "neutral")] = (Nn.ncols(), inertia(Rn))
        wr = raw_functional(box, witness)
        assert wr.transpose()*box.G == fmpq_mat(1, box.G.ncols())           # 全理论规范不变
        assert all(0 <= box.var[k][0][0] and box.var[k][0][0]+step(box.var[k][1])[0] <= 4
                   for k in range(box.D) if wr[k, 0] != 0)
        ws = P*(Vm.transpose()*wr)
        out[(pi, "witness")] = (th(ws).transpose()*Cr*ws)[0, 0]
        if witness is F1:
            wr1 = raw_functional(box, [((x[0]+1, x[1]+1), S, cf) for x, S, cf in F1])
            assert wr1.transpose()*box.G == fmpq_mat(1, box.G.ncols())
            assert all(1 <= box.var[k][0][0] and box.var[k][0][0]+step(box.var[k][1])[0] <= 3
                       for k in range(box.D) if wr1[k, 0] != 0)
            ws1 = P*(Vm.transpose()*wr1)
            out[(pi, "strict witness")] = (th(ws1).transpose()*Cr*ws1)[0, 0]
    # PT 型（e，n(w)=0，半空间 [0,4]）：反周期 4|M 时符号差平衡（推论的实例）。
    if not periodic and M % 4 == 0:
        Pe = box.theta((0, 1, 2, 3), (0, 0, 0, 0))
        tau = box.mshift(M//2)
        assert Pe.transpose()*box.Q*Pe == box.Q and Pe*tau*Pe == -tau
        th = lambda Fx: Vm.transpose()*(Pe.transpose()*(Vp.transpose()*Fx))
        Fb = colbasis(P.transpose()*half(B, 0, 4))
        R = th(Fb).transpose()*Cr*Fb
        p_, q_, z_ = inertia(R)
        out["PT"] = (p_, q_, z_)
    return B, P, out


# -------------------------------------------------------------- Arb 谱旋转 ----
def acb_from(d, jt, qt, ls, qs):
    z = acb(0)
    for (p, q), c in d.items():
        z += acb(fmpq(c.numerator, c.denominator))*acb.exp_pi_i(acb(fmpq(p*jt, qt)+fmpq(q*ls, qs)))
    return z


def herm(M):
    return acb_mat([[M[j, i].conjugate() for j in range(M.nrows())] for i in range(M.ncols())])


def chat(K, Gd, W, k, rot, want_v=False):
    """单个动量处旋转后协方差（规范不变观测量上）。rot ∈ {trace, spec}。"""
    Km = acb_mat([[acb_from(K[i][j], *k) for j in range(8)] for i in range(8)])
    Gm = acb_mat([[acb_from(Gd[i][j], *k) for j in range(2)] for i in range(8)])
    if rot == "trace":
        w = acb_mat([[acb_from(W[i], *k)] for i in range(8)])
        Kw = Km*w; q = (herm(w)*Kw)[0, 0].real
        assert q < 0
        Kr = Km-acb_mat([[2*Kw[i, 0]*Kw[j, 0].conjugate()/q for j in range(8)] for i in range(8)])
        return (Kr+Gm*herm(Gm)).inv()
    dh = [arb(d).sqrt() for d in DD]
    Kt = acb_mat([[Km[i, j]/(dh[i]*dh[j]) for j in range(8)] for i in range(8)])
    Gt = acb_mat([[Gm[i, j]*dh[i] for j in range(2)] for i in range(8)])
    Mm = Kt+Gt*herm(Gt)
    # Mm 恰有一个负本征值（K 的惯性 (5,1,2)，Gt Gt† 在核上正定）。
    Mn = np.array([[complex(float(Mm[i, j].real.mid()), float(Mm[i, j].imag.mid())) for j in range(8)] for i in range(8)])
    e, U = np.linalg.eigh(Mn)
    assert (e < 0).sum() == 1
    I8 = acb_mat([[int(i == j) for j in range(8)] for i in range(8)])
    p = lambda x: (Mm-acb(x)*I8).det().real
    lo, hi = arb(e[0])-arb(1e-9), arb(e[0])+arb(1e-9)
    assert hi < 0 and ((p(lo) < 0 and p(hi) > 0) or (p(lo) > 0 and p(hi) < 0))   # 唯一负根在其中
    lam = arb(lo.mid()).union(arb(hi.mid()))
    Ml = Mm-acb(lam)*I8
    col = int(np.argmax(np.abs(U[:, 0])))
    v = [((-1)**(i+col))*acb_mat([[Ml[r, c] for c in range(8) if c != i] for r in range(8) if r != col]).det()
         for i in range(8)]
    vv = sum((x*x.conjugate()).real for x in v)
    assert vv > 0
    if want_v:
        return [x/d for x, d in zip(v, dh)], Km, Gm     # D 度量负本征向量（扇区坐标）
    Minv = Mm.inv()
    Ct = acb_mat([[Minv[i, j]+2*v[i]*v[j].conjugate()/((-lam)*vv) for j in range(8)] for i in range(8)])
    return acb_mat([[Ct[i, j]/(dh[i]*dh[j]) for j in range(8)] for i in range(8)])


def rank_two(K, Gd, W, k):
    """w、v 与两个规范向量线性无关（某个 4×4 子式的球不含 0），
    所以两种旋转之差 2ww†/|w†Kw|-2vv†/|v†Kv| 在规范不变观测量上的秩恰为 2。"""
    v, Km, Gm = chat(K, Gd, W, k, "spec", want_v=True)
    w = [acb_from(W[i], *k) for i in range(8)]
    cols = [w, v, [Gm[i, 0] for i in range(8)], [Gm[i, 1] for i in range(8)]]
    from itertools import combinations
    for rows in combinations(range(8), 4):
        d = acb_mat([[cols[j][r] for j in range(4)] for r in rows]).det()
        if abs(d) > 0:
            return rows
    raise AssertionError("rank-two minor not certified")


def arb_reflection(B, P, terms, momenta, K, Gd, W, rot, pi=(2, 3, 0, 1)):
    """R(f,f)=<Θf, C f>，C 由盒上块内动量逐个组装：C(x,y)=(1/|Λ|)Σ_k e^{ik(x-y)} Ĉ(k)。"""
    box, Vm, Vp = B["box"], B["Vm"], B["Vp"]
    wr = raw_functional(box, terms)
    fa = P*(Vm.transpose()*wr); raw = Vp.transpose()*fa
    Pt = box.theta(pi, (0, 0, 0, 0))
    tha = Vm.transpose()*(Pt.transpose()*raw)
    nz = lambda v: [(i//8, i % 8, v[i, 0]) for i in range(v.nrows()) if v[i, 0] != 0]
    A, Bv = nz(tha), nz(fa)
    tot = acb(0)
    for k in momenta:
        C = chat(K, Gd, W, k, rot); jt, qt, ls, qs = k
        for xi, a, ca in A:
            x = box.sites[xi]
            for yi, b, cb in Bv:
                y = box.sites[yi]
                ph = acb.exp_pi_i(acb(fmpq(jt*(x[0]-y[0]), qt)+fmpq(ls*(x[1]-y[1]), qs)))
                tot += acb(ca)*acb(cb)*ph*C[a, b]
    return tot/len(box.sites)


# ------------------------------------------------------ 无限体积（仅数值） ----
def numeric_infinite(K, Gd, W):
    def ev(d, t, s):
        return sum(float(c)*np.exp(1j*(p*t+q*s)) for (p, q), c in d.items())
    Bm = np.zeros((15, 8))
    for j, (_, o) in enumerate(BASIS):
        for k_, v in o.items():
            Bm[MASKS.index(k_), j] = float(v)
    D = np.diag(DD).astype(float)

    def C_at(t, s, rot):
        k = np.array([[ev(K[i][j], t, s) for j in range(8)] for i in range(8)])
        g = np.array([[ev(Gd[i][j], t, s) for j in range(2)] for i in range(8)])
        if rot == "trace":
            w = np.array([ev(W[i], t, s) for i in range(8)]); Kw = k@w
            k = k-2*np.outer(Kw, Kw.conj())/np.real(w.conj()@Kw)
        else:
            dh = np.sqrt(np.diag(D)); e, U = np.linalg.eigh(k/dh[:, None]/dh[None, :])
            k = (U*np.abs(e))@U.conj().T*dh[:, None]*dh[None, :]
        return Bm@np.linalg.inv(k+g@g.conj().T)@Bm.T

    def theta_terms(terms, pi):
        out = []
        for (n, m), S, c in terms:
            dn, dm = step(S)
            out.append(((-n-dn, m+dm), perm_mask(pi, S), c))
        return out

    def Rf(terms_list, s, rot, Nt=384, pi=(2, 3, 0, 1)):
        ths = [theta_terms(f, pi) for f in terms_list]
        hat = lambda terms, t: np.array([sum(float(c)*np.exp(1j*(s*m+t*n)) for (n, m), S2, c in terms if S2 == S)
                                         for S in MASKS])
        Mx = np.zeros((len(terms_list),)*2, complex)
        for t in 2*np.pi*(np.arange(Nt)+0.5)/Nt:
            C = C_at(t, s, rot)
            A = np.array([hat(th, t) for th in ths]); Bv = np.array([hat(f, t) for f in terms_list])
            Mx += A@C@Bv.conj().T
        return Mx/Nt

    f1 = [((n, m), S, float(Fr(int(c.p), int(c.q)))) for (n, m), S, c in F1]
    vals = {}
    for s in [0.1, 0.3, 0.6, 0.9]:
        for rot in ["trace", "spec"]:
            v = Rf([f1], s, rot)[0, 0]
            assert abs(v.imag) < 1e-9 and v.real < 0
            vals[(s, rot)] = v.real
    for rot in ["trace", "spec"]:
        v = Rf([f1], np.pi/2, rot)[0, 0]
        assert abs(v.imag) < 1e-9 and v.real > 0     # s=π/2 处 f1 不再是见证
        vals[("pi/2", rot)] = v.real
    return vals


def verify():
    t0 = time.time()
    K, Gd, W, Ks, G2 = symbols()
    plucker = check_compound(K, Gd)
    A, Bc, F, ab = q_in_uv()
    nbox = positivity_certificate(A, F, ab)
    assert nbox == 17 and len(QCOEF) == 37 and len(G_PI2) == 25
    for wit in (F1, G_PI2):   # 在 <(12),(34)> 下不变：只看到 (+,+) 扇区
        dct = {}
        for x, S, cf in wit:
            dct[(x, S)] = dct.get((x, S), 0)+cf
        for g in [(1, 0, 2, 3), (0, 1, 3, 2)]:
            assert {(x, perm_mask(g, S)): v for (x, S), v in dct.items()} == dct
    c2, t, s = leading_orders()
    Kw, q = check_conformal(K, W)
    K0 = zero_momentum(Ks)
    lines = limit_lines(K, Gd)
    t1 = time.time()
    print(f"[symbol] 规范不变行列式 c=-3Q/1280000，三组主子式恒等式成立；Q>0 于去心环面（{nbox} 个有理盒）")
    print("[symbol] c ≈ -(9/640)(t²+s²/5)，w†Kw ≈ -24(t²+s²/5)；w†Kw=(6/125)(Σc_pq cos-220)<0 于去心环面")
    print("[symbol] 零动量惯性 (3,0,5)；负本征线沿 (1,0),(0,1),(1,1) 的极限两两不同：", {k: [round(x, 3) for x in v] for k, v in lines.items()})
    print(f"  用时 {t1-t0:.1f} 秒")

    blocks = {
        "pi/4": (8, 4, False, 0, [], F1),
        "pi/6": (8, 6, False, 1, [-2], F1),
        "pi/2": (8, 6, False, -2, [1], G_PI2),
        "pi/3": (8, 6, True, -1, [2], F1),
    }
    expect = {
        "pi/4": dict(m=16, R=(14, 10, 0), strict=(2, 2, 0), neu=(14, (8, 6, 0)), w=fmpq(-270721840955, 3088612034704)),
        "pi/6": dict(m=16, R=(16, 8, 0), strict=(2, 2, 0), neu=(14, (8, 6, 0)), w=fmpq(-30700157525965, 155183768488104)),
        "pi/2": dict(m=8, R=(7, 5, 0), strict=(1, 1, 0), neu=(7, (4, 3, 0)), w=fmpq(-8016671, 7114338)),
        "pi/3": dict(m=16, R=(14, 10, 0), strict=(2, 2, 0), neu=(14, (8, 6, 0)), w=fmpq(-182051589335, 15708820029048)),
    }
    res, built = {}, {}
    for name, (N, M, per, c, oth, wit) in blocks.items():
        B, P, out = box_block(N, M, per, c, oth, wit)
        built[name] = (B, P, wit)
        e = expect[name]
        assert out["momenta"] == e["m"]
        for pi in SITE_REFL:
            assert out[(pi, 0, 4)][1] == e["R"], (name, out[(pi, 0, 4)])
            assert out[(pi, 1, 3)][1] == e["strict"], (name, out[(pi, 1, 3)])
            assert out[(pi, "neutral")] == e["neu"], (name, out[(pi, "neutral")])
            wv = out[(pi, "witness")]
            assert wv < 0
            assert wv == e["w"]
            if name == "pi/4":
                assert out[(pi, "strict witness")] == fmpq(-10214340412033487, 101933462981336112)
        res[name] = out
    assert res["pi/2"][(SITE_REFL[0], "witness")] == res["pi/2"][(SITE_REFL[1], "witness")]
    assert res["pi/4"]["PT"] == (12, 12, 0)
    _, _, o2 = box_block(8, 6, False, -2, [1], F1)         # s=π/2 块上 f1 不是见证
    assert all(o2[(pi, "witness")] == fmpq(2461, 404976) for pi in SITE_REFL)
    t2 = time.time()
    print("[box, 迹旋转, 精确] 各横向动量块：站点反射 [0,4] 惯性、严格半空间 [1,3] 惯性、中性观测量惯性、见证值")
    for name, out in res.items():
        pi = SITE_REFL[0]
        print("  s=", name, out[(pi, 0, 4)], out[(pi, 1, 3)], out[(pi, "neutral")], out[(pi, "witness")])
    print("  8×4 盒严格半空间见证（f1 平移到 n=1）：", res["pi/4"][(SITE_REFL[0], "strict witness")])
    print("  s=π/2 见证值 g：", res["pi/2"][(SITE_REFL[0], "witness")], "；8×4 盒 PT 型符号差", res["pi/4"]["PT"])
    print(f"  用时 {t2-t1:.1f} 秒")

    moms = {
        "pi/4": [(2*j, 8, 2*l+1, 4) for j in range(8) for l in range(2)],
        "pi/6": [(2*j, 8, l, 6) for j in range(8) for l in (1, 5)],
        "pi/2": [(2*j, 8, 3, 6) for j in range(8)],
        "pi/3": [(2*j, 8, l, 3) for j in range(8) for l in (1, 2)],
    }
    arbres, arbres2 = {}, {}
    for name, (B, P, wit) in built.items():
        vt = arb_reflection(B, P, wit, moms[name], K, Gd, W, "trace")
        exact = res[name][(SITE_REFL[0], "witness")]
        assert vt.real.overlaps(arb(exact.p)/arb(exact.q)) and vt.real.rad() < arb(1e-40)
        assert abs(vt.imag) < arb(1e-30)
        vs = arb_reflection(B, P, wit, moms[name], K, Gd, W, "spec", pi=SITE_REFL[0])
        vs2 = arb_reflection(B, P, wit, moms[name], K, Gd, W, "spec", pi=SITE_REFL[1])
        assert vs.real < 0 and abs(vs.imag) < arb(1e-4) and vs2.real.overlaps(vs.real) and vs2.real < 0
        arbres[name] = vs.real
        arbres2[name] = vs2.real
    # 单侧界（Arb 球的上端点）。
    for name, ub in {"pi/4": -0.137, "pi/6": -0.286, "pi/2": -1.90, "pi/3": -0.0183}.items():
        assert arbres[name] < arb(ub) and arbres2[name] < arb(ub)
    rows = rank_two(K, Gd, W, (2, 8, 1, 4))
    print("[box, 谱旋转（平方边长平直度量）, Arb 200 比特] 同一见证的反射期望（两个站点型反射）：",
          {k: (v.str(5, radius=True), arbres2[k].str(5, radius=True)) for k, v in arbres.items()})
    print("  两种旋转之差的秩为 2：(t,s)=(π/4,π/4) 处 [w,v,G] 的子式行", rows)
    t3 = time.time()
    print(f"  用时 {t3-t2:.1f} 秒")

    vals = numeric_infinite(K, Gd, W)
    assert vals[(0.1, "trace")] < -1.6 and vals[(0.1, "spec")] < -1.6
    assert vals[(0.9, "trace")] < -0.03 and vals[(0.9, "spec")] < -0.025
    print("[无限体积, 仅数值] f1 的反射期望：", {k: round(v, 4) for k, v in vals.items()})
    print(f"总用时 {time.time()-t0:.1f} 秒")


if __name__ == "__main__":
    verify()

"""辅助电流重叠 c=g_L^dag g_D：Q-L 扇区接触项的静态单点 Berezin 不变量的精确多项式，
以及与 kappa、|F|、双线性归一都无关的 Delta B=0 读数比 X(rho)。
静态不变量是人为定义的读数（零跳跃、无其余单点项、只对 Q,L,Qbar,Lbar 积分），不是模型配分函数的极限。

H = sum_a [K11 Dbar_a D_a + K22 (QL)_a T_a + K21 (QL)_a D_a + K12 Dbar_a T_a]，
K = G^dag G，G = [g_D, conj(F) g_L]。Z = Berezin 顶阶系数 of exp(H)（32 个 Grassmann 生成元）。
全部为精确整数/有理运算，不经过玻色化；末尾另用 det(I-KA)^{-3} 的生成函数独立复核。
"""
from collections import defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial

import sympy as sp

EPS2 = {(0, 1): 1, (1, 0): -1}
FULL = (1 << 16) - 1


def perm_sign(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def q_index(a, i, al):
    return 4*a+2*i+al


def l_index(j, be):
    return 12+2*j+be


def word(idx):
    """有序生成元乘积 -> (位掩码, 符号)；重复生成元给 None。"""
    m, s = 0, 0
    for x in idx:
        if m >> x & 1:
            return None
        s += bin(m >> (x+1)).count("1")
        m |= 1 << x
    return m, (-1)**s


def mul(A, B):
    out = defaultdict(int)
    for m1, c1 in A.items():
        for m2, c2 in B.items():
            if m1 & m2:
                continue
            s, mm = 0, m2
            while mm:
                y = (mm & -mm).bit_length()-1
                mm &= mm-1
                s += bin(m1 >> (y+1)).count("1")
            out[m1 | m2] += c1*c2*(-1)**s
    return {k: v for k, v in out.items() if v}


def bilinears(barred):
    """十六个生成元上的两种色(反)三重态双线性。未加横线：D_a=eps QQ, (QL)_a；
    加横线（另一份十六生成元）：Dbar_a 为 D_a 的反序共轭，T_a=((QL)_a)^dag 的反序。"""
    first = [defaultdict(int) for _ in range(3)]
    second = [defaultdict(int) for _ in range(3)]
    for a in range(3):
        for b, c in permutations(range(3), 2):
            if a in (b, c):
                continue
            e = perm_sign((a, b, c))
            for (i, j), ei in EPS2.items():
                for (al, be), ea in EPS2.items():
                    pair = [q_index(b, i, al), q_index(c, j, be)]
                    m, s = word(pair[::-1] if barred else pair)
                    first[a][m] += e*ei*ea*s
        for (i, j), ei in EPS2.items():
            for (al, be), ea in EPS2.items():
                pair = [q_index(a, i, al), l_index(j, be)]
                m, s = word(pair[::-1] if barred else pair)
                second[a][m] += ei*ea*s
    clean = lambda xs: [{k: v for k, v in d.items() if v} for d in xs]
    return clean(first), clean(second)


def top_table(first, second):
    pf = [_pow(B) for B in first]
    ps = [_pow(B) for B in second]
    table = {}
    for d in product(range(5), repeat=3):
        for l in product(range(5), repeat=3):
            if sum(d) != 4 or sum(l) != 4:
                continue
            r = {0: 1}
            for a in range(3):
                r = mul(r, pf[a][d[a]])
            for a in range(3):
                r = mul(r, ps[a][l[a]])
            table[d, l] = r.get(FULL, 0)
    return table


def _pow(B):
    out = [{0: 1}]
    for _ in range(4):
        out.append(mul(out[-1], B))
    return out


def full_partition(top_u, top_b):
    """Z 作为 K 四个元素的多项式：逐色枚举四种四次项的重数，不预设任何选择律。"""
    per_color = [t for t in product(range(5), repeat=4)]
    Z = defaultdict(Fraction)

    def rec(a, acc, du, lu, db, lb):
        if a == 3:
            if (sum(du), sum(lu), sum(db), sum(lb)) != (4, 4, 4, 4):
                return
            u, b = top_u.get((du, lu), 0), top_b.get((db, lb), 0)
            if u == 0 or b == 0:
                return
            w = Fraction(u*b)
            e = [0, 0, 0, 0]
            for n11, n22, n21, n12 in acc:
                w /= factorial(n11)*factorial(n22)*factorial(n21)*factorial(n12)
                e = [e[0]+n11, e[1]+n12, e[2]+n21, e[3]+n22]
            Z[tuple(e)] += w
            return
        for n11, n22, n21, n12 in per_color:
            nd, nl, ndb, nlb = n11+n21, n22+n21, n11+n12, n22+n12
            if (sum(du)+nd > 4 or sum(lu)+nl > 4 or sum(db)+ndb > 4 or sum(lb)+nlb > 4):
                continue
            rec(a+1, acc+[(n11, n22, n21, n12)], du+(nd,), lu+(nl,), db+(ndb,), lb+(nlb,))

    rec(0, [], (), (), (), ())
    return {k: v for k, v in Z.items() if v}


def verify():
    D, QL = bilinears(False)
    Db, T = bilinears(True)
    assert [len(x) for x in D+QL+Db+T] == [4]*12
    top_u, top_b = top_table(D, QL), top_table(Db, T)
    # 单侧顶阶：仅当每色 D 数等于 QL 数时非零，且精确等于 384*prod(d_a!)。
    for table in (top_u, top_b):
        for (d, l), v in table.items():
            expected = 384*factorial(d[0])*factorial(d[1])*factorial(d[2]) if d == l else 0
            assert v == expected, (d, l, v)
    Z = full_partition(top_u, top_b)
    # 指数 (e11,e12,e21,e22)。每项四种双线性各恰出现 4 次 => 对 D,Dbar,QL,T 任意重标度只乘整体常数。
    for e11, e12, e21, e22 in Z:
        assert e11+e21 == 4 and e22+e21 == 4 and e11+e12 == 4 and e22+e12 == 4
        assert e12 == e21 and e11 == e22  # 只经 K12*K21 依赖 c
    poly = {e12: Z[(4-e12, e12, e12, 4-e12)] for e12 in range(5)}
    scale = poly[0]
    assert scale == 2211840
    ratio = [poly[p]/scale for p in range(5)]
    assert ratio == [1, 24, 90, 80, 15]
    rho = sp.symbols("rho")
    p = sum(sp.Integer(int(r))*rho**k for k, r in enumerate(ratio))
    # Bernstein 形式：系数 C(6+j,6)，严格递增 => p 在 [0,1] 上严格递增；p(1)=210=C(10,6)。
    bern = sum(comb(4, j)*comb(6+j, 6)*rho**j*(1-rho)**(4-j) for j in range(5))
    assert sp.expand(bern-p) == 0
    assert [comb(6+j, 6) for j in range(5)] == [1, 7, 28, 84, 210]
    assert p.subs(rho, 1) == 210 and p.subs(rho, 0) == 1
    assert all(c > 0 for c in sp.Poly(sp.diff(p, rho), rho).all_coeffs())
    # 独立复核：玻色化生成函数 [s^4 r^4] det(I-KA)^{-3}，A=[[0,s],[r,0]]。
    k11, k12, k21, k22, s, r = sp.symbols("k11 k12 k21 k22 s r")
    det = 1-k12*r-k21*s+(k12*k21-k11*k22)*r*s
    gen = sp.expand(sp.series(sp.series(det**-3, s, 0, 5).removeO(), r, 0, 5).removeO())
    coef = gen.coeff(s, 4).coeff(r, 4)
    Zsym = sum(sp.Rational(v.numerator, v.denominator)*k11**e[0]*k12**e[1]*k21**e[2]*k22**e[3]
               for e, v in Z.items())
    assert sp.simplify(Zsym/coef) == 384**2
    # 逐色恒等式：sum_{d,e} C(d,m)C(e,m) w^m x^d y^e = 1/[(1-x)(1-y)-w x y]，截到 d,e<=4 逐项核对。
    x, y, w = sp.symbols("x y w")
    lhs = sum(comb(d, m)*comb(e, m)*w**m*x**d*y**e
              for d in range(5) for e in range(5) for m in range(min(d, e)+1))
    rhs = sp.expand(sp.series(sp.series(1/((1-x)*(1-y)-w*x*y), x, 0, 5).removeO(), y, 0, 5).removeO())
    rhs = sum(rhs.coeff(x, d).coeff(y, e)*x**d*y**e for d in range(5) for e in range(5))
    assert sp.expand(lhs-rhs) == 0
    closed = sum(sp.binomial(10-k, 2)*sp.factorial(8-k)/(sp.factorial(4-k)**2*sp.factorial(k))
                 * rho**(4-k)*(1-rho)**k for k in range(5))
    assert sp.expand(closed-15*p) == 0
    # Gram 矩阵：j_a=g_D D_a + conj(F) g_L T_a，j^dag j=Y^dag G^dag G Y。K12=conj(F c)，K21=F c。
    F, c = sp.Rational(3, 7)+sp.I/5, sp.Rational(2, 5)-sp.I*sp.Rational(1, 3)
    gD = sp.Matrix([1, 0])
    gL = sp.Matrix([sp.conjugate(c), sp.sqrt(1-c*sp.conjugate(c))])  # 任意 |c|<=1 由二维辅助空间实现
    G = sp.Matrix.hstack(gD, sp.conjugate(F)*gL)
    K = sp.simplify(G.H*G)
    assert sp.simplify((gL.H*gD)[0]-c) == 0
    assert sp.simplify(K[1, 0]-F*c) == 0 and sp.simplify(K[0, 1]-sp.conjugate(F*c)) == 0
    assert sp.simplify(K[0, 0]-1) == 0 and sp.simplify(K[1, 1]-F*sp.conjugate(F)) == 0
    tau = sp.simplify(K[0, 1]*K[1, 0]/(K[0, 0]*K[1, 1]))
    assert sp.simplify(tau-c*sp.conjugate(c)) == 0
    assert sp.simplify(K.det()-F*sp.conjugate(F)*(1-c*sp.conjugate(c))) == 0
    # Cauchy-Schwarz：|c|^2 = |gL^dag gD|^2 <= |gL|^2|gD|^2 = 1，差为 Gram 行列式（Lagrange 恒等式，符号向量）。
    a1, a2, b1, b2 = sp.symbols("a1 a2 b1 b2")
    A1, A2, B1, B2 = sp.symbols("A1 A2 B1 B2")  # 共轭分量作独立符号
    inner = A1*b1+A2*b2  # conj(gL).gD
    lagrange = (a1*A1+a2*A2)*(b1*B1+b2*B2)-inner*(a1*B1+a2*B2)
    assert sp.expand(lagrange-(a1*b2-a2*b1)*(A1*B2-A2*B1)) == 0
    assert sp.simplify(1-c*sp.conjugate(c)-(gL.H*gL)[0]*(gD.H*gD)[0]+(gL.H*gD)[0]*(gD.H*gL)[0]) == 0
    # 一维辅助空间：gD=e^{i alpha}, gL=e^{i beta}，|c|=1。
    al, be = sp.symbols("alpha beta", real=True)
    c1 = sp.exp(-sp.I*be)*sp.exp(sp.I*al)
    assert sp.simplify(c1*sp.conjugate(c1)) == 1
    # 关于 K12 一次的 Delta B=+1 读数 <sum_a (QL)^a D_a> = kappa^{-1} d log Z/dK21 = K12 p'(tau)/(K11 K22 p(tau))：仅在 c=0 为零。
    Kp = {k11: 1, k22: 1}
    logZ = sp.log(Zsym.subs(Kp))
    o6 = sp.simplify(sp.diff(logZ, k21))
    assert sp.simplify(o6 - k12*sp.diff(p, rho).subs(rho, k12*k21)/p.subs(rho, k12*k21)) == 0
    assert o6.subs({k12: 0, k21: 0}) == 0
    assert sp.diff(p, rho).subs(rho, 1)/p.subs(rho, 1) == sp.Rational(12, 5)
    # 与 kappa、|F| 及双线性归一都无关的 Delta B=0 读数比
    #   X = <A B>/(<A><B>)，A=sum_a Dbar_a D_a，B=sum_a (QL)^a T_a，
    # <A>=kappa^{-1} d_{K11} log Z，<B>=kappa^{-1} d_{K22} log Z，<AB>=kappa^{-2} d_{K11}d_{K22} Z / Z。
    lam = sp.symbols("lam", positive=True)
    zf = Zsym
    X = sp.cancel(sp.diff(zf, k11, k22)*zf/(sp.diff(zf, k11)*sp.diff(zf, k22)))
    tsym = sp.symbols("t")
    Xt = sp.cancel(X.subs({k11: 1, k22: 1, k21: 1, k12: tsym}))
    # 只经 tau 依赖：对四个 K 元任意独立重标度后，X 仍只取决于 tau。
    for sc in [{k11: lam*k11, k12: lam*k12}, {k22: lam*k22, k21: lam*k21},
               {k11: lam*k11, k21: lam*k21}, {k22: lam*k22, k12: lam*k12}]:
        assert sp.simplify(X.subs(sc, simultaneous=True)-X) == 0
    assert sp.simplify(X-Xt.subs(tsym, k12*k21/(k11*k22))) == 0
    pt = p.subs(rho, tsym)
    closedX = pt*(16*pt-7*tsym*sp.diff(pt, tsym)+tsym**2*sp.diff(pt, tsym, 2))/(4*pt-tsym*sp.diff(pt, tsym))**2
    assert sp.simplify(Xt-closedX) == 0
    assert sp.expand(4*pt-tsym*sp.diff(pt, tsym)-4*(20*tsym**3+45*tsym**2+18*tsym+1)) == 0
    assert Xt.subs(tsym, 0) == 1 and Xt.subs(tsym, 1) == sp.Rational(5, 4)
    cub = 20*tsym**3+45*tsym**2+18*tsym+1
    assert sp.simplify(Xt-(10*tsym**3+45*tsym**2+27*tsym+2)*pt/(2*cub**2)) == 0
    num, den = sp.fraction(sp.cancel(sp.diff(Xt, tsym)))
    nc = sp.Poly(num, tsym).all_coeffs()
    # 约分后：X'(t) = N(t)/[2 (20t^3+45t^2+18t+1)^3]，N 系数全正 => X 在 [0,1] 上严格增。
    assert nc == [3000, 20250, 48150, 68775, 52920, 22365, 6090, 675, 36, 3]
    assert sp.expand(den-2*cub**3) == 0
    assert all(v > 0 for v in nc)
    nX, dX = sp.fraction(Xt)
    print("归一无关读数 X(rho) =", nX, "/", dX, "；X(0)=1, X(1)=", Xt.subs(tsym, 1))
    print("X' 分子系数", nc[:4], "... 全正")
    print("一代单点顶阶系数 384*prod(d!)；Z/(K11K22)^4 =", scale, "* (1+24t+90t^2+80t^3+15t^4)")
    print("Bernstein 系数", [comb(6+j, 6) for j in range(5)], "；Z(平行)/Z(正交) =", p.subs(rho, 1))
    print("生成函数比例常数", 384**2, "；|c|^2 首阶相对系数 24")


if __name__ == "__main__":
    verify()

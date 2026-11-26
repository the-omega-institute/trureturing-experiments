"""螺旋度正负一向量扇区的格点有隙模破坏对角切片站点反射正性，正偏振所在扇区因含共形模不定。

裸雷吉作用在投影动量族 (t+s,t+s,t-s,t-s) 上的全部偏振扇区与对角时间 n=Σx_i 的反射正性证书
（说明见同目录《螺旋度正负一向量扇区的格点有隙模破坏对角切片站点反射正性而正偏振所在扇区因含共形模不定》）。

格：平直四维坐标链剖分，度量 g=I-J/5；场为平方边长扰动，只依赖 (p,q)=(x1+x2,x3+x4)，
即 n=p+q（对角时间）、m=p-q（横向），基点标号，场点构成棋盘 n≡m (mod 2)。
X=e^{it}、Y=e^{is}，平移 (Δn,Δm) 的傅里叶相位为 X^Δn Y^Δm。
扇区按点式对称 (12)、(34) 的特征 (ε12,ε34) 划分：(-,-) 交叉，(+,-)、(-,+) 向量，(+,+) 含正偏振。
"""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
from runpy import run_path
import time

import sympy as sp
from flint import fmpq, fmpq_mat


ROOT = Path(__file__).parent
FLAT = run_path(str(ROOT/"ch39-regge-flat-cell.py"))
CROSS = run_path(str(ROOT/"ch39-regge-transverse-polarization.py"))
PREV = run_path(str(ROOT/"ch39-cross-branch-reflection-classification.py"))
charpoly_psd = PREV["charpoly_psd"]
EDGES = FLAT["EDGES"]
MASKS = list(range(1, 16))
X, Y, T, V = sp.symbols("X Y T V")


def bits(S):
    return tuple((S >> k) & 1 for k in range(4))


def step(S):
    b = bits(S)
    return sum(b), b[0]+b[1]-b[2]-b[3]


def perm_mask(pi, S):
    return sum(1 << pi[k] for k in range(4) if (S >> k) & 1)


def stencil():
    """裸雷吉二次型在 (Δn,Δm) 上的实空间模板：(S,S') -> {(Δn,Δm): 系数}（基点差）。"""
    local = sum((m for _, _, m in FLAT["canonical"]()), sp.zeros(10))
    out = defaultdict(lambda: defaultdict(lambda: sp.Integer(0)))
    for order in permutations(range(4)):
        vs = [sum(1 << i for i in order[:j]) for j in range(5)]
        for (i, j), val in local.todok().items():
            a, b = [vs[k] for k in EDGES[i]]
            c, d = [vs[k] for k in EDGES[j]]
            da = [x-y for x, y in zip(bits(c), bits(a))]
            out[(b-a, d-c)][(sum(da), da[0]+da[1]-da[2]-da[3])] += val
    return out


ST = stencil()


def conj(e):
    return e.subs({X: 1/X, Y: 1/Y}, simultaneous=True)


def full_symbol():
    H = sp.zeros(15)
    for (S, S2), d in ST.items():
        for (dn, dm), val in d.items():
            H[MASKS.index(S), MASKS.index(S2)] += val*X**dn*Y**dm
    G = sp.zeros(15, 4)
    g = sp.eye(4)-sp.ones(4)/5
    for S in MASKS:
        dn, dm = step(S)
        G[MASKS.index(S), :] = 2*(sp.Matrix(bits(S)).T*g)*(X**dn*Y**dm-1)
    assert (H-H.T.applyfunc(conj)).applyfunc(sp.cancel) == sp.zeros(15)
    assert (H*G).applyfunc(sp.cancel) == sp.zeros(15, 4)
    return H, G


def orbit_vector(S, chars):
    e1, e2 = chars
    orb = defaultdict(int)
    for g, ch in [((0, 1, 2, 3), 1), ((1, 0, 2, 3), e1), ((0, 1, 3, 2), e2), ((1, 0, 3, 2), e1*e2)]:
        orb[perm_mask(g, S)] += ch
    return {k: v for k, v in orb.items() if v}


def sector_basis(chars):
    seen, vecs = set(), []
    for S in MASKS:
        o = orbit_vector(S, chars)
        if not o or frozenset(o) in seen:
            continue
        seen.add(frozenset(o))
        top = max(abs(v) for v in o.values())
        vecs.append((S, {k: sp.Rational(v, top) for k, v in o.items()}))
    return vecs


def check_symbols():
    H, G = full_symbol()
    dims = {}
    for chars in [(-1, -1), (1, -1), (-1, 1), (1, 1)]:
        B = sector_basis(chars)
        Bm = sp.zeros(15, len(B))
        for j, (_, o) in enumerate(B):
            for k, v in o.items():
                Bm[MASKS.index(k), j] = v
        Gs = ((Bm.T*Bm).inv()*Bm.T*G).applyfunc(sp.cancel)
        rank = Gs.subs({X: sp.Rational(2), Y: sp.Rational(3)}).rank()
        dims[chars] = (len(B), rank)
    assert dims == {(-1, -1): (1, 0), (1, -1): (3, 1), (-1, 1): (3, 1), (1, 1): (8, 2)}
    # 交叉扇区复现已有符号 (2T^2+TV-12)/20。
    Bc = sp.zeros(15, 1)
    for k, v in orbit_vector(5, (-1, -1)).items():
        Bc[MASKS.index(k), 0] = v/2  # 归一化交叉模式：系数 ±1/2
    cross = sp.cancel((Bc.T*H*Bc)[0])
    target = CROSS["cross_symbol"]().subs({T: X+1/X, V: Y+1/Y})
    assert sp.cancel(cross-target) == 0

    # 向量扇区 (+,-)：边型 e3（步 (1,-1)）、e13（步 (2,0)）、e123（步 (3,1)）。
    B = sector_basis((1, -1))
    assert [S for S, _ in B] == [4, 5, 7]
    Bm = sp.zeros(15, 3)
    for j, (_, o) in enumerate(B):
        for k, v in o.items():
            Bm[MASKS.index(k), j] = v
    Hs = (Bm.T*H*Bm).applyfunc(sp.cancel)
    gam = sp.Matrix([X/Y-1, X**2-1, X**3*Y-1])
    assert (Hs*gam).applyfunc(sp.cancel) == sp.zeros(3, 1)
    c = (4*X**4*Y+X**3*Y**2-X**3+22*X**2*Y-X*Y**2+X+4*Y)/(50*X**2*Y)
    A = Hs.adjugate().applyfunc(sp.cancel)
    for i in range(3):
        for j in range(3):
            assert sp.cancel(A[i, j]-c*gam[i]*conj(gam[j])) == 0
    # c = (2T^2+7-2 sin t sin s)/25：(X-1/X)(Y-1/Y) = -4 sin t sin s。
    assert sp.cancel(c-(4*(X**2+X**-2)+22+(X-1/X)*(Y-1/Y))/50) == 0
    # (-,+) 扇区的同一不变量为 c(t,-s)。
    B2 = sector_basis((-1, 1))
    assert [S for S, _ in B2] == [1, 5, 13]
    Bm2 = sp.zeros(15, 3)
    for j, (_, o) in enumerate(B2):
        for k, v in o.items():
            Bm2[MASKS.index(k), j] = v
    Hs2 = (Bm2.T*H*Bm2).applyfunc(sp.cancel)
    gam2 = sp.Matrix([X*Y-1, X**2-1, X**3/Y-1])
    assert (Hs2*gam2).applyfunc(sp.cancel) == sp.zeros(3, 1)
    A2 = Hs2.adjugate().applyfunc(sp.cancel)
    c2 = c.subs(Y, 1/Y)
    assert all(sp.cancel(A2[i, j]-c2*gam2[i]*conj(gam2[j])) == 0 for i in range(3) for j in range(3))
    # 规范不变局域场强 f1=(γ̄2,-γ̄1,0)、f2=(0,γ̄3,-γ̄2)：det(F†F)=|γ2|^2|γ|^2，
    # det(F†HF)=c|γ2|^2|γ|^4，于是两点函数矩阵 M=N^{-1} 满足 det M=|γ2|^2/c。
    f1 = sp.Matrix([conj(gam[1]), -conj(gam[0]), 0])
    f2 = sp.Matrix([0, conj(gam[2]), -conj(gam[1])])
    F = sp.Matrix.hstack(f1, f2)
    Fd = F.T.applyfunc(conj)
    assert (Fd*gam).applyfunc(sp.cancel) == sp.zeros(2, 1)
    gn = sp.expand(sum(gam[i]*conj(gam[i]) for i in range(3)))
    g2 = sp.expand(gam[1]*conj(gam[1]))
    FHF = (Fd*Hs*F).applyfunc(sp.expand)
    FF = (Fd*F).applyfunc(sp.expand)
    assert sp.cancel(FF.det()-g2*gn) == 0
    assert sp.cancel(FHF[0, 0]*FHF[1, 1]-FHF[0, 1]*FHF[1, 0]-c*g2*gn**2) == 0
    # 对 γ⊥ 中的 f、g：f†H⁺g=(F†f)†(F†HF)^{-1}(F†g)。M=(F†F)(F†HF)^{-1}(F†F) 的元素与第三个
    # Koszul 场强 f3=(γ̄3,0,-γ̄1) 的自关联，其分母只含 c 的分子多项式与单项式（不含 |γ|^2）。
    cnum = sp.Poly(4*X**4*Y+X**3*Y**2-X**3+22*X**2*Y-X*Y**2+X+4*Y, X, Y)
    adj = sp.Matrix([[FHF[1, 1], -FHF[0, 1]], [-FHF[1, 0], FHF[0, 0]]])
    detF = sp.expand(FHF[0, 0]*FHF[1, 1]-FHF[0, 1]*FHF[1, 0])
    Mm = (FF*adj*FF/detF).applyfunc(lambda e: sp.factor(sp.cancel(e)))
    f3 = sp.Matrix([conj(gam[2]), 0, -conj(gam[0])])
    assert sp.cancel((f3.T.applyfunc(conj)*gam)[0]) == 0
    u3 = (Fd*f3).applyfunc(sp.expand)
    m33 = sp.factor(sp.cancel((u3.T.applyfunc(conj)*adj*u3)[0]/detF))

    def denominator_only_c(e):
        den = sp.Poly(sp.fraction(sp.cancel(e))[1], X, Y)
        q, r = sp.div(den, cnum)
        if r.is_zero and q.is_monomial:
            return True
        return den.is_monomial  # 分母为单项式（无极点）亦可

    assert all(denominator_only_c(Mm[i, j]) for i in range(2) for j in range(2))
    assert denominator_only_c(m33)
    assert any(not sp.Poly(sp.fraction(sp.cancel(Mm[i, j]))[1], X, Y).is_monomial for i in range(2) for j in range(2))
    assert sp.cancel(Mm.det()-g2/c) == 0
    # (-,+) 扇区的同样恒等式：det(F'†F')=|γ'2|^2|γ'|^2，det(F'†H'F')=c(t,-s)|γ'2|^2|γ'|^4。
    F2 = sp.Matrix.hstack(sp.Matrix([conj(gam2[1]), -conj(gam2[0]), 0]), sp.Matrix([0, conj(gam2[2]), -conj(gam2[1])]))
    F2d = F2.T.applyfunc(conj)
    gn2 = sp.expand(sum(gam2[i]*conj(gam2[i]) for i in range(3)))
    g22 = sp.expand(gam2[1]*conj(gam2[1]))
    FHF2 = (F2d*Hs2*F2).applyfunc(sp.expand)
    assert sp.cancel((F2d*F2).applyfunc(sp.expand).det()-g22*gn2) == 0
    assert sp.cancel(FHF2[0, 0]*FHF2[1, 1]-FHF2[0, 1]*FHF2[1, 0]-c2*g22*gn2**2) == 0
    # 零动量：H_V=(3/5)[[-1,2,-1],[2,-4,2],[-1,2,-1]]，秩 1；核恰为常度规模式 u2⊗a、u2⊗b
    # （a=(1,1,0,0)、b=(0,0,1,1)，坐标 (0,1,2)、(1,1,1)）；唯一非零方向 (1,-2,1)，本征值 -18/5，
    # 与度规模式正交，是格点非度规边模。c(0,0)=3/5，故小动量处另一本征值 = c|γ|^2/(-18/5+o(1)) = O(k^2)。
    H0 = Hs.subs({X: 1, Y: 1})
    assert H0 == sp.Rational(3, 5)*sp.Matrix([[-1, 2, -1], [2, -4, 2], [-1, 2, -1]]) and H0.rank() == 1
    u2 = sp.Matrix([0, 0, 1, -1])
    metric = []
    for a in (sp.Matrix([1, 1, 0, 0]), sp.Matrix([0, 0, 1, 1])):
        h = u2*a.T+a*u2.T
        coords = []
        for _, o in B:
            vec = {k: v for k, v in o.items()}
            num = sum(v*(sp.Matrix(bits(k)).T*h*sp.Matrix(bits(k)))[0] for k, v in vec.items())
            coords.append(num/sum(v*v for v in vec.values()))
        metric.append(sp.Matrix(coords)/2)
    assert metric[0] == sp.Matrix([0, 1, 2]) and metric[1] == sp.Matrix([1, 1, 1])
    assert all(H0*v == sp.zeros(3, 1) for v in metric)
    lat = sp.Matrix([1, -2, 1])
    assert H0*lat == -sp.Rational(18, 5)*lat and all((lat.T*v)[0] == 0 for v in metric)
    assert c.subs({X: 1, Y: 1}) == sp.Rational(3, 5)
    return c


def check_vector_roots():
    """c 的零点：X=iZ 化为 4W^2+2σW-30=0（W=Z+1/Z，σ=sin s），故 X 纯虚、X^2=-Z^2<0。"""
    Z, W, sig, t, s = sp.symbols("Z W sigma t s", real=True)
    c = (4*X**4*Y+X**3*Y**2-X**3+22*X**2*Y-X*Y**2+X+4*Y)/(50*X**2*Y)
    # 4X^4+(Y-1/Y)X^3+22X^2-(Y-1/Y)X+4，再以 Y-1/Y=2iσ 代入
    p2 = sp.expand(4*X**4+2*sp.I*sig*X**3+22*X**2-2*sp.I*sig*X+4)
    assert sp.expand(sp.expand(c*50*X**2)-sp.expand(4*X**4+(Y-1/Y)*X**3+22*X**2-(Y-1/Y)*X+4)) == 0
    real = sp.expand(p2.subs(X, sp.I*Z))
    assert sp.expand(real-(4*Z**4+2*sig*Z**3-22*Z**2+2*sig*Z+4)) == 0
    assert sp.expand(real/Z**2-(4*W**2+2*sig*W-30).subs(W, Z+1/Z)) == 0
    r = sp.sqrt(sig**2+120)
    Wp, Wm = (-sig+r)/4, (-sig-r)/4
    # 两支 cosh E = |W|/2 = (r∓σ)/8；对 σ∈[-1,1] 落在 [5/4,3/2]，故 Z 实、|Z|<1 的根存在。
    assert sp.expand(4*Wp**2+2*sig*Wp-30) == 0 and sp.expand(4*Wm**2+2*sig*Wm-30) == 0
    for val in (-1, 1):
        assert {sp.nsimplify(Wp.subs(sig, val)), sp.nsimplify(-Wm.subs(sig, val))} == {sp.Rational(5, 2), 3}
    assert sp.simplify(Wp.subs(sig, 0)-sp.sqrt(30)/2) == 0
    # (r-σ)/8 单调减、(r+σ)/8 单调增：导数 (σ/r∓1)/8，|σ|<r。
    assert sp.simplify(sp.diff((r-sig)/8, sig)-(sig/r-1)/8) == 0
    # 能隙：cosh E=5/4 即 E=ln2（Z=1/2），cosh E=3/2 即 E=arcosh(3/2)=2 ln φ。
    assert sp.simplify(sp.cosh(sp.log(2))-sp.Rational(5, 4)) == 0
    phi = (1+sp.sqrt(5))/2
    assert sp.radsimp(sp.expand((phi**2+phi**-2)/2)) == sp.Rational(3, 2)  # cosh(2 ln φ)
    assert (Z+1/Z).subs(Z, sp.Rational(1, 2)) == sp.Rational(5, 2)
    # |Z|=e^{-E} 时 |W|/2=(e^{E}+e^{-E})/2=cosh E。
    Ev = sp.symbols("E", positive=True)
    assert sp.simplify(((Z+1/Z)/2).subs(Z, sp.exp(-Ev))-sp.cosh(Ev)) == 0
    # 环面上 c ≥ 1/5：2T^2+7-2 sin t sin s ≥ 5，在 t=s=±π/2 取到。
    cts = (2*(2*sp.cos(t))**2+7-2*sp.sin(t)*sp.sin(s))/25
    for v in (sp.pi/2, -sp.pi/2):
        assert sp.simplify(cts.subs({t: v, s: v})-sp.Rational(1, 5)) == 0
    return Wp, Wm


def check_tensor_signs():
    """常度规模式上：站点型反射把正偏振 u1u1-u2u2 变号、保持交叉 u1u2+u2u1；PT 型保持正偏振。"""
    u1 = sp.Matrix([1, -1, 0, 0])
    u2 = sp.Matrix([0, 0, 1, -1])
    plus = u1*u1.T-u2*u2.T
    cross = u1*u2.T+u2*u1.T
    out = {}
    for pi in [(2, 3, 0, 1), (3, 2, 1, 0), (0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2)]:
        P = sp.Matrix(4, 4, lambda i, j: int(pi[j] == i))
        # 边 e_S -> e_{πS}，平方边长不变：h' = P h P^T。
        signs = []
        for h in (plus, cross):
            img = P*h*P.T
            signs.append(1 if img == h else -1 if img == -h else 0)
        # 与边层次一致：q_S(h) = e_S^T h e_S，q'_{πS} = q_S。
        for h in (plus, cross):
            for S in MASKS:
                e = sp.Matrix(bits(S))
                assert ((P*e).T*(P*h*P.T)*(P*e))[0] == (e.T*h*e)[0]
        out[pi] = tuple(signs)
    assert out[(2, 3, 0, 1)] == out[(3, 2, 1, 0)] == (-1, 1)
    assert out[(0, 1, 2, 3)] == out[(1, 0, 3, 2)] == (1, 1)
    assert out[(1, 0, 2, 3)] == out[(0, 1, 3, 2)] == (1, -1)
    return out


class Box:
    """棋盘盒：n 以 N 周期，m 以 M 反周期，场点 n≡m (mod 2)，每点 15 个边分量。"""

    def __init__(self, N, M):
        self.N, self.M = N, M
        self.sites = [(n, m) for n in range(N) for m in range(M) if (n-m) % 2 == 0]
        self.var = [(x, S) for x in self.sites for S in MASKS]
        self.idx = {v: i for i, v in enumerate(self.var)}
        D = self.D = len(self.var)
        Q = fmpq_mat(D, D)
        for x in self.sites:
            for (S, S2), d in ST.items():
                for (dn, dm), val in d.items():
                    if val == 0:
                        continue
                    y, sg = self.red(x[0]+dn, x[1]+dm)
                    r = sp.Rational(val)
                    Q[self.idx[(x, S)], self.idx[(y, S2)]] += sg*fmpq(int(r.p), int(r.q))
        assert Q == Q.transpose()
        self.Q = Q
        g = sp.eye(4)-sp.ones(4)/5
        G = fmpq_mat(D, 4*len(self.sites))
        for a, x in enumerate(self.sites):
            for j in range(4):
                for S in MASKS:
                    coef = sp.Rational(2*(sp.Matrix(bits(S)).T*g)[j])
                    dn, dm = step(S)
                    for base, sgn in (((x[0]-dn, x[1]-dm), 1), (x, -1)):
                        y, sg = self.red(*base)
                        G[self.idx[(y, S)], 4*a+j] += fmpq(int(coef.p), int(coef.q))*sgn*sg
        self.G = G

    def red(self, n, m):
        n %= self.N
        sg = 1
        while m < 0:
            m += self.M
            sg = -sg
        while m >= self.M:
            m -= self.M
            sg = -sg
        return (n, m), sg

    def theta(self, pi, w):
        """x -> w - πx。边 (x,S) 映为基点 w-πx-e_{πS}、边型 πS 的边，场符号 +1。"""
        nw, mw = sum(w), w[0]+w[1]-w[2]-w[3]
        keeps = (pi[0] < 2) == (pi[1] < 2) and pi[0] < 2
        P = fmpq_mat(self.D, self.D)
        for (x, S), i in self.idx.items():
            S2 = perm_mask(pi, S)
            n2, m2 = step(S2)
            mpx = x[1] if keeps else -x[1]
            y, sg = self.red(nw-x[0]-n2, mw-mpx-m2)
            P[self.idx[(y, S2)], i] = sg
        return P

    def pointwise(self, pi):
        P = fmpq_mat(self.D, self.D)
        for (x, S), i in self.idx.items():
            P[self.idx[(x, perm_mask(pi, S))], i] = 1
        return P

    def mshift(self, dm):
        P = fmpq_mat(self.D, self.D)
        for (x, S), i in self.idx.items():
            y, sg = self.red(x[0], x[1]+dm)
            P[self.idx[(y, S)], i] = sg
        return P


def colbasis(A):
    R, rk = A.transpose().rref()
    return fmpq_mat([[R[j, i] for j in range(rk)] for i in range(A.nrows())])


def nullspace(B):
    R, rk = B.rref()
    n = B.ncols()
    piv, r = [], 0
    for j in range(n):
        if r < rk and R[r, j] != 0:
            piv.append(j)
            r += 1
    free = [j for j in range(n) if j not in piv]
    Xm = fmpq_mat(n, len(free))
    for c, f in enumerate(free):
        Xm[f, c] = 1
        for r, p in enumerate(piv):
            Xm[p, c] = -R[r, f]
    assert B*Xm == fmpq_mat(B.nrows(), len(free))
    return Xm


def inertia(M):
    co = M.charpoly().coeffs()
    z = next(i for i, c in enumerate(co) if c != 0)
    nz = [c for c in co if c != 0]
    changes = lambda L: sum((a > 0) != (b > 0) for a, b in zip(L, L[1:]))
    return changes(nz), changes([(-1)**i*c for i, c in enumerate(co) if c != 0]), z


class Sector:
    def __init__(self, box, chars):
        cols = []
        for x in box.sites:
            for _, o in sector_basis(chars):
                cols.append({box.idx[(x, k)]: fmpq(int(v.p), int(v.q)) for k, v in o.items()})
        Vm = fmpq_mat(box.D, len(cols))
        for j, c in enumerate(cols):
            for k, v in c.items():
                Vm[k, j] = v
        self.box, self.V = box, Vm
        self.Vp = (Vm.transpose()*Vm).inv()*Vm.transpose()
        self.K = Vm.transpose()*(-box.Q)*Vm          # 使交叉扇区可归一的整体符号
        self.Ga = colbasis(self.Vp*box.G)
        self.gauge_free = self.K.rank()+self.Ga.ncols() == self.K.nrows()
        if self.gauge_free and charpoly_psd(self.K):
            self.C = (self.K+self.Ga*self.Ga.transpose()).inv()  # 规范不变观测量上等于 K^+

    def half(self, lo, hi):
        """端点全在 [lo,hi] 内的边上的规范不变线性观测量（以原始场坐标的泛函表示）。"""
        hcols = []
        for j in range(self.V.ncols()):
            k = next(k for k in range(self.box.D) if self.V[k, j] != 0)
            (x, S) = self.box.var[k]
            if x[0] >= lo and x[0]+step(S)[0] <= hi:
                hcols.append(j)
        E = fmpq_mat(self.V.ncols(), len(hcols))
        for c, j in enumerate(hcols):
            E[j, c] = 1
        F = E*nullspace(self.Ga.transpose()*E)
        return self.Vp.transpose()*F

    def cov(self, w1, w2):
        return (self.V.transpose()*w1).transpose()*self.C*(self.V.transpose()*w2)

    def invariant(self, w):
        return (self.V.transpose()*w).transpose()*self.Ga == fmpq_mat(w.ncols(), self.Ga.ncols())


def raw(box, terms):
    w = fmpq_mat(box.D, 1)
    for x, S, cf in terms:
        y, sg = box.red(*x)
        w[box.idx[(y, S)], 0] += sg*cf
    return w


def check_box():
    box = Box(8, 4)
    assert box.D == 240 and box.Q.rank() == 176 and box.G.rank() == 64
    assert box.Q*box.G == fmpq_mat(240, 64)
    sec = {ch: Sector(box, ch) for ch in [(-1, -1), (1, -1), (-1, 1), (1, 1)]}
    inert = {ch: inertia(s.K) for ch, s in sec.items()}
    assert inert == {(-1, -1): (16, 0, 0), (1, -1): (32, 0, 16), (-1, 1): (32, 0, 16), (1, 1): (80, 16, 32)}
    assert all(s.gauge_free for s in sec.values())
    A, B = sec[(1, -1)], sec[(-1, 1)]
    # 站点型反射：交换 (+,-) 与 (-,+)，交叉扇区半正定（复现上一结论）。
    site = {pi: box.theta(pi, (0, 0, 0, 0)) for pi in [(2, 3, 0, 1), (3, 2, 1, 0)]}
    U12, U34 = box.pointwise((1, 0, 2, 3)), box.pointwise((0, 1, 3, 2))
    witness = {}
    for pi, P in site.items():
        assert P.transpose()*box.Q*P == box.Q and P*P == box.pointwise((0, 1, 2, 3))
        assert P*U12*P == U34
        Xc = sec[(-1, -1)]
        wc = Xc.half(0, 4)
        Rc = Xc.cov(P.transpose()*wc, wc)
        assert Rc == Rc.transpose() and charpoly_psd(Rc) and Rc.nrows() == 6
        wA, wB = A.half(0, 4), B.half(0, 4)
        assert wA.ncols() == wB.ncols() == 8
        assert B.invariant(P.transpose()*wA) and A.invariant(P.transpose()*wB)
        RAB = B.cov(P.transpose()*wA, wB)
        assert RAB == A.cov(P.transpose()*wB, wA).transpose() and RAB.rank() == 8
        # 整个向量扇区的反射矩阵 [[0,R],[R^T,0]]：符号差 (8,8)。
        full = fmpq_mat([[(RAB[i, j-8] if j >= 8 else 0) if i < 8 else (RAB[j, i-8] if j < 8 else 0)
                          for j in range(16)] for i in range(16)])
        assert inertia(full) == (8, 8, 0)
        # 显式局域见证：(+,-) 环 L_A 与 (-,+) 环 L_B。
        LA = raw(box, [((0, 0), 4, 1), ((1, -1), 7, 1), ((2, 0), 5, -1), ((0, 0), 5, -1)])
        LB = raw(box, [((0, 0), 1, 1), ((1, 1), 13, 1), ((2, 0), 5, -1), ((0, 0), 5, -1)])
        assert A.invariant(LA) and B.invariant(LB)
        r = B.cov(P.transpose()*LA, LB)[0, 0]
        assert r == A.cov(P.transpose()*LB, LA)[0, 0]
        assert A.cov(LA, LA)[0, 0] == fmpq(1643, 752)
        lam = -1 if r > 0 else 1
        value = 2*lam*r   # <Θ(L_A+λL_B), L_A+λL_B>，对角项因扇区正交为零
        witness[pi] = (r, value)
    assert witness[(2, 3, 0, 1)] == (fmpq(-17, 564), fmpq(-17, 282))
    assert witness[(3, 2, 1, 0)] == (fmpq(17, 564), fmpq(-17, 282))
    # 一致性检查（不是独立检验）：U12 在 (+,-) 上为 +1、在 (-,+) 上为 -1，R(Uf,Ug)=-R 等价于特征符号；
    # 引理的实质条件 Θ(12)Θ=(34) 由上面的 P*U12*P == U34 断言。
    P = site[(2, 3, 0, 1)]
    wA = A.half(0, 4)
    RAB = B.cov(P.transpose()*wA, B.half(0, 4))
    wBU = U12.transpose()*B.half(0, 4)
    assert B.cov(P.transpose()*(U12.transpose()*wA), wBU) == -RAB
    # PT 型反射：在 4|M 的反周期盒上，τ=m 平移 M/2 满足 τ^2=-1、ΘτΘ=τ^{-1}=-τ，反射矩阵符号差平衡。
    tau = box.mshift(2)
    assert tau*tau == -box.pointwise((0, 1, 2, 3))
    pt = {}
    cases = [((0, 1, 2, 3), (0, 0, 0, 0), (0, 4)), ((0, 1, 2, 3), (1, 0, 0, 0), (1, 4)),
             ((1, 0, 2, 3), (0, 0, 0, 0), (0, 4)), ((1, 0, 2, 3), (0, 0, 1, 0), (1, 4)),
             ((0, 1, 3, 2), (0, 0, 0, 0), (0, 4)), ((0, 1, 3, 2), (1, 0, 0, 0), (1, 4)),
             ((1, 0, 3, 2), (0, 0, 0, 0), (0, 4)), ((1, 0, 3, 2), (1, 1, 0, 0), (1, 5))]
    for pi, w, (lo, hi) in cases:
        Pt = box.theta(pi, w)
        assert Pt.transpose()*box.Q*Pt == box.Q and Pt*tau*Pt == -tau
        for ch in [(-1, -1), (1, -1), (-1, 1)]:
            s = sec[ch]
            wh = s.half(lo, hi)
            R = s.cov(Pt.transpose()*wh, wh)
            assert R == R.transpose()
            p, q, z = inertia(R)
            assert p == q and p > 0 and z == 0
            pt[(pi, sum(w), ch)] = (p, q)
    assert sorted(set(pt.values())) == [(2, 2), (3, 3), (4, 4)]
    return box, witness, pt


def verify():
    t0 = time.time()
    c = check_symbols()
    Wp, Wm = check_vector_roots()
    signs = check_tensor_signs()
    box, witness, pt = check_box()
    print("扇区维数（边模式, 规范秩）：交叉 (1,0)，向量 (3,1)+(3,1)，(+,+) (8,2)")
    print("向量扇区规范不变量 c(t,s)=(2T^2+7-2 sin t sin s)/25 ≥ 1/5；adj H = c γγ†；det M=|γ2|^2/c；M 与 f3 自关联的分母只含 c")
    print("零动量 H_V=(3/5)[[-1,2,-1],[2,-4,2],[-1,2,-1]]：核=常度规 u2⊗a、u2⊗b，非零方向 (1,-2,1) 本征值 -18/5")
    print("c 的零点 X=iZ，Z+1/Z=W，4W^2+2σW-30=0；两支 cosh E=(√(σ^2+120)∓σ)/8 ∈ [5/4,3/2]，E∈[ln2, arcosh(3/2)]")
    print("常度规张量符号（正偏振, 交叉）：", signs)
    print("8×4 反周期棋盘：扇区惯性（正,负,零）交叉 (16,0,0)，向量 (32,0,16)×2，(+,+) (80,16,32)")
    print("站点型反射：交叉半正定；向量扇区反射矩阵符号差 (8,8)；局域见证", witness)
    print("PT 型反射（τ^2=-1 的平衡符号差）：", sorted(set(pt.values())))
    print(f"用时 {time.time()-t0:.1f} 秒")


if __name__ == "__main__":
    verify()

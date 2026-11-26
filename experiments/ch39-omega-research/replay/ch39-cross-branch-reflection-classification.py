"""交叉偏振负比率支的 (-1)^n 规范本性、有隙第二支的合法消除与对角时间切片可实现反射的分类证书。

扇区：平直四维坐标链剖分 g=I-J/5 上的双奇交叉偏振，动量族 (t+s,t+s,t-s,t-s)，
即沿 (1,-1,0,0)、(0,0,1,-1) 平移不变的构型；格点坐标 n=x1+x2+x3+x4，m=x1+x2-x3-x4。
"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import fmpq, fmpq_mat


ROOT = Path(__file__).parent
T, V = sp.symbols("T V")


def bare_symbol():
    regge = run_path(str(ROOT/"ch39-regge-transverse-polarization.py"))
    K = sp.expand(-regge["cross_symbol"]())
    assert K == sp.expand((12-2*T**2-T*V)/20)
    return K


# 实空间模板：(Δn,Δm) -> 系数；T=z+1/z，V=w+1/w。
BARE = {(0, 0): fmpq(2, 5), (2, 0): fmpq(-1, 10), (-2, 0): fmpq(-1, 10),
           (1, 1): fmpq(-1, 20), (1, -1): fmpq(-1, 20),
           (-1, 1): fmpq(-1, 20), (-1, -1): fmpq(-1, 20)}
# 合法局域修正 K'=(6+V^2)/10-VT/4：(0,0) 处 4/5，(0,±2) 处 1/10，(±1,±1) 处 -1/4。
MODIFIED = {(0, 0): fmpq(4, 5), (0, 2): fmpq(1, 10), (0, -2): fmpq(1, 10),
            (1, 1): fmpq(-1, 4), (1, -1): fmpq(-1, 4), (-1, 1): fmpq(-1, 4), (-1, -1): fmpq(-1, 4)}
STENCIL = BARE


def symbol_of(stencil):
    t, s = sp.symbols("t s", real=True)
    expr = sum(sp.Rational(int(c.p), int(c.q))*sp.cos(a*t+b*s) for (a, b), c in stencil.items())
    return t, s, expr


def check_symbol_and_branches(K):
    t, s = sp.symbols("t s", real=True)
    symbol = sum(sp.Rational(int(c.p), int(c.q))*sp.cos(a*t+b*s) for (a, b), c in STENCIL.items())
    target = K.subs({T: 2*sp.cos(t), V: 2*sp.cos(s)})
    assert sp.simplify(sp.expand_trig(symbol-target)) == 0
    # 格点合法性：模板偏移同奇偶，等价于 K(-T,-V)=K(T,V)。
    assert all((a-b) % 2 == 0 for a, b in STENCIL)
    assert sp.expand(K.subs({T: -T, V: -V}, simultaneous=True)-K) == 0
    root = sp.sqrt(V**2+96)
    a, b = (-V+root)/4, (-V-root)/4
    assert sp.expand(K+(T-a)*(T-b)/10) == 0
    # 负比率支恰为 -V 处正支的像：b(V)=-a(-V)，权重 d(V)=c(-V)。
    assert sp.simplify(b+a.subs(V, -V)) == 0
    # 2cosh E_+(s)=a(V)、2cosh E_-(s)=-b(V)=a(-V)=2cosh E_+(s+π)；长波处 -b(2)=3，即 E_-(0)=arcosh(3/2)。
    assert sp.simplify(-b-a.subs(V, -V)) == 0
    assert sp.simplify(b.subs(V, 2)) == -3 and sp.simplify(a.subs(V, 2)) == 2
    c = 10/((a-b)*sp.sqrt(a**2-4))
    d = 10/((a-b)*sp.sqrt(b**2-4))
    assert sp.simplify(d-c.subs(V, -V)) == 0
    # a>2 <=> V<2，b<-2 <=> V>-2：(8+V)^2-(V^2+96)=16V-32。
    assert sp.expand((8+V)**2-(V**2+96)-(16*V-32)) == 0
    assert sp.expand((8-V)**2-(V**2+96)-(-16*V-32)) == 0
    # 非法的“消负支”反项：K*=(20-(8+V)T)/20 只有一个正根，但违背合法性。
    star = (20-(8+V)*T)/20
    assert sp.expand(star-K-(T-2)**2/10) == 0
    assert sp.expand(star.subs({T: -T, V: -V}, simultaneous=True)-star) != 0
    return a, b


def point_group():
    vectors = [tuple(int(i in S) for i in range(4)) for k in range(1, 5) for S in combinations(range(4), k)]
    vectors += [tuple(-x for x in v) for v in vectors]
    vset = set(vectors)
    g5 = sp.Matrix(4, 4, lambda i, j: 5*int(i == j)-1)
    short = [v for v in vectors if sum(x*x for x in v) in (1, 4)]
    group = []
    for cols in product(short, repeat=4):
        L = sp.Matrix(4, 4, lambda i, j: cols[j][i])
        if L.T*g5*L != g5:
            continue
        if all(tuple(L*sp.Matrix(v)) in vset for v in vectors):
            group.append(L)
    assert len(group) == 240
    n = sp.Matrix([[1, 1, 1, 1]])
    m = sp.Matrix([[1, 1, -1, -1]])
    u = [sp.Matrix([1, -1, 0, 0]), sp.Matrix([0, 0, 1, -1])]
    span = sp.Matrix.hstack(*u)
    results = []
    for L in group:
        if n*L != -n:
            continue
        # 扇区不变：L 保持 span(u1,u2)。
        if sp.Matrix.hstack(span, L*span).rank() != 2:
            continue
        if L*L != sp.eye(4):
            continue
        # 线性部分必为负坐标置换 -π，且 π 保持分块 {{1,2},{3,4}}。
        P = -L
        assert all(sorted(P.row(i)) == [0, 0, 0, 1] for i in range(4))
        # 扇区内 m 的线性部分：m∘L = μ m (模 span(u) 的零化子)。
        mu = (m*L)[0]
        assert m*L == mu*m and mu in (1, -1)
        # 平移 w 须满足 (I+L)w=0；读出 c=n(w) 的可能奇偶。
        # L=-π 时整数解 w 恰为在 π 的各轮换上取常值的向量，由 0/1 轮换指示向量生成。
        values = [int((n*sp.Matrix(w))[0]) for w in product((0, 1), repeat=4)
                  if (sp.eye(4)+L)*sp.Matrix(w) == sp.zeros(4, 1)]
        lattice_n = sp.gcd(values)
        # 场符号：交叉模式在边型 e1+e3,e2+e3,e1+e4,e2+e4 上的系数 (+,-,-,+)。
        mode = {(1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1, (0, 1, 0, 1): 1}
        signs = set()
        for e, val in mode.items():
            image = tuple(-x for x in L*sp.Matrix(e))  # 边无向，取 -L e 的正向代表
            signs.add(mode[image]*val)
        assert len(signs) == 1
        results.append((mu, int(lattice_n), signs.pop(), L))
    return group, results



def check_modified():
    """合法、局域、保持长波主项的修正 K'：有隙第二支被消去。"""
    t, s, expr = symbol_of(MODIFIED)
    Kp = (6+V**2)/10-V*T/4
    assert sp.simplify(sp.expand_trig(expr-Kp.subs({T: 2*sp.cos(t), V: 2*sp.cos(s)}))) == 0
    assert all((a-b) % 2 == 0 for a, b in MODIFIED) and sum(MODIFIED.values()) == 0
    assert sp.expand(Kp.subs({T: -T, V: -V}, simultaneous=True)-Kp) == 0
    # 长波主项与裸雷吉相同：t^2/2+s^2/10。
    eps = sp.symbols("eps")
    sub = {t: eps*t, s: eps*s}
    for e in (expr, symbol_of(BARE)[2]):
        ser = sp.series(e.subs(sub, simultaneous=True), eps, 0, 4).removeO()
        assert sp.expand(ser-eps**2*(t**2/2+s**2/10)) == 0
    # 环面非负：对 |T|≤2，K'≥(6+V^2)/10-|V|/2=(|V|-2)(|V|-3)/10≥0，仅 |V|=2、T=2sgnV（长波点）为零。
    Va = sp.symbols("Va", nonnegative=True)
    assert sp.expand((6+Va**2)/10-Va/2-(Va-2)*(Va-3)/10) == 0
    # 每个 0<V<2：K' 对 T 仿射、斜率 -V/4<0、唯一根 r=2(6+V^2)/(5V)，r-2=2(V-2)(V-3)/(5V)>0。
    r = sp.solve(Kp, T)[0]
    assert sp.simplify(r-2*(6+V**2)/(5*V)) == 0
    assert sp.simplify(r-2-2*(V-2)*(V-3)/(5*V)) == 0
    # 精确有理检验：八片链 V=±1，8×4 反周期棋盘。
    N = 8
    chains = {v: chain(MODIFIED, N, fmpq(v)) for v in (1, -1)}
    for k in range(N):
        assert chains[-1][k] == (-1)**k*chains[1][k]
    for v in (1, -1):
        H = fmpq_mat([[chains[v][(i+j) % N] for j in range(5)] for i in range(5)])
        assert charpoly_psd(H)
    link1 = fmpq_mat([[chains[1][(i+j+1) % N] for j in range(4)] for i in range(4)])
    assert charpoly_psd(link1) and not charpoly_psd(-link1)
    # V=1 的周期单链：唯一比率 x 满足 x+1/x=r=14/5，C_n 正比于 x^n+x^{N-n}。
    assert sp.simplify(r.subs(V, 1)) == sp.Rational(14, 5)
    x = (sp.Rational(14, 5)-sp.sqrt(sp.Rational(196, 25)-4))/2
    ratio = [sp.Rational(int(chains[1][k].p), int(chains[1][k].q)) for k in range(N)]
    for k in range(N):
        assert sp.simplify(ratio[k]*(1+x**N)-ratio[0]*(x**k+x**(N-k))) == 0
    sites, reduce, C = checkerboard(MODIFIED, 8, 4)
    site, _ = reflection_matrix(sites, reduce, C, lambda n, m: (-n, m), {0, 1, 2, 3, 4})
    assert charpoly_psd(site)
    pt, _ = reflection_matrix(sites, reduce, C, lambda n, m: (1-n, 1-m), {1, 2, 3, 4})
    assert not charpoly_psd(pt) and not charpoly_psd(-pt)
    return r


def chain(stencil, N, Vval):
    """固定横向 V 的周期单链：把模板沿 m 求和为 V 的多项式系数。"""
    w = {}
    for (a, b), c in stencil.items():
        # cos(b s) 对应 Chebyshev：2cos(b s)=2T_b(V/2)，按 e^{ibs} 拆为两半
        val = 2*sp.chebyshevt(abs(b), sp.Rational(int(Vval.p), int(Vval.q))/2)/2 if b != 0 else 1
        w[a] = w.get(a, 0)+sp.Rational(int(c.p), int(c.q))*val
    Q = fmpq_mat(N, N)
    for i in range(N):
        for a, val in w.items():
            Q[i, (i+a) % N] += fmpq(int(sp.fraction(val)[0]), int(sp.fraction(val)[1]))
    assert Q == Q.transpose()
    Cm = Q.inv()
    return [Cm[0, j] for j in range(N)]


def charpoly_psd(M):
    """对称有理矩阵半正定 <=> 特征多项式系数符号交错（允许零）。"""
    n = M.nrows()
    coeffs = M.charpoly().coeffs()  # 低次到高次
    return all((-1)**(n-k)*coeffs[k] >= 0 for k in range(n+1))


def check_chains():
    N = 8
    chains = {v: chain(BARE, N, fmpq(v)) for v in (1, -1, 0)}
    assert chains[0][:5] == [fmpq(35, 12), 0, fmpq(5, 6), 0, fmpq(5, 12)]
    for k in range(N):
        assert chains[-1][k] == (-1)**k*chains[1][k]
    for v, C in chains.items():
        H = fmpq_mat([[C[(i+j) % N] for j in range(5)] for i in range(5)])
        assert charpoly_psd(H)
    link = {v: fmpq_mat([[C[(i+j+1) % N] for j in range(4)] for i in range(4)]) for v, C in chains.items()}
    D = fmpq_mat([[(-1)**i*int(i == j) for j in range(4)] for i in range(4)])
    # 规范依赖性：link(-V) = -D link(V) D。
    assert link[-1] == -(D*link[1]*D)
    minors = {}
    for v in (1, -1):
        assert not charpoly_psd(link[v]) and not charpoly_psd(-link[v])
        minors[v] = link[v][0, 0]*link[v][1, 1]-link[v][0, 1]**2
        assert minors[v] == fmpq(-123725, 138384)
    return chains, minors


def checkerboard(stencil, N, M):
    """n mod N，m 以 M 反周期（φ(n,m+M)=-φ(n,m)），n≡m (mod 2)。"""
    sites = [(n, m) for n in range(N) for m in range(M) if (n-m) % 2 == 0]
    index = {x: i for i, x in enumerate(sites)}

    def reduce(n, m):
        n %= N
        sign = 1
        while m < 0:
            m += M; sign = -sign
        while m >= M:
            m -= M; sign = -sign
        return index[(n, m)], sign
    Q = fmpq_mat(len(sites), len(sites))
    for (n, m), i in index.items():
        for (a, b), c in stencil.items():
            j, sg = reduce(n+a, m+b)
            Q[i, j] += sg*c
    assert Q == Q.transpose()
    C = Q.inv()
    # C 只依赖差（带反周期符号），且在 (Δn,Δm)->(-Δn,-Δm)、(Δn,Δm)->(-Δn,Δm) 下不变。
    for x in sites:
        for y in sites:
            i, _ = reduce(*x); j, _ = reduce(*y)
            k, sg = reduce(x[0]-y[0], x[1]-y[1])
            o, _ = reduce(0, 0)
            assert C[i, j] == sg*C[k, o]
            k2, sg2 = reduce(y[0]-x[0], x[1]-y[1])
            assert C[k2, o]*sg2 == C[i, j]
    return sites, reduce, C


def reflection_matrix(sites, reduce, C, theta, half):
    pts = [x for x in sites if x[0] in half]
    R = fmpq_mat(len(pts), len(pts))
    for i, x in enumerate(pts):
        k, sg = reduce(*theta(*x))
        for j, y in enumerate(pts):
            l, _ = reduce(*y)
            R[i, j] = sg*C[k, l]
    assert R == R.transpose()
    return R, pts


def check_checkerboard():
    N, M = 8, 4
    # 反周期 m 盒且 4|M：全部允许动量 s=π(2j+1)/M 满足 4s≢0 (mod 2π)，定理 2 的前提成立。
    assert all((4*(2*j+1)) % (2*M) != 0 for j in range(M))
    sites, reduce, C = checkerboard(BARE, N, M)
    assert len(sites) == 16
    site, _ = reflection_matrix(sites, reduce, C, lambda n, m: (-n, m), {0, 1, 2, 3, 4})
    assert site.nrows() == 10 and charpoly_psd(site)
    pt, pts = reflection_matrix(sites, reduce, C, lambda n, m: (1-n, 1-m), {1, 2, 3, 4})
    assert pt.nrows() == 8
    assert not charpoly_psd(pt) and not charpoly_psd(-pt)
    # 显式见证：场点 (1,3) 的对角元 -17/12，(1,1) 的对角元 +17/12；故两种场符号都不正。
    assert pt[pts.index((1, 3)), pts.index((1, 3))] == fmpq(-17, 12)
    assert pt[pts.index((1, 1)), pts.index((1, 1))] == fmpq(17, 12)
    # 固定 m 的单步反射 (n,m)->(1-n,m) 把 n≡m 的站点送出场格点。
    assert all((1-n-m) % 2 == 1 for n, m in sites)
    return site, pt


def checkerboard_involutions():
    """有效二维棋盘上反转 n 的仿射映射 (n,m)->(c-n, ±m+d)：对合与保格点条件。"""
    out = set()
    for c in range(-3, 4):
        for d in range(-3, 4):
            for eps in (1, -1):
                f = lambda n, m: (c-n, eps*m+d)
                pts = [(n, m) for n in range(-4, 5) for m in range(-4, 5) if (n-m) % 2 == 0]
                if not all(f(*f(*x)) == x for x in pts):
                    continue
                if not all((f(*x)[0]-f(*x)[1]) % 2 == 0 for x in pts):
                    continue
                out.add((eps, c % 2, d % 2 if eps == -1 else d))
    # ε=+1：仅 d=0 且 c 偶（纯时间站点反射）；ε=-1：c≡d (mod 2) 的 PT 型。
    assert out == {(1, 0, 0), (-1, 0, 0), (-1, 1, 1)}
    return out


def verify():
    K = bare_symbol()
    a, b = check_symbol_and_branches(K)
    group, refl = point_group()
    assert len(refl) == 6
    site_type = [r for r in refl if r[0] == 1]
    pt_type = [r for r in refl if r[0] == -1]
    assert len(site_type) == 2 and len(pt_type) == 4
    assert all(r[1] % 2 == 0 and r[2] == 1 for r in site_type)
    assert sorted(r[1] for r in pt_type) == [1, 1, 1, 2]
    inv2d = checkerboard_involutions()
    chains, minors = check_chains()
    check_checkerboard()
    r = check_modified()
    print("点群阶", len(group), "；反转 n 且保持扇区的对合（μ, n(w) 的 gcd, 场符号）：",
          sorted((x[0], x[1], x[2]) for x in refl))
    print("二维棋盘上反转 n 的对合类型 (ε, c mod 2, d)：", sorted(inv2d))
    print("裸雷吉 V=0 八片 C0..C4：", chains[0][:5])
    print("V=±1 八片：C_n(-1)=(-1)^n C_n(1)，link(-V)=-D link(V) D；固定 m 单步 Hankel 二阶子式", minors[1])
    print("8×4 反周期棋盘：站点反射半正定；PT 型不定，对角见证 -17/12")
    print("合法修正 K'=(6+V^2)/10-VT/4：长波主项 t^2/2+s^2/10，唯一传递根", r, "；有隙第二支消失")


if __name__ == "__main__":
    verify()

"""交叉偏振负传递支的规范本性与对角时间切片可实现反射的完全分类证书。

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
STENCIL = {(0, 0): fmpq(2, 5), (2, 0): fmpq(-1, 10), (-2, 0): fmpq(-1, 10),
           (1, 1): fmpq(-1, 20), (1, -1): fmpq(-1, 20),
           (-1, 1): fmpq(-1, 20), (-1, -1): fmpq(-1, 20)}


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
    # 负支恰为 -V 处正支的像：b(V)=-a(-V)，权重 d(V)=c(-V)。
    assert sp.simplify(b+a.subs(V, -V)) == 0
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


def periodic_chain(N, Vval):
    Q = fmpq_mat(N, N)
    for i in range(N):
        Q[i, i] += fmpq(12, 20)-fmpq(4, 20)
        Q[i, (i+2) % N] += fmpq(-2, 20)
        Q[i, (i-2) % N] += fmpq(-2, 20)
        Q[i, (i+1) % N] += -Vval/20
        Q[i, (i-1) % N] += -Vval/20
    C = Q.inv()
    return [C[0, j] for j in range(N)]


def charpoly_psd(M):
    """对称有理矩阵半正定 <=> 特征多项式系数符号交错（允许零）。"""
    n = M.nrows()
    coeffs = M.charpoly().coeffs()  # 低次到高次
    return all((-1)**(n-k)*coeffs[k] >= 0 for k in range(n+1))


def check_chains():
    N = 8
    chains = {v: periodic_chain(N, fmpq(v)) for v in (1, -1, 0)}
    assert chains[0][:5] == [fmpq(35, 12), 0, fmpq(5, 6), 0, fmpq(5, 12)]
    for k in range(N):
        assert chains[-1][k] == (-1)**k*chains[1][k]
    site = {}
    for v, C in chains.items():
        H = fmpq_mat([[C[(i+j) % N] for j in range(5)] for i in range(5)])
        assert charpoly_psd(H)
        site[v] = H
    link = {v: fmpq_mat([[C[(i+j+1) % N] for j in range(4)] for i in range(4)]) for v, C in chains.items()}
    D = fmpq_mat([[(-1)**i*int(i == j) for j in range(4)] for i in range(4)])
    assert link[-1] == -(D*link[1]*D)
    minors = {}
    for v in (1, -1):
        assert not charpoly_psd(link[v]) and not charpoly_psd(-link[v])
        minors[v] = link[v][0, 0]*link[v][1, 1]-link[v][0, 1]**2
        assert minors[v] < 0
    return chains, minors


def checkerboard(N, M):
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
        for (a, b), c in STENCIL.items():
            j, sg = reduce(n+a, m+b)
            Q[i, j] += sg*c
    assert Q == Q.transpose()
    C = Q.inv()
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
    return R


def check_checkerboard():
    N, M = 8, 4
    sites, reduce, C = checkerboard(N, M)
    assert len(sites) == 16
    # 平移不变与 (Δn,Δm) 偶性：C 只依赖差，且在 Δn→-Δn、Δm→-Δm 下不变。
    site = reflection_matrix(sites, reduce, C, lambda n, m: (-n, m), {0, 1, 2, 3, 4})
    assert site.nrows() == 10 and charpoly_psd(site)
    pt = reflection_matrix(sites, reduce, C, lambda n, m: (1-n, 1-m), {1, 2, 3, 4})
    assert pt.nrows() == 8
    assert not charpoly_psd(pt) and not charpoly_psd(-pt)
    # 固定 m 的单步反射 (n,m)->(1-n,m) 把 n≡m 的站点送出场格点。
    assert all((1-n-m) % 2 == 1 for n, m in sites)
    # PT 型二次型的显式负见证：f=±e_x 组合，取最小对角/2×2 子式。
    witness = None
    for i in range(pt.nrows()):
        if pt[i, i] < 0:
            witness = ("diag", i, pt[i, i]); break
    if witness is None:
        for i in range(pt.nrows()):
            for j in range(i+1, pt.nrows()):
                det = pt[i, i]*pt[j, j]-pt[i, j]**2
                if det < 0:
                    witness = ("minor", (i, j), det); break
            if witness:
                break
    assert witness is not None
    return site, pt, witness


def verify():
    K = bare_symbol()
    a, b = check_symbol_and_branches(K)
    group, refl = point_group()
    # 分类：反转 n、保持扇区的对合恰六个线性部分。
    assert len(refl) == 6
    site_type = [r for r in refl if r[0] == 1]
    pt_type = [r for r in refl if r[0] == -1]
    assert len(site_type) == 2 and len(pt_type) == 4
    # 站点型：可行平移的 n(w) 全为偶数；场符号为 +1。
    assert all(r[1] % 2 == 0 and r[2] == 1 for r in site_type)
    # PT 型：恒等、(12)、(34) 允许 n(w) 为奇数，(12)(34) 只允许偶数；
    # 故任何半片（c 奇）反射都只能是 PT 型。
    assert sorted(r[1] for r in pt_type) == [1, 1, 1, 2]
    chains, minors = check_chains()
    site, pt, witness = check_checkerboard()
    print("点群阶", len(group), "；反转 n 且保持扇区的对合（μ, n(w) 的生成 gcd, 场符号）：",
          sorted((r[0], r[1], r[2]) for r in refl))
    print("V=0 八片协方差 C0..C4：", chains[0][:5])
    print("V=±1 八片：C_n(-1)=(-1)^n C_n(1)；不可实现的单步半片 Hankel 首二阶子式",
          minors[1], minors[-1])
    print("8×4 反周期棋盘：站点反射 10 阶矩阵半正定；PT 型 8 阶矩阵（任一场符号）不定，见证", witness)


if __name__ == "__main__":
    verify()

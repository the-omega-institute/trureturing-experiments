"""Exact lower bound for coordinate-permutation invariant completions at cutoff 56.

Every positive semidefinite completion W = Delta - S of the cutoff-56 problem that is
invariant under the six coordinate permutations (real symmetric or complex Hermitian)
has rank at least 12.  The proof splits W into its trivial, sign and standard isotypic
blocks, rank W = a + b + 2c.  A fixed matrix Gamma with inertia (4, 2, 0) satisfies
Std = Gamma + (2/9) Triv_ABCD (+) 6 Sign, which forces c >= 4; fixed entries force
a >= 2 and b >= 1; an exact resultant excludes (a, b, c) = (2, 1, 4).
"""
import sys
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import fmpq, fmpq_mat

sys.set_int_max_str_digits(0)
prev = run_path(str(Path(__file__).with_name("ch39-local-memory-rank.py")))
K1, K2 = fmpq(2, 9), fmpq(6)


def identity(n):
    return fmpq_mat([[int(i == j) for j in range(n)] for i in range(n)])


def column(values):
    return fmpq_mat([[v] for v in values])


def dot(u, v):
    return (u.transpose()*v)[0, 0]


def frame(vectors):
    return fmpq_mat([[v[i, 0] for v in vectors] for i in range(vectors[0].nrows())])


def inertia(matrix):
    """Exact (positive, negative, zero) counts of a rational symmetric matrix.

    The characteristic polynomial of a real symmetric matrix has only real roots,
    so Descartes' rule counts its positive roots exactly."""
    coeffs = matrix.charpoly().coeffs()
    zero = 0
    while coeffs[zero] == 0:
        zero += 1
    def changes(cs):
        cs = [c for c in cs if c != 0]
        return sum(1 for x, y in zip(cs, cs[1:]) if (x > 0) != (y > 0))
    rest = coeffs[zero:]
    pos = changes(rest)
    neg = changes([c*(-1)**k for k, c in enumerate(rest)])
    assert pos+neg+zero == matrix.nrows()
    return pos, neg, zero


def setup():
    _, _, _, g0, g1, g2 = prev["MEM"]["coefficients"]()
    delta = g2-g1*g0.solve(g1)
    points = sorted(prev["LOCAL"]["shell_grid"](1) | {(17, 17, 17)})
    n = len(points)
    assert n == 19 and delta == delta.transpose()
    dist = [[prev["LOCAL"]["squared"](p, q) for q in points] for p in points]
    index = {p: i for i, p in enumerate(points)}
    perms = list(permutations(range(3)))
    act = [[index[tuple(p[k] for k in g)] for p in points] for g in perms]
    pis = []
    for a in act:
        m = fmpq_mat(n, n)
        for i in range(n):
            m[a[i], i] = 1
        pis.append(m)
    assert all(delta[a[i], a[j]] == delta[i, j] for a in act for i in range(n) for j in range(n))
    orbit = lambda p: sorted(i for i, q in enumerate(points) if sorted(q) == sorted(p))
    A, B, C, D = orbit((1, 5, 29)), orbit((7, 17, 23)), orbit((11, 11, 25)), orbit((13, 13, 23))
    E = [index[(17, 17, 17)]]
    assert sorted(A+B+C+D+E) == list(range(n)) and [len(o) for o in (A, B, C, D)] == [6, 6, 3, 3]
    near = [(i, j) for i, j in combinations(range(n), 2) if dist[i][j] <= 56]
    assert len(near) == 18
    # Orbits of free positions; real symmetric generators and allowed antisymmetric ones.
    free = [(i, i) for i in range(n)]+near
    keys = sorted({min(tuple(sorted((a[i], a[j]))) for a in act) for i, j in free})
    sym, anti = [], []
    for i0, j0 in keys:
        m = fmpq_mat(n, n)
        for a in act:
            m[a[i0], a[j0]] = m[a[j0], a[i0]] = 1
        sym.append(m)
        if i0 != j0:
            h = fmpq_mat(n, n)
            ok = True
            for a in act:
                i, j = a[i0], a[j0]
                ok &= h[j, i] != 1
                h[i, j], h[j, i] = 1, -1
            if ok:
                anti.append(h)
    assert len(keys) == 9 and len(anti) == 3
    support = {(i, j) for i, j in free} | {(j, i) for i, j in free}
    assert all(m[i, j] == 0 for m in sym+anti for i in range(n) for j in range(n)
               if (i, j) not in support)
    sign = [1 if g in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1 for g in perms]
    chi = [2 if g == (0, 1, 2) else (0 if s < 0 else -1) for g, s in zip(perms, sign)]
    ptriv = sum((p for p in pis), fmpq_mat(n, n))*fmpq(1, 6)
    psgn = sum((p*s for p, s in zip(pis, sign)), fmpq_mat(n, n))*fmpq(1, 6)
    pstd = sum((p*c for p, c in zip(pis, chi)), fmpq_mat(n, n))*fmpq(1, 3)
    tau = pis[perms.index((1, 0, 2))]
    pu = pstd*(identity(n)+tau)*fmpq(1, 2)
    assert ptriv+psgn+pstd == identity(n)
    assert [ptriv.rank(), psgn.rank(), pstd.rank(), pu.rank()] == [5, 2, 12, 6]
    e = lambda i: column([int(k == i) for k in range(n)])
    def adjacency(left, right):
        m = fmpq_mat(n, n)
        for i, j in near:
            if (i in left and j in right) or (i in right and j in left):
                m[i, j] = m[j, i] = 1
        return m
    ra, nbc, nbd, mcd = adjacency(A, A), adjacency(B, C), adjacency(B, D), adjacency(C, D)
    va, vb = pu*e(A[0]), pu*e(A[1])
    aplus, aminus = va+ra*va, vb-ra*vb
    assert va-ra*va == fmpq_mat(n, 1)
    vc = pu*e(index[(11, 11, 25)])
    vd = mcd*vc
    b1 = nbc*vc
    assert nbd*vd == b1
    rb = pu*e(B[1])
    assert pu*e(B[0]) == b1*fmpq(1, 2)
    b2 = (rb-b1*(dot(b1, rb)/dot(b1, b1)))*2
    std = [aplus, b1, vc, vd, aminus, b2]
    triv = [column([int(k in o) for k in range(n)]) for o in (A, B, C, D, E)]
    sgn = [psgn*e(A[0]), psgn*e(B[0])]
    fs, ft, fg = frame(std), frame(triv), frame(sgn)
    assert pu*fs == fs and ptriv*ft == ft and psgn*fg == fg
    assert [fs.rank(), ft.rank(), fg.rank()] == [6, 5, 2]
    for vs in (std, triv, sgn):
        assert all(dot(vs[i], vs[j]) == 0 for i in range(len(vs)) for j in range(i))
    norms = [[dot(v, v) for v in vs] for vs in (std, triv, sgn)]
    assert norms == [[fmpq(4, 3), fmpq(4, 3), fmpq(2, 3), fmpq(2, 3), 1, 1], [6, 6, 3, 3, 1],
                     [fmpq(1, 6), fmpq(1, 6)]]
    assert [norms[0][i]/norms[1][i] for i in range(4)] == [K1]*4
    assert [norms[0][4+i]/norms[2][i] for i in range(2)] == [K2]*2
    return dict(g0=g0, delta=delta, points=points, dist=dist, act=act, sym=sym, anti=anti,
                keys=keys, fs=fs, ft=ft, fg=fg, A=A, B=B, C=C, D=D, E=E)


def blocks(data, w):
    fs, ft, fg = data["fs"], data["ft"], data["fg"]
    return ft.transpose()*w*ft, fg.transpose()*w*fg, fs.transpose()*w*fs


def gamma(data, w):
    t, s, l = blocks(data, w)
    out = fmpq_mat(l)
    for i in range(4):
        for j in range(4):
            out[i, j] -= K1*t[i, j]
    for i in range(2):
        for j in range(2):
            out[4+i, 4+j] -= K2*s[i, j]
    return out


def excluded_214(gm, t0, s0):
    """Exact algebra: no PSD rank-4 standard block with trivial rank 2 and sign rank 1."""
    q = lambda x: sp.Rational(int(x.p), int(x.q))
    hq = gm.inv()
    H = sp.Matrix(6, 6, lambda i, j: q(hq[i, j]))
    a, tau = sp.symbols("a tau", positive=True)
    r = [q(t0[0, j]) for j in (1, 2, 3)]
    h = [q(t0[j, 4]) for j in (1, 2, 3)]
    e0, sigma = q(t0[0, 4]), q(s0[0, 1])
    k1, k2 = sp.Rational(2, 9), sp.Integer(6)
    u = sp.Matrix([a]+r+[0, 0])
    g = sp.Matrix([0]+[h[j]-e0/a*r[j] for j in range(3)]+[0, 0])
    e5, e6 = sp.Matrix([0]*4+[1, 0]), sp.Matrix([0]*5+[1])
    form = lambda x, y: sp.expand((x.T*H*y)[0])
    n11, n12 = sp.expand(a/k1+form(u, u)), form(u, g)
    u5, u6, g5, g6 = form(u, e5), form(u, e6), form(g, e5), form(g, e6)
    # tau = s1^2; E1 = tau*(N11 N33 - N13^2), E2 = s1*(N11 N23 - N12 N13).
    e1 = sp.expand(n11*(tau/k2+tau**2*H[4, 4]+2*sigma*tau*H[4, 5]+sigma**2*H[5, 5])
                   - (tau*u5+sigma*u6)**2)
    e2 = sp.expand(n11*(tau*g5+sigma*g6)-n12*(tau*u5+sigma*u6))
    assert sp.Poly(e2, tau).degree() == 1 and sp.Poly(e1, tau).degree() == 2
    c1, c0 = e2.coeff(tau, 1), e2.coeff(tau, 0)
    d2, d1, d0 = e1.coeff(tau, 2), e1.coeff(tau, 1), e1.coeff(tau, 0)
    res = sp.together(d2*c0**2-d1*c0*c1+d0*c1**2)
    num, den = sp.fraction(res)
    assert len(sp.Poly(den, a).terms()) == 1          # denominator is c * a^k, nonzero for a > 0
    p = sp.Poly(sp.expand(num), a)
    n11num, n11den = sp.fraction(sp.together(n11))
    n12num, n12den = sp.fraction(sp.together(n12))
    assert len(sp.Poly(n11den, a).terms()) == 1 and len(sp.Poly(n12den, a).terms()) == 1
    f11, f12 = sp.Poly(sp.expand(n11num), a), sp.Poly(sp.expand(n12num), a)
    quotient, remainder = sp.div(p, f11)
    assert remainder.is_zero and p.degree() == 6 and f11.degree() == 2
    assert quotient.count_roots() == 0                  # no real root at all (exact Sturm count)
    assert sp.gcd(f11, f12).degree() == 0              # N11 = 0 forces N12 != 0
    return p.degree(), f11.count_roots(0, sp.oo)


def verify():
    data = setup()
    delta, sym, anti = data["delta"], data["sym"], data["anti"]
    zero = fmpq_mat(6, 6)
    # Std - (2/9) Triv_ABCD (+) 6 Sign does not depend on the free correction.
    assert all(gamma(data, m) == zero for m in sym+anti)
    gm = gamma(data, delta)
    assert gm.det() != 0 and inertia(gm) == (4, 2, 0)
    t0, s0, _ = blocks(data, delta)
    tri = [blocks(data, m)[0] for m in sym+anti]
    sg = [blocks(data, m)[1] for m in sym+anti]
    # Fixed trivial block rows {A, E} x columns {B, C, D} of rank two: a >= 2.
    rows, cols = [0, 4], [1, 2, 3]
    assert all(m[i, j] == 0 for m in tri for i in rows for j in cols)
    assert fmpq_mat([[t0[i, j] for j in cols] for i in rows]).rank() == 2
    # Fixed nonzero sign off-diagonal: b >= 1.
    assert all(m[0, 1] == 0 for m in sg) and s0[0, 1] != 0
    # Only the diagonal entries of rows A and E of the trivial block are free.
    assert all(m[0, j] == 0 and m[4, j] == 0 for m in tri for j in range(1, 4))
    assert all(m[0, 4] == 0 for m in tri)
    degree, n11_positive_roots = excluded_214(gm, t0, s0)
    assert (degree, n11_positive_roots) == (6, 1)
    # Exact rank additivity rank W = a + b + 2c at a correction with singular blocks.
    x0_shift = [blocks(data, m) for m in sym]
    w = fmpq_mat(delta)
    key = data["keys"]
    k_a, k_e = key.index((0, 0)), key.index((11, 11))
    sa = x0_shift[k_a][1][0, 0]
    target = s0[0, 1]**2/s0[1, 1]
    xa = (s0[0, 0]-target)/sa
    w = w-sym[k_a]*xa
    t1, _, _ = blocks(data, w)
    tm = fmpq_mat([[t1[i, j] for j in range(4)] for i in range(4)])
    col = fmpq_mat([[t1[i, 4]] for i in range(4)])
    xe = t1[4, 4]-(col.transpose()*tm.solve(col))[0, 0]
    w = w-sym[k_e]*xe
    tt, ss, ll = blocks(data, w)
    assert (tt.rank(), ss.rank()) == (4, 1)
    assert w.rank() == tt.rank()+ss.rank()+2*ll.rank() == 4+1+2*6
    print("Gamma inertia", inertia(gm), "so the standard block has rank c >= 4")
    print("a >= 2, b >= 1; resultant degree", degree, "with only N11 roots, positive:",
          n11_positive_roots)
    print("Every coordinate-permutation invariant PSD completion has rank >= 12;",
          "invariant realizations need at least", 19+12, "auxiliary states")


if __name__ == "__main__":
    verify()

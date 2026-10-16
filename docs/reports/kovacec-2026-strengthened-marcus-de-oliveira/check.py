"""Exact check of a counterexample to Kovacec's Conjecture 7 (DMUC Preprint 26-11, 2026),
the strengthened Marcus-de Oliveira conjecture.

a = (-3-i, -5-2i, 4-3i), b = (-2+i, 4+3i, -2-5i), Q = (1/7)[[-3,-2,6],[-2,-6,-3],[-6,3,-2]].
Checks (exact rational arithmetic, Gaussian rationals as pairs):
  1. Q Q^T = I, so U = Q is unitary;
  2. det(D_a + Q D_b Q^T) = (2864 + 1840 i)/49;
  3. the six sigma-points and the nine transposition edges of the graph G(a, b);
  4. the polygon P with vertices z3, z5, z2, r, z6, z4 (r = (7672/155, 4424/155)) is x-monotone
     with lower boundary L and upper boundary U as listed, L < U strictly inside, and every
     edge of G lies in P (piecewise-affine endpoint checks);
  5. the determinant lies strictly below P: 12 Re d - 280 - 11 Im d = 408/49 > 0 with
     7672/155 < Re d < 280 (on that interval L(x) = (12x - 280)/11);
  6. the family Q_s (s in [17/2, 9]) stays unitary and escapes (sampled exactly at rational s
     with sqrt(s) rational, and symbolically through the stated rational functions).
Prediction registered before the run: every check holds.
Usage: python3 check.py   (standard library only)
"""
import sys
from fractions import Fraction as F
from itertools import permutations


def cmul(x, y):
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def cadd(x, y):
    return (x[0] + y[0], x[1] + y[1])


def det3(M):
    def t(i, j, k):
        return cmul(cmul(M[0][i], M[1][j]), M[2][k])
    s = (F(0), F(0))
    for (i, j, k), sg in [((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1), ((0, 2, 1), -1), ((1, 0, 2), -1), ((2, 1, 0), -1)]:
        v = t(i, j, k)
        s = cadd(s, (sg * v[0], sg * v[1]))
    return s


a = [(F(-3), F(-1)), (F(-5), F(-2)), (F(4), F(-3))]
b = [(F(-2), F(1)), (F(4), F(3)), (F(-2), F(-5))]


def check_unitary_det(Q):
    n = 3
    ok = all(sum(Q[i][k] * Q[j][k] for k in range(n)) == (1 if i == j else 0) for i in range(n) for j in range(n))
    # M = D_a + Q D_b Q^T (Q real)
    M = [[(F(0), F(0))] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            s = (F(0), F(0))
            for k in range(n):
                c = Q[i][k] * Q[j][k]
                s = cadd(s, (c * b[k][0], c * b[k][1]))
            if i == j:
                s = cadd(s, a[i])
            M[i][j] = s
    return ok, det3(M)


def main():
    ok = True
    Q = [[F(x, 7) for x in row] for row in [[-3, -2, 6], [-2, -6, -3], [-6, 3, -2]]]
    unit, d = check_unitary_det(Q)
    print("Q Q^T = I:", unit, "  det =", d)
    ok &= unit and d == (F(2864, 49), F(1840, 49))
    perms = list(permutations(range(3)))  # lexicographic: 123,132,213,231,312,321
    z = []
    for p in perms:
        v = (F(1), F(0))
        for i in range(3):
            v = cmul(v, cadd(a[i], b[p[i]]))
        z.append(v)
    print("sigma-points:", [(int(x), int(y)) for x, y in z])
    expected = [(-30, -50), (280, 280), (-130, 10), (-28, -56), (232, 376), (24, -20)]
    ok &= [(int(x), int(y)) for x, y in z] == expected and all(x.denominator == 1 and y.denominator == 1 for x, y in z)
    edges = []
    for i, p in enumerate(perms):
        for j, q in enumerate(perms):
            if i < j and sum(1 for k in range(3) if p[k] != q[k]) == 2:
                edges.append((i, j))
    print("edges (1-based):", [(i + 1, j + 1) for i, j in edges])
    ok &= sorted((i + 1, j + 1) for i, j in edges) == sorted([(1, 3), (1, 6), (1, 2), (2, 5), (2, 4), (3, 5), (3, 4), (4, 6), (5, 6)])
    r = (F(7672, 155), F(4424, 155))

    def L(x):
        if x <= -28: return (-11 * x - 1260) / 17
        if x <= 24: return (9 * x - 476) / 13
        if x <= F(7672, 155): return (99 * x - 3416) / 52
        return (12 * x - 280) / 11

    def U(x):
        return (183 * x + 25600) / 181 if x <= 232 else -2 * x + 840

    # boundary pieces pass through the listed vertices
    for v in [z[2], z[3], z[5], r, z[1]]:
        ok &= L(v[0]) == v[1]
    for v in [z[2], z[4], z[1]]:
        ok &= U(v[0]) == v[1]
    breaks = sorted({F(-130), F(-28), F(24), F(7672, 155), F(232), F(280)})
    ok &= all(L(x) < U(x) for x in [(breaks[i] + breaks[i + 1]) / 2 for i in range(len(breaks) - 1)] + breaks[1:-1])
    # every edge lies in P: split at breakpoints, check endpoints of each affine piece
    def in_P(p):
        return F(-130) <= p[0] <= 280 and L(p[0]) <= p[1] <= U(p[0])
    edge_ok = True
    for i, j in edges:
        p, q = z[i], z[j]
        lo, hi = sorted([p[0], q[0]])
        pts = [lo] + [x for x in breaks if lo < x < hi] + [hi]
        for x in pts:
            t = (x - p[0]) / (q[0] - p[0]) if q[0] != p[0] else F(0)
            y = p[1] + t * (q[1] - p[1])
            edge_ok &= in_P((x, y))
    print("all nine edges inside P:", edge_ok)
    ok &= edge_ok
    X, Y = d
    esc = 12 * X - 280 - 11 * Y
    print("Re d =", X, " 12 Re d - 280 - 11 Im d =", esc, " in (7672/155, 280):", F(7672, 155) < X < 280, " L(Re d) - Im d =", L(X) - Y)
    ok &= esc == F(408, 49) and F(7672, 155) < X < 280 and Y < L(X)
    # family: s with rational square root, s in [17/2, 9]
    fam_ok = True
    for sq in [F(3), F(29, 10), F(59, 20), F(117, 40)]:
        s = sq * sq
        if not (F(17, 2) <= s <= 9):
            continue
        c = 1 / (s + 5)
        Qs = [[c * (-2 * sq), c * (-4), c * (s + 3)], [c * (5 - s), c * (-4 * sq), c * (-2 * sq)], [c * (-4 * sq), c * (s - 3), c * (-4)]]
        u, ds = check_unitary_det(Qs)
        Xs = 8 * (29 * s * s - 102 * s + 1) / (s + 5) ** 2
        Ys = 8 * (47 * s * s - 370 * s + 443) / (s + 5) ** 2
        good = u and ds == (Xs, Ys) and F(7672, 155) < Xs < 280 and Ys < L(Xs)
        fam_ok &= good
        print(f"   s = {s}: unitary {u}, det = ({float(ds[0]):.4f}, {float(ds[1]):.4f}), below P: {Ys < L(Xs)}")
    ok &= fam_ok
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

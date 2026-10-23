#!/usr/bin/env python3
"""Checks for O. Coskun and A. Eden, arXiv:2511.09395v1, Remark 5.9 (octonion clause), exact arithmetic (sympy).

V7 = R^7 with the coordinate product of Definition 5.1: c1 = a1 b1 - sum_{i>=2} ai bi + (a2 b3 - a3 b2) + (a4 b5 - a5 b4)
+ (a6 b7 - a7 b6), ci = a1 bi + b1 ai (i >= 2). Octonions: pairs of quaternions with Baez (2002), (2):
(a, b)(c, d) = (ac - d b*, a* d + c b). A push-forward is a pair of real-linear maps nu: V7 -> O, pi: O -> V7 with
a * b = pi(nu(a) nu(b)).
Fidelity: Section 3.1 (3)-(4) for V3 (the quaternion push-forward with the source's pi and nu). Lemmas: the octonion
norm is multiplicative; L_{e2} on V7 has rank 2 and the kernel basis e2 + e3, e4, ..., e7; L_a has rank 2 for random
nonzero a orthogonal to e1 in V3, V5, V7, V9. Exits 1 on failure.
"""
import random
import sys
import sympy as sp

ok = True

def qmul(p, q):
    a1, b1, c1, d1 = p; a2, b2, c2, d2 = q
    return (a1*a2 - b1*b2 - c1*c2 - d1*d2, a1*b2 + b1*a2 + c1*d2 - d1*c2, a1*c2 - b1*d2 + c1*a2 + d1*b2, a1*d2 + b1*c2 - c1*b2 + d1*a2)

def qconj(p):
    return (p[0], -p[1], -p[2], -p[3])

def omul(x, y):
    a, b = x; c, d = y
    left = tuple(s - t for s, t in zip(qmul(a, c), qmul(d, qconj(b))))
    right = tuple(s + t for s, t in zip(qmul(qconj(a), d), qmul(c, b)))
    return (left, right)

def vid(n, a, b):
    """Coordinate product of the Vidinli algebra of odd dimension n (indices 0..n-1, unit e_1 = index 0)."""
    om = sum(a[2*k + 1]*b[2*k + 2] - a[2*k + 2]*b[2*k + 1] for k in range((n - 1) // 2))
    c1 = a[0]*b[0] - sum(a[i]*b[i] for i in range(1, n)) + om
    return [sp.expand(c1)] + [sp.expand(a[0]*b[i] + b[0]*a[i]) for i in range(1, n)]

A = sp.symbols('a1:4'); B = sp.symbols('b1:4')
abar = (A[0], -A[1], -A[2])
cross = (abar[1]*B[2] - abar[2]*B[1], abar[2]*B[0] - abar[0]*B[2], abar[0]*B[1] - abar[1]*B[0])
Rc = (cross[0], -cross[2], cross[1])
dot = sum(x*y for x, y in zip(abar, B))
geo = [sp.expand(dot - Rc[0]), sp.expand(-Rc[1]), sp.expand(-Rc[2])]
coord_ok = all(sp.expand(g - c) == 0 for g, c in zip(geo, vid(3, A, B)))
nu = lambda a: (a[0], 0, a[1], a[2])
pi = lambda x: (x[0] + x[1], x[2], x[3])
push = pi(qmul(nu(A), nu(B)))
push_ok = all(sp.expand(p - c) == 0 for p, c in zip(push, vid(3, A, B)))
print(f"fidelity Section 3.1 (3): the geometric formula equals the coordinate product on V3: {coord_ok}")
print(f"fidelity Section 3.1 (4): pi(nu(a) nu(b)) equals the V3 product for the source's pi and nu: {push_ok}")
ok &= coord_ok and push_ok
X = sp.symbols('x0:8'); Y = sp.symbols('y0:8')
ox = (X[:4], X[4:]); oy = (Y[:4], Y[4:])
pr = omul(ox, oy)
n2 = lambda o: sum(v**2 for v in o[0]) + sum(v**2 for v in o[1])
norm_ok = sp.expand(n2(pr) - n2(ox) * n2(oy)) == 0
print(f"Baez (2) octonions: |xy|^2 = |x|^2 |y|^2 identically in 16 real variables: {norm_ok}")
ok &= norm_ok
e2 = [0, 1, 0, 0, 0, 0, 0]
Bv = sp.symbols('b1:8')
L = sp.Matrix([[sp.diff(c, bj) for bj in Bv] for c in vid(7, e2, Bv)])
ker = L.nullspace()
kb = [sp.Matrix([0, 1, 1, 0, 0, 0, 0])] + [sp.Matrix([0]*k + [1] + [0]*(6 - k)) for k in range(3, 7)]
ker_ok = L.rank() == 2 and len(ker) == 5 and all(L * v == sp.zeros(7, 1) for v in kb) and sp.Matrix.hstack(*kb).rank() == 5
print(f"V7: L_e2 has rank {L.rank()}, kernel dimension {len(ker)}, kernel basis e2 + e3, e4, e5, e6, e7: {ker_ok}")
ok &= ker_ok
random.seed(20261023)
for n in (3, 5, 7, 9):
    Bn = sp.symbols(f'b1:{n + 1}')
    ranks = set()
    for _ in range(20):
        a = [0] + [sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for _ in range(n - 1)]
        if all(v == 0 for v in a):
            continue
        Ln = sp.Matrix([[sp.diff(c, bj) for bj in Bn] for c in vid(n, a, Bn)])
        ranks.add(Ln.rank())
    print(f"V{n}: rank of L_a for 20 random nonzero a orthogonal to e1: {sorted(ranks)}; necessary bound N >= 2d - 2 = {2*n - 2}")
    ok &= ranks == {2}
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)

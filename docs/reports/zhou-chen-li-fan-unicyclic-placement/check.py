#!/usr/bin/env python3
"""Exact checks for Conjecture 1 (odd branch), arXiv:2610.38539v1.

Uses only Python integers and SymPy exact arithmetic. No floating point.
The symbolic identities support the family argument; the rational certificates
independently certify the a=0 strict comparison. This is not a Lean proof.
"""
from itertools import combinations
from functools import reduce
from math import gcd, lcm
import sympy as S


def graph(a, high, leaf):
    edges = {(i, (i + 1) % 18) for i in range(18)}
    n = 18
    for i in range(6):
        for _ in range(a + (i in high)):
            edges.update({(3 * i, n), (n, n + 1)})
            n += 2
    edges.add((3 * leaf, n))
    n += 1
    return n, sorted(tuple(sorted(e)) for e in edges)


def matrix(n, edges):
    m = S.zeros(n)
    for u, v in edges:
        m[u, v] = m[v, u] = 1
    return m


def connected(n, edges):
    seen = {0}
    while True:
        more = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if more == seen:
            return len(seen) == n
        seen = more


def is_dissociation(n, edges, chosen):
    degrees = [0] * n
    for u, v in edges:
        if u in chosen and v in chosen:
            degrees[u] += 1
            degrees[v] += 1
    return max(degrees, default=0) <= 1


def brute_beta(n, edges):
    for k in range(n, -1, -1):
        if any(is_dissociation(n, edges, set(c)) for c in combinations(range(n), k)):
            return k
    raise AssertionError('empty set must be admissible')


def cycle_beta(m):
    full = (1 << m) - 1
    best = 0
    for mask in range(1 << m):
        left = ((mask << 1) | (mask >> (m - 1))) & full
        right = (mask >> 1) | ((mask & 1) << (m - 1))
        if not mask & left & right:
            best = max(best, mask.bit_count())
    return best


def primitive_integer_vector(v):
    den = lcm(*(x.q for x in v))
    ints = [int(x * den) for x in v]
    common = reduce(gcd, ints)
    return S.Matrix([x // common for x in ints])


def exact_certificate(label, a, lower, c):
    n, edges = graph(a, {0, 3} if label == 'A' else {0, 2}, 1 if label == 'A' else 4)
    adj = matrix(n, edges)
    sign = 1 if lower else -1
    v = primitive_integer_vector((sign * (adj - c * S.eye(n))).inv() * S.ones(n, 1))
    assert min(v) > 0
    residual = sign * (c.q * adj * v - c.p * v)
    assert min(residual) > 0
    print(f'{label}_{a} n={n} edges={edges}')
    print(f'{label}_{a} positive_integer_vector={list(v)}')
    print(f'{label}_{a} {"100*Adj*v-213*v" if lower else "213*v-100*Adj*v"}={list(residual)}')
    return adj, v


# Fidelity: explicit base case in the proof of source Theorem 1, main.tex line 270.
base_edges = [(0, 1), (1, 2), (0, 2), (0, 3)]
assert brute_beta(4, base_edges) == 3
print('FIDELITY Theorem 1 proof: beta(C3 with one pendant leaf)=3>floor(8/3)=2')
for m in range(3, 19):
    b = cycle_beta(m)
    assert b == (2 * m) // 3
    print(f'FIDELITY cycle order={m} beta={b}; source Theorem 1 bounds hold')
core_bound = cycle_beta(18)
assert core_bound == 12

# Actual graphs, not just reduced matrices: membership and exact dissociation number.
for a in range(7):
    for label, high, leaf in [('A', {0, 3}, 1), ('B', {0, 2}, 4)]:
        n, edges = graph(a, high, leaf)
        assert n == 12 * a + 23 and len(edges) == n and connected(n, edges)
        selected = set(range(n)) - {3 * i for i in range(6)}
        assert is_dissociation(n, edges, selected)
        # Any dissociation set restricts to one in the induced C18, so it
        # has at most 12+(n-18) elements; the explicit set attains this.
        phi = len(selected)
        assert phi == core_bound + n - 18 == 12 * a + 17
        assert (2 * n) // 3 < phi < n - 5
        assert phi % 2 == 1 and ((n - 18) // 2) % 6 == 2
        print(f'MEMBERSHIP {label}_{a}: n={n}, edges={len(edges)}, phi={phi}, strict_range=True')

# Every allowed low-hub leaf in the opposite-high graph is equivalent.
for leaf in [1, 2, 4, 5]:
    assert any({(s * i + b) % 6 for i in {0, 3}} == {0, 3}
               and (s * 1 + b) % 6 == leaf for s in [-1, 1] for b in range(6))
print('SHAPE all four low-hub leaf placements in A are dihedrally equivalent')

# Principal-minor determinant identity, independent variables.
x, h, ell, t, a, u, v, w, q = S.symbols('x h ell t a u v w q')
cycle = S.zeros(6)
for i in range(6):
    cycle[i, (i + 1) % 6] = cycle[(i + 1) % 6, i] = 1
ma = cycle + S.diag(h, ell, 0, h, 0, 0)
mb = cycle + S.diag(h, 0, h, 0, ell, 0)
pa = (x * S.eye(6) - ma).det(method='domain-ge')
pb = (x * S.eye(6) - mb).det(method='domain-ge')
identity = h * (h + ell - ell * x * (x - h))
assert S.expand(pa - pb - identity) == 0
print('IDENTITY det(xI-M_A)-det(xI-M_B)=h*(h+ell-ell*x*(x-h))')

# Elimination: two subdivision vertices, a pendant P3, and a single leaf.
y0, y1 = S.symbols('y0 y1')
p1 = (t * y0 + y1) / (t**2 - 1)
p2 = (y0 + t * y1) / (t**2 - 1)
assert S.cancel(t * p1 - y0 - p2) == 0
assert S.cancel(t * p2 - p1 - y1) == 0
pmid, pend = t * y0 / (t**2 - 1), y0 / (t**2 - 1)
assert S.cancel(t * pmid - y0 - pend) == 0
assert S.cancel(t * pend - pmid) == 0
z = t**3 - (a + 3) * t
assert S.cancel(t * (t**2 - 1) - (a + 2) * t - z) == 0
assert S.cancel((t**2 - 1) / t - (t - 1 / t)) == 0
print('ELIMINATION all non-hub equations and z=t^3-(a+3)*t checked symbolically')

# Reflection of B: Perron vector (u,v,u,w,q,w).
zz = S.symbols('z')
reflection = (zz * S.eye(6) - mb) * S.Matrix([u, v, u, w, q, w])
assert reflection[0] == (zz - h) * u - v - w
assert reflection[1] == zz * v - 2 * u
assert reflection[3] == zz * w - u - q
# Combining those three zero equations gives z(z-h)=3+q/u>3.
assert S.expand(zz * reflection[0] + reflection[1] + reflection[3]
                - ((zz * (zz - h) - 3) * u - q)) == 0
assert S.cancel((h * (h - 2 * ell)).subs({h: t, ell: t - 1 / t}) - (2 - t**2)) == 0
print('REFLECTION zero equations imply (z*(z-h)-3)*u=q; threshold 2-t^2<0 for t>2')

c = S.Rational(213, 100)
exact_certificate('A', 0, True, c)
exact_certificate('B', 0, False, c)
print('STRICT rho(B_0)<213/100<rho(A_0), by positive integer residual certificates')
print('ALL_CHECKS_OK')

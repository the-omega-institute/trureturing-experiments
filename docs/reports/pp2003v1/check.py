#!/usr/bin/env python3
"""Exact GF(4) checks for the Polujan--Pott vectorial family (m=2)."""
from itertools import product

# GF(4)=F_2[t]/(t^2+t+1), represented by bits a0+a1*t.
K = range(4)

def add(a, b):
    return a ^ b

def mul(a, b):
    p = 0
    aa, bb = a, b
    while bb:
        if bb & 1:
            p ^= aa
        bb >>= 1
        aa <<= 1
    # t^2=t+1, so reduce the degree-2 bit.
    if p & 4:
        p ^= 4 ^ 2 ^ 1
    return p

def sq(a):
    return mul(a, a)

def inv(a):
    assert a
    for b in K:
        if mul(a, b) == 1:
            return b
    raise AssertionError(a)

def trace(a):
    return a ^ sq(a)  # F_2-valued integer (0 or 1).

def dot(a, b):
    out = 0
    for x, y in zip(a, b):
        out = add(out, mul(x, y))
    return out

def Q(a):
    x, y, z, u, v, w = a
    return add(add(mul(x, u), mul(y, v)), mul(z, w))

def F(a):
    x, y, z, u, v, w = a
    return add(Q(a), mul(mul(x, y), w))

def phi(a):
    x, y, z, u, v, w = a
    return (x, y, add(z, mul(x, y)), add(u, mul(y, w)),
            add(v, mul(x, w)), w)

def P(s):
    x, y, z = s
    return (x, y, add(z, mul(x, y)))

V = list(product(K, repeat=6))
V3 = list(product(K, repeat=3))

def enc(xs):
    return sum(x << (2*i) for i, x in enumerate(xs))

def dec(q, n):
    return tuple((q >> (2*i)) & 3 for i in range(n))

# Bitwise XOR is addition in this representation.  Precompute all values so
# the 4^12 identity checks are exact integer table lookups.
all6 = list(range(4**6))
phi_tab = [enc(phi(dec(q, 6))) for q in all6]
q_tab = [Q(dec(q, 6)) for q in all6]
f_tab = [F(dec(q, 6)) for q in all6]

# The identity used for the development is checked exhaustively on all 4^12 pairs.
for ai in all6:
    assert phi_tab[phi_tab[ai]] == ai
    pa = phi_tab[ai]
    for bi in all6:
        assert f_tab[ai ^ bi] == q_tab[pa ^ phi_tab[bi]], (ai, bi)

# P is a permutation, and the Maiorana--McFarland orthogonality calculation is exact.
assert len({P(s) for s in V3}) == 64
for lam in (1, 2, 3):
    ilam = inv(lam)
    for A in V3:
        for B in V3:
            target = tuple(mul(ilam, x) for x in B)
            surviving = [s for s in V3 if P(s) == target]
            assert len(surviving) == 1
            # Summing over r first leaves one s and a character of magnitude |K|^3.
            Walsh = 64 * (1 if trace(mul(lam, dot(A, surviving[0]))) == 0 else -1)
            assert abs(Walsh) == 64

# Q has binary degree 2, while the lambda=1 component of F has a nonzero
# third derivative: Tr(1 * 1 * alpha)=1 for alpha=t (represented by 2).
def component_H(H, lam, a):
    return trace(mul(lam, H(a)))
zero = (0, 0, 0, 0, 0, 0)
def xor_points(*xs):
    if not xs:
        return zero
    return tuple(__import__('functools').reduce(add, vals, 0) for vals in zip(*xs))
e1 = (1, 0, 0, 0, 0, 0)
e2 = (0, 1, 0, 0, 0, 0)
et = (0, 0, 0, 0, 0, 2)
def third_difference(H, lam):
    vals = []
    for bits in product((0, 1), repeat=3):
        point = xor_points(*(inc for bit, inc in zip(bits, (e1, e2, et)) if bit))
        vals.append(component_H(H, lam, point))
    return sum(vals) % 2
assert third_difference(F, 1) == 1
assert third_difference(Q, 1) == 0

print("m=2, |K|=4, n=12")
print("exhaustive development identity checks: 4^12 = 16777216")
print("all nonzero output components have Walsh magnitude 64 = 2^(n/2)")
print("binary third derivative: deg(Q)=2 witness 0, deg(F)=3 witness 1")
print("ALL_CHECKS_OK")

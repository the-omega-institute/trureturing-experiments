#!/usr/bin/env python3
"""Exact checks for the Rajchel--Gasiorowski--Zyczkowski Question 3 draft.
All arithmetic is integer, Fraction, or SymPy exact algebra; no floats are used.
"""
from fractions import Fraction
from itertools import combinations
import sympy as sp

MOD = (1 << 5) | (1 << 2) | 1  # X^5 + X^2 + 1, irreducible over F_2
Q = 1 << 5


def gf_add(a, b):
    return a ^ b


def gf_mul(a, b):
    out = 0
    aa, bb = a, b
    while bb:
        if bb & 1:
            out ^= aa
        bb >>= 1
        aa <<= 1
        if aa & Q:
            aa ^= MOD
    return out


def gf_pow(a, k):
    out = 1
    while k:
        if k & 1:
            out = gf_mul(out, a)
        a = gf_mul(a, a)
        k >>= 1
    return out


def gf_inv(a):
    assert a != 0
    return gf_pow(a, Q - 2)


def gf_cube(a):
    return gf_mul(gf_mul(a, a), a)


def hadamard_skew(n):
    if n == 1:
        return [[1]]
    h = hadamard_skew(n // 2)
    ht = [list(row) for row in zip(*h)]
    return [h[i] + h[i] for i in range(n // 2)] + [
        [-ht[i][j] for j in range(n // 2)] + ht[i] for i in range(n // 2)
    ]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def check_robust_hadamard():
    for k in range(0, 8):
        n = 1 << k
        h = hadamard_skew(n)
        hh = matmul(h, [list(row) for row in zip(*h)])
        assert hh == [[n if i == j else 0 for j in range(n)] for i in range(n)]
        assert all(h[i][j] + h[j][i] == (2 if i == j else 0) for i in range(n) for j in range(n))
    return "skew Hadamard recursion verified at orders 1..128"


def check_unistochastic_triangle_fidelity():
    # Source Lemma 4.4: a strong-complementarity triangle at n=4.
    # Q swaps (0,1) and (2,3); H is the exact skew-Hadamard recursion at order 4.
    h = hadamard_skew(4)
    a, b, c = Fraction(1, 2), Fraction(1, 4), Fraction(1, 8)
    assert a + b + 2 * c == 1
    q = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
    B = [[(b if i == j else a if q[i][j] else c) for j in range(4)] for i in range(4)]
    assert all(sum(row) == 1 for row in B)
    assert all(sum(B[i][j] for i in range(4)) == 1 for j in range(4))
    # SymPy keeps sqrt(Fraction) exact. U = entrywise sqrt(B) times normalized H.
    U = sp.Matrix([[sp.sqrt(sp.Rational(B[i][j].numerator, B[i][j].denominator)) * sp.Integer(h[i][j])
                    for j in range(4)] for i in range(4)])
    assert sp.simplify(U * U.T - sp.eye(4)) == sp.zeros(4)
    assert all(sp.simplify(U[i, j] ** 2 - sp.Rational(B[i][j].numerator, B[i][j].denominator)) == 0
               for i in range(4) for j in range(4))
    return "n=4 strongly-complementary triangle has an exact orthogonal lift"


def check_sidon_and_clifford_bound():
    # S={(x,x^3)} in (F_32)^2. Verify all inverses and unique unordered pair sums.
    nonzero = range(1, Q)
    assert all(gf_mul(x, gf_inv(x)) == 1 for x in nonzero)
    S = [(x, gf_cube(x)) for x in range(Q)]
    assert len(set(S)) == Q and (0, 0) in S
    multiplicities = {}
    for u, v in combinations(S, 2):
        key = (gf_add(u[0], v[0]), gf_add(u[1], v[1]))
        assert key != (0, 0)
        multiplicities[key] = multiplicities.get(key, 0) + 1
    assert max(multiplicities.values()) == 1
    assert len(multiplicities) == Q * (Q - 1) // 2
    m = 5
    r = 4 * m + 2
    q_minus_one = Q - 1
    matrix_dim = (Q * Q) ** 2
    word_count = 1 << r
    assert q_minus_one >= r
    assert word_count > matrix_dim
    return (f"F_32 cubic Sidon support: {len(multiplicities)} distinct nonzero pair sums; "
            f"q-1={q_minus_one} >= r={r}, words={word_count} > dim M_n(C)={matrix_dim}")


if __name__ == "__main__":
    messages = [check_unistochastic_triangle_fidelity(), check_sidon_and_clifford_bound(), check_robust_hadamard()]
    for msg in messages:
        print(msg)
    print("ALL_CHECKS_OK")

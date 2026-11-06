#!/usr/bin/env python3
"""Exact finite-field checks for arXiv:2611.04863v1, Section 8.

Only Python integers are used. Elements of GF(2^d) are bit polynomials,
reduced modulo an irreducible polynomial checked below by trial division.
The lower-bound test covers a whole vector space by testing a basis of
linear moment functionals; it does not enumerate or sample multipliers.
"""

from functools import reduce
from itertools import permutations
from math import comb
from operator import xor


def psum(values):
    return reduce(xor, values, 0)


def remainder(a, b):
    assert b > 0
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


class Field:
    def __init__(self, degree, modulus):
        self.degree, self.q, self.modulus = degree, 1 << degree, modulus
        assert modulus.bit_length() == degree + 1
        for d in range(1, degree // 2 + 1):
            for divisor in range(1 << d, 1 << (d + 1)):
                assert remainder(modulus, divisor) != 0
        self.table = [[self.raw_mul(a, b) for b in range(self.q)]
                      for a in range(self.q)]
        for a in range(1, self.q):
            assert self.pow(a, self.q - 1) == 1
        for a in range(self.q):
            assert self.pow(a, self.q) == a
        print(f"FIELD GF(2^{degree}) modulus={modulus:#x} irreducible=OK")

    def raw_mul(self, a, b):
        answer = 0
        while b:
            if b & 1:
                answer ^= a
            b >>= 1
            a <<= 1
            if a & self.q:
                a ^= self.modulus
        return answer

    def mul(self, a, b):
        return self.table[a][b]

    def prod(self, values):
        return reduce(self.mul, values, 1)

    def pow(self, a, exponent):
        assert exponent >= 0
        result = 1
        while exponent:
            if exponent & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            exponent >>= 1
        return result


def det3(field, matrix):
    # All permutation signs coincide in characteristic two.
    return psum(field.prod(matrix[i][p[i]] for i in range(3))
                for p in permutations(range(3)))


def gram(field, points, values):
    columns = [(1, x, y) for x, y in points]
    return [[psum(field.prod((field.pow(f, 2), v[i], v[j]))
                  for v, f in zip(columns, values))
             for j in range(3)] for i in range(3)]


def basis(r, m):
    return [(i, j) for j in range(r) for i in range(m // r + 1)
            if r * i + (r + 1) * j <= m]


def points_on_curve(field, r):
    return [(x, y) for x in range(field.q) for y in range(field.q)
            if field.pow(y, r) ^ y == field.pow(x, r + 1)]


def powers(field, value, top):
    result = [1]
    for _ in range(top):
        result.append(field.mul(result[-1], value))
    return result


def lower_bound(field, r, points, m0):
    """Every f in L((m0-1)Q) has an alternating scaled 3x3 Gram.

    For u=1,x,y, the diagonal is (sum_D f*u)^2. Vanishing of
    sum_D monomial*u on the complete basis is equivalent to its
    vanishing for every coefficient tuple. Nestedness covers all m<m0.
    The 3x3 symmetric zero-diagonal determinant is identically zero.
    """
    rows = basis(r, m0 - 1)
    xp = [powers(field, x, max(i for i, _ in rows) + 1) for x, _ in points]
    yp = [powers(field, y, r) for _, y in points]
    for i, j in rows:
        for di, dj in [(0, 0), (1, 0), (0, 1)]:
            assert psum(field.mul(xx[i + di], yy[j + dj])
                        for xx, yy in zip(xp, yp)) == 0
    # Independent formal-polynomial calculation of the determinant:
    # the two abc terms cancel and every other term has a zero diagonal.
    terms = {}
    entry = {(0, 1): 'a', (1, 0): 'a', (0, 2): 'b',
             (2, 0): 'b', (1, 2): 'c', (2, 1): 'c'}
    for p in permutations(range(3)):
        if any(i == p[i] for i in range(3)):
            continue
        monomial = ''.join(sorted(entry[i, p[i]] for i in range(3)))
        terms[monomial] = terms.get(monomial, 0) ^ 1
    assert not any(terms.values())
    print(f"LOWER r={r} all_m<{m0} basis_size={len(rows)} "
          f"linear_moments={3 * len(rows)} zero; det_alternating_3=0")


def source_elliptic():
    field, r = Field(2, 0b111), 2
    points = points_on_curve(field, r)
    assert len(points) == 8
    values = [1 ^ field.pow(x, 3) ^ field.mul(x, y) for x, y in points]
    assert all(values)
    matrix = gram(field, points, values)
    assert matrix == [[0, 1, 1], [1, 0, 0], [1, 0, 1]]
    assert det3(field, matrix) == 1
    assert max(0, 2 * 3, 2 + 3) == 6
    lower_bound(field, r, points, 6)
    print(f"FIDELITY Theorem_5.2 n=8 multiplier=1+x^3+xy "
          f"pole_order=6 Gram={matrix} determinant=1 threshold=6")


def check_family(e, modulus):
    field, r = Field(2 * e, modulus), 1 << e
    points = points_on_curve(field, r)
    n, g = r ** 3, r * (r - 1) // 2
    m0 = r ** 3 + r ** 2 - 2 * r - 2
    s = n + 2 * g - 2
    assert len(points) == n
    assert all(sum(x == u for x, _ in points) == r for u in range(field.q))
    assert basis(r, r + 1) == [(0, 0), (1, 0), (0, 1)]
    assert m0 == n + 2 * g - r - 2
    assert n + 2 * g - 3 - m0 == r - 1
    assert gram(field, points, [1] * n) == [[0] * 3 for _ in range(3)]
    lower_bound(field, r, points, m0)
    # Independent direct check of the source's residue-sum identity.
    sr = basis(r, s)
    xp = {x: powers(field, x, max(i for i, _ in sr)) for x, _ in points}
    yp = {y: powers(field, y, r - 1) for _, y in points}
    for i, j in sr:
        assert psum(field.mul(xp[x][i], yp[y][j]) for x, y in points) == 0

    subfield = [a for a in range(field.q) if field.pow(a, r) == a]
    betas = [b for b in range(field.q) if field.pow(b, r) ^ b == 1]
    assert len(subfield) == len(betas) == r
    a = next(a for a in subfield if a not in (0, 1))
    beta = betas[0]
    special = [(0, 0), (1, beta), (a, field.mul(field.pow(a, 2), beta))]
    weights = [field.mul(beta, a ^ 1), field.mul(beta, a), beta]
    assert len(set(special)) == 3 and all(p in points for p in special)
    assert psum(weights) == 0
    assert psum(field.mul(w, x) for w, (x, _) in zip(weights, special)) == 0
    ymoment = psum(field.mul(w, y) for w, (_, y) in zip(weights, special))
    assert ymoment == field.prod((field.pow(beta, 2), a, a ^ 1)) != 0

    coefficients = {(0, 0): 1}
    def add(i, j, c):
        coefficients[i, j] = coefficients.get((i, j), 0) ^ c
    for (u, v), w in zip(special, weights):
        # Expand (1+(x-u)^(q-1))(1+(y-v)^(r-1)) exactly.
        aa = [field.pow(u, field.q - 1 - i) if comb(field.q - 1, i) % 2 else 0
              for i in range(field.q)]
        bb = [field.pow(v, r - 1 - j) if comb(r - 1, j) % 2 else 0
              for j in range(r)]
        aa[0] ^= 1
        bb[0] ^= 1
        for i, c in enumerate(aa):
            for j, d in enumerate(bb):
                add(i, j, field.prod((w, c, d)))
        indicator = [field.mul(1 ^ field.pow(x ^ u, field.q - 1),
                               1 ^ field.pow(y ^ v, r - 1)) for x, y in points]
        assert indicator == [int(p == (u, v)) for p in points]
    coefficients = {ij: c for ij, c in coefficients.items() if c}
    assert max(r * i + (r + 1) * j for i, j in coefficients) == m0
    assert coefficients.get((field.q - 1, r - 1), 0) == 0
    assert coefficients.get((field.q - 2, r - 1), 0) == 0
    assert coefficients[field.q - 1, r - 2] == ymoment
    values = [1 ^ psum(w for p, w in zip(special, weights) if p == (x, y))
              for x, y in points]
    assert all(values)
    for (x, y), value in zip(points, values):
        assert psum(field.prod((c, xp[x][i], yp[y][j]))
                    for (i, j), c in coefficients.items()) == value
    matrix = gram(field, points, values)
    expected_det = field.prod((field.pow(beta, 8), field.pow(a, 4),
                               field.pow(a ^ 1, 4)))
    assert det3(field, matrix) == expected_det != 0
    information_matrix = [[1, 1, 1], [p[0] for p in special], [p[1] for p in special]]
    assert det3(field, information_matrix) == field.prod((beta, a, a ^ 1)) != 0
    print(f"UPPER r={r} n={n} g={g} a={a} beta={beta} "
          f"weights={weights} support={len(coefficients)} pole_order={m0} "
          f"nonzero_evaluations={len(values)} Gram={matrix} determinant={expected_det}")
    print(f"THRESHOLD r={r} exact={m0}; source_Theorem_5.1_bounds="
          f"[{m0},{n + 2 * g - 3}]; residue_basis_checked={len(sr)}")


if __name__ == '__main__':
    source_elliptic()
    check_family(2, 0b10011)
    check_family(3, 0b1000011)
    print('ALL_CHECKS_OK')

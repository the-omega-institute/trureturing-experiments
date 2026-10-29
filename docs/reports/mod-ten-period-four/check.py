#!/usr/bin/env python3
"""Exact finite checks for Hanna's first OEIS A357398 COMMENTS conjecture.

The independent source equation is solved over Fraction, without the Stirling
substitution. A second construction solves the integral ordinary-series equation
of A357397 and applies the factorial-weighted Stirling transform. No finite
check establishes the all-index conjecture.
"""
from fractions import Fraction
from math import comb, factorial

N = 32
DIRECT_N = 16
# OEIS A357398 revision 13, %S/%T/%U, offset 0.
SOURCE = [
    1, 1, 3, 37, 1083, 53701, 3934443, 395502997, 51998075643,
    8643190760101, 1770707733052683, 438247927304923957,
    128926370847248904603, 44477192002157773868101,
    17787359321176954021105323, 8164879801560415752441320917,
    4264616618784184682736855291963,
]
# OEIS A357397 revision 13, %S/%T/%U, offset 0.
SIBLING = [
    1, 1, 1, 5, 37, 367, 4463, 63797, 1043961, 19208815,
    392278493, 8802891869, 215335062049, 5704017709585,
    162695460126735, 4972552233126827, 162156046298476305,
    5620675413587870585, 206382551428754263839, 8003189847508668434429,
]


def mul(a, b, d):
    c = [0] * (d + 1)
    for i, x in enumerate(a[:d + 1]):
        if x:
            for j, y in enumerate(b[:d + 1 - i]):
                if y:
                    c[i + j] += x * y
    return c


def power(a, r, d):
    out = [1] + [0] * d
    while r:
        if r % 2:
            out = mul(out, a, d)
        a = mul(a, a, d)
        r //= 2
    return out


def ordinary_summand(f, r, d):
    inv = [(-1)**j * comb(r + j - 1, j) for j in range(d + 1)]
    quotient = mul(f, inv, d)
    bracket = [-v for v in quotient]
    bracket[0] += 1
    assert bracket[0] == 0
    return power(bracket, r, d)


def ordinary_residual(f, d):
    return [sum(ordinary_summand(f, r, d)[j] for r in range(1, d + 1))
            for j in range(d + 1)]


def construct_integral(d):
    f = [1]
    for j in range(1, d + 1):
        f.append(0)
        f[j] = sum(ordinary_summand(f, r, j)[j] for r in range(1, j + 1))
        assert isinstance(f[j], int)
    return f


def exp_scaled(q, d):
    return [Fraction(q**j, factorial(j)) for j in range(d + 1)]


def source_residual(a, d):
    total = [0] * (d + 1)
    for r in range(1, d + 1):
        exp_rx = exp_scaled(r, d)
        bracket = [exp_rx[j] - a[j] for j in range(d + 1)]
        assert bracket[0] == 0
        term = mul(exp_scaled(-r*r, d), power(bracket, r, d), d)
        total = [x + y for x, y in zip(total, term)]
    return total


def construct_source(d):
    a = [Fraction(1)]
    for j in range(1, d + 1):
        a.append(Fraction(0))
        a[j] = source_residual(a, j)[j]
    assert source_residual(a, d) == [0] * (d + 1)
    return a


def stirling_rows(d):
    rows = [[1]]
    for n in range(1, d + 1):
        old = rows[-1]
        rows.append([0] + [(old[k - 1] if k - 1 < len(old) else 0)
                              + k * (old[k] if k < len(old) else 0)
                          for k in range(1, n + 1)])
    return rows


def main():
    f = construct_integral(N)
    assert ordinary_residual(f, N) == [0] * (N + 1)
    assert f[:len(SIBLING)] == SIBLING
    assert f[:5] == [1, 1, 1, 5, 37]
    jet = [1, 1, 1, 5, 37]
    jet_terms = [ordinary_summand(jet, r, 4) for r in range(1, 5)]
    assert jet_terms == [
        [0, 0, -1, -4, -33], [0, 0, 1, -4, 0],
        [0, 0, 0, 8, -48], [0, 0, 0, 0, 81],
    ]
    s = stirling_rows(N)
    a = [sum(f[k] * factorial(k) * s[n][k] for k in range(n + 1))
         for n in range(N + 1)]
    assert a[:len(SOURCE)] == SOURCE
    direct = construct_source(DIRECT_N)
    direct_scaled = [factorial(n) * direct[n] for n in range(DIRECT_N + 1)]
    assert direct_scaled == a[:DIRECT_N + 1]
    assert all(x.denominator == 1 for x in direct_scaled)
    # Independently compose F(exp(x)-1) in the rational series ring.
    e_minus_one = exp_scaled(1, DIRECT_N)
    e_minus_one[0] = 0
    composed = [sum(f[k] * power(e_minus_one, k, DIRECT_N)[j]
                    for k in range(DIRECT_N + 1))
                for j in range(DIRECT_N + 1)]
    assert composed == direct
    # Polynomial in y=exp(x): 1+(y-1)+(y-1)^2+2(y-1)^4 modulo 5.
    poly = [sum(c * comb(k, j) * (-1)**(k-j)
                for k, c in [(0, 1), (1, 1), (2, 1), (4, 2)] if j <= k) % 5
            for j in range(5)]
    assert poly == [3, 1, 3, 2, 2]
    expected = [3, 1, 3, 7]
    for n in range(1, N + 1):
        assert a[n] % 2 == 1
        assert a[n] % 5 == (1 + 3*2**n + 2*3**n + 2*4**n) % 5
        assert a[n] % 10 == expected[n % 4]
    assert [a[n] % 10 for n in range(1, 9)] == [1, 3, 7, 3, 1, 3, 7, 3]
    print('A357397 fidelity: coefficients 0..19 match revision 13')
    print('A357398 fidelity: coefficients 0..16 match revision 13')
    print('F implicit equation: exact zero residual through degree', N)
    print('A original NAME equation: exact zero residual through degree', DIRECT_N)
    print('Independent rational construction = F(exp(x)-1): degrees 0..16')
    print('F[0..4] =', f[:5])
    print('Degree-four summands r=1..4 =', jet_terms)
    print('A factorial coefficients a[0..8] =', a[:9])
    print('Polynomial coefficients modulo 5 =', poly)
    print('Oddness, exponential residue formula, and last digits: indices 1..32')
    print('a[1..16] modulo 10 =', [a[n] % 10 for n in range(1, 17)])
    print('Evidence scope: finite exact arithmetic; no all-index or Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact checks for Konoike's Question 5.5; no numerical root finding.

Polynomials are lists of fractions in increasing monomial order. All
threshold checks are finite evidence, separate from the uniform paper argument.
"""
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial


def add(p, q):
    out = [Q(0)] * max(len(p), len(q))
    for i, c in enumerate(p):
        out[i] += c
    for i, c in enumerate(q):
        out[i] += c
    return out


def mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def scale(p, a):
    return [a * c for c in p]


def evaluate(p, x):
    y = Q(0)
    for c in reversed(p):
        y = y * x + c
    return y


def shifted_choose(r):
    """binom(X+r,r) = product_{i=1}^r (X+i)/r!, also r=0."""
    p = [Q(1)]
    for i in range(1, r + 1):
        p = mul(p, [Q(i), Q(1)])
    return scale(p, Q(1, factorial(r)))


def s_polynomial(n, k):
    m = n - k
    out = [Q(0)]
    for j in range(k):
        out = add(out, scale(shifted_choose(j), Q(comb(m - 1 + j, j))))
    return out


def enumerator(n, k):
    assert 1 <= k <= n // 2
    m = n - k
    return scale(mul(shifted_choose(m), s_polynomial(n, k)),
                 Q(1, comb(n - 1, k - 1)))


def dilate(p, s):
    return [c * s ** i for i, c in enumerate(p)]


def magic_coefficients(p):
    """Triangular inversion of c_j = sum_{i<=j} a_i binom(d-i,j-i)."""
    d = len(p) - 1
    a = []
    for i, c in enumerate(p):
        a.append(c - sum(a[j] * comb(d - j, i - j) for j in range(i)))
    return a


def reconstruct(a):
    """Independent forward expansion of sum a_i X^i (X+1)^(d-i)."""
    d = len(a) - 1
    return [sum(a[i] * comb(d - i, j - i) for i in range(j + 1))
            for j in range(d + 1)]


def check_magic(p, s):
    q = dilate(p, s)
    a = magic_coefficients(q)
    assert reconstruct(a) == q
    return all(c >= 0 for c in a)


def main():
    # Konoike, Proposition 5.4: D_{2,4} is not magic positive.
    # Its displayed expansion has an erroneous final +X^3; the enumerator
    # gives final coefficient 0, and the negative coefficient -5/6 is unchanged.
    p = enumerator(4, 2)
    target = [Q(1), Q(-5, 6), Q(1, 6), Q(0)]
    assert magic_coefficients(p) == target
    assert reconstruct(target) == p
    assert target[1] < 0
    print('FIDELITY Proposition 5.4: D_{2,4}(X) is not magic positive.')
    print('D_{2,4}(X) exact magic coefficients:', ', '.join(map(str, target)))
    print('SOURCE TYPO: the displayed +X^3 in the rank-two proof must be omitted; D_{2,4}(-1)=0.')

    # Independent lattice count for Ferroni's literal polytope inequalities:
    # 0<=x_i<=t, sum x_i=tk, sum_{i=k+1}^n x_i<=t.
    lattice_cases = 0
    for n in range(2, 8):
        for k in range(1, n // 2 + 1):
            p = enumerator(n, k)
            for t in range(4):
                count = sum(1 for x in product(range(t + 1), repeat=n)
                            if sum(x) == t * k and sum(x[k:]) <= t)
                assert evaluate(p, Q(t)) == count, (n, k, t)
                lattice_cases += 1
    print('FIDELITY exact lattice counts: n=2..7, 1<=k<=floor(n/2), t=0..3;',
          lattice_cases, 'cases')

    # Exhaust every candidate smaller integer, together with the threshold and
    # two following integers. This is a finite check, not a uniform proof.
    pairs = dilations = known_pairs = 0
    for n in range(2, 29):
        for k in range(1, n // 2 + 1):
            m = n - k
            p = enumerator(n, k)
            successes = []
            for s in range(1, m + 3):
                ok = check_magic(p, Q(s))
                assert ok == (s >= m), (n, k, s)
                if ok:
                    successes.append(s)
                dilations += 1
            assert min(successes) == m
            pairs += 1
            if k <= 2:
                known_pairs += 1
    print('FIDELITY Proposition 5.4 ranks 1 and 2 through n=28;',
          known_pairs, 'parameter pairs')
    print('THRESHOLD n=2..28, every admissible k and integer s=1..n-k+2;',
          pairs, 'parameter pairs;', dilations, 'dilations')

    # Larger exact check of the key identity, including balanced n=2k and k=1.
    pairs = evaluations = rational_dilations = 0
    for n in range(2, 61):
        for k in range(1, n // 2 + 1):
            m = n - k
            s_poly = s_polynomial(n, k)
            assert len(s_poly) - 1 == k - 1
            assert s_poly[-1] == Q(comb(m + k - 2, k - 1), factorial(k - 1))
            for r in range(1, k + 1):
                assert evaluate(s_poly, Q(-r)) == (-1) ** (r - 1) * comb(m - 1, r - 1)
                evaluations += 1
            p = enumerator(n, k)
            assert len(p) - 1 == n - 1 and p[-1] > 0
            assert evaluate(p, Q(0)) == 1
            assert evaluate(p, Q(1)) == k * m + 1
            assert evaluate(p, Q(-m)) == 0
            # Test S directly as well as D, so Q cannot conceal a failed S sign.
            assert check_magic(s_poly, Q(m))
            for s in [Q(2 * m - 1, 2), Q(m), Q(2 * m + 1, 2)]:
                assert check_magic(p, s) == (s >= m), (n, k, s)
                rational_dilations += 1
            pairs += 1
    print('IDENTITY S(-r)=(-1)^(r-1) binom(m-1,r-1): n=2..60;',
          pairs, 'parameter pairs;', evaluations, 'evaluations')
    print('DEGREE, leading coefficient, D(0), D(1), D(-m), and S(mX) checked;',
          pairs, 'parameter pairs')
    print('RATIONAL THRESHOLD s=m-1/2,m,m+1/2;', rational_dilations, 'dilations')
    print('Arithmetic: Python integers and fractions.Fraction only; no floating point.')
    print('Finite checks support the proposed uniform argument; no Lean proof is claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

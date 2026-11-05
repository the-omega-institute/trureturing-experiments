#!/usr/bin/env python3
"""Exact checks for Au, arXiv:2412.17196v4, Conjecture 5.6.

Fidelity: Example 5.5's polynomial reversal and skipped pole classes,
plus the explicit g_n recurrence in arXiv:2412.17196v1, section 4.
Settlement diagnostics: the proposed denominator/parity identity on finite
ranges, including every conjecture parameter in -299 <= k <= 1.
Only Python integers and fractions.Fraction are used. Sampling is not a proof.
"""
from fractions import Fraction as Q
from math import comb, factorial, prod


def rising(a, n):
    return prod((a + r for r in range(n)), start=Q(1))


def gauss_terminating(a, m, c, z):
    terms = []
    for j in range(m + 1):
        den = rising(c, j) * factorial(j)
        assert den != 0, (a, m, c, j)
        terms.append(rising(a, j) * rising(Q(-m), j) * z**j / den)
    return sum(terms, Q(0))


def cleared(m, k):
    a = [3*k + 10*r for r in range(m)]
    b = [3*k + 5*r for r in range(m)]
    assert all(x != 0 for x in b)
    d = 8**m * prod(b)
    terms = [(-1)**j * comb(m, j) * 3**j * 8**(m-j)
             * prod(a[:j]) * prod(b[j:]) for j in range(m+1)]
    n = sum(terms)
    assert all(t % 8 == 0 for t in terms[:-1])
    assert n % 8 == ((-3)**m * prod(a)) % 8
    assert n % 2 == 1
    assert d != 0 and n != 0
    return d, n


def check_pair(m, k):
    assert k % 2 == 1 and k % 5 != 0
    h = gauss_terminating(Q(3*k, 10), m, Q(3*k, 5), Q(3, 4))
    d, n = cleared(m, k)
    assert d * h == n, (m, k, h, d, n)
    assert h != 0
    return h


def multiply(p, q, degree):
    out = [Q(0)] * min(len(p) + len(q) - 1, degree + 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            if i+j < len(out):
                out[i+j] += a*b
    return out


def coefficient_g(n):
    """[x^n] prod(1+alpha*x)^((n-1)/5), from v1's definition.

    Compute by the generalized binomial theorem and convolution,
    independently of both the published recurrence and H_{m,k}.
    """
    exponent = Q(n-1, 5)
    series = [Q(1)]
    for alpha in (Q(1), Q(2), Q(3), Q(3, 2)):
        factor = [Q(1)]
        for j in range(1, n+1):
            factor.append(factor[-1] * (exponent-j+1) * alpha / j)
        series = multiply(series, factor, n)
    return series[n]


def main():
    f = [Q(1)]
    for factor in ([1, 1], [1, 2], [1, 3], [2, 3]):
        f = multiply(f, factor, 4)
    assert list(reversed(f)) == [9*a / 3**i for i, a in enumerate(f)]
    print('FIDELITY v4 Example 5.5: x^4 f(1/x) = 9 f(x/3); f =', f)

    g = [coefficient_g(n) for n in range(61)]
    assert g[3] == g[13] == 0
    skipped = [n for n in range(61) if n % 10 in (3, 5, 7, 9)]
    assert all(g[n] == 0 for n in skipped)
    for n in range(41):
        residual = (729*(n+4)*(n+9)*(n+14)*(n+19)*(6*n+103)*g[n]
                    - 1875*(n+14)*(n+19)
                    *(216*n**3+7596*n**2+83694*n+290603)*g[n+10]
                    + 50000*(n+12)*(n+16)*(n+18)*(n+20)*(6*n+43)*g[n+20])
        assert residual == 0, (n, residual)
    print('FIDELITY v1 section 4: g_0,...,g_10 =', ', '.join(map(str, g[:11])))
    print('FIDELITY: 41 published recurrence identities; skipped classes:', len(skipped))

    pairs = [(m, k) for m in range(25) for k in range(-49, 50)
             if k % 2 == 1 and k % 5 != 0]
    for m, k in pairs:
        check_pair(m, k)
    print('PARITY IDENTITY: m=0,...,24; odd k=-49,...,49 with 5 not dividing k;', len(pairs), 'pairs')

    ks = [k for k in range(-299, 2) if k % 10 in (1, 3, 7, 9)]
    for k in ks:
        m = (1-k)//2
        assert 2*m == 1-k
        assert Q(5)*Q(k, 5)-1 == -2*m
        check_pair(m, k)
    print('CONJECTURE 5.6: all', len(ks), 'admissible integers -299 <= k <= 1 checked exactly')
    for k in (1, -1, -3, -7, -9):
        m = (1-k)//2
        print('EXACT H: k =', k, 'm =', m, 'H =', check_pair(m, k))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

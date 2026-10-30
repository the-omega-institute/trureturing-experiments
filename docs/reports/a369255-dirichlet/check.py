#!/usr/bin/env python3
"""Exact arithmetic checks for OEIS A369255's two convolution conjectures.

Only Python integer arithmetic is used. Source definitions are evaluated by
finite divisor enumeration, not by the claimed exponent classification.
OEIS data prefixes are attributed to their entries under the OEIS EULA.
"""
from itertools import combinations_with_replacement
from math import comb, gcd, isqrt, prod

LIMIT = 2**20


def convolution(f, g):
    """Return sum_{d|n} f(d)g(n/d) for 1 <= n <= LIMIT."""
    ans = [0] * (LIMIT + 1)
    for d in range(1, LIMIT + 1):
        if f[d]:
            fd = f[d]
            for k in range(1, LIMIT // d + 1):
                ans[d * k] += fd * g[k]
    return ans


def prime_data():
    spf = [0] * (LIMIT + 1)
    primes = []
    mu = [0] * (LIMIT + 1)
    mu[1] = 1
    for n in range(2, LIMIT + 1):
        if spf[n] == 0:
            spf[n] = n
            primes.append(n)
            mu[n] = -1
        for p in primes:
            if p * n > LIMIT:
                break
            spf[p * n] = p
            if n % p == 0:
                mu[p * n] = 0
                break
            mu[p * n] = -mu[n]
    return spf, mu


def exponents(n, spf):
    out = []
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        out.append(e)
    return out


def classification(es):
    return int(all(e % 4 in (0, 1) for e in es)
               and sum(e % 4 == 1 for e in es) % 2 == 0)


def main():
    spf, mu = prime_data()
    # Literal A038548: each d contributes to multiples n >= d*d.
    g = [0] * (LIMIT + 1)
    # Literal A007875: enumerate the unordered factor pair (d, n//d).
    a = [0] * (LIMIT + 1)
    for d in range(1, isqrt(LIMIT) + 1):
        for n in range(d * d, LIMIT + 1, d):
            g[n] += 1
            a[n] += int(gcd(d, n // d) == 1)
    h = convolution([0] + [1] * LIMIT, g)
    # Literal A219009's published PARI formula: sum over fourth-power divisors.
    b = [0] * (LIMIT + 1)
    d = 1
    while d**4 <= LIMIT:
        d4 = d**4
        for m in range(1, LIMIT // d4 + 1):
            b[d4 * m] += mu[m]
        d += 1
    # Literal A355837: multiplicative, with its published prime-power values.
    c = [0] * (LIMIT + 1)
    for n in range(1, LIMIT + 1):
        es = exponents(n, spf)
        c[n] = prod(2 * (e % 2) * (-1)**((e + 1) // 2) for e in es)
        assert g[n] == (prod(e + 1 for e in es) + 1) // 2, ('g', n)
        assert a[n] == (1 if n == 1 else 2**(len(es) - 1)), ('A', n)
        assert b[n] == prod((-1)**e if e % 4 < 2 else 0 for e in es), ('B', n)

    # Fidelity: the actual %S/%T/%U published prefixes, read from the entries.
    source_prefixes = {'A038548': [1, 1, 1, 2, 1, 2, 1, 2, 2, 2, 1, 3, 1, 2, 2, 3, 1, 3, 1, 3, 2, 2, 1, 4, 2, 2, 2, 3, 1, 4, 1, 3, 2, 2, 2, 5, 1, 2, 2, 4, 1, 4, 1, 3, 3, 2, 1, 5, 2, 3, 2, 3, 1, 4, 2, 4, 2, 2, 1, 6, 1, 2, 3, 4, 2, 4, 1, 3, 2, 4, 1, 6, 1, 2, 3, 3, 2, 4, 1, 5, 3, 2, 1, 6, 2, 2, 2, 4, 1, 6, 2, 3, 2, 2, 2, 6, 1, 3, 3, 5, 1, 4, 1, 4, 4], 'A140773': [1, 2, 2, 4, 2, 5, 2, 6, 4, 5, 2, 10, 2, 5, 5, 9, 2, 10, 2, 10, 5, 5, 2, 16, 4, 5, 6, 10, 2, 14, 2, 12, 5, 5, 5, 20, 2, 5, 5, 16, 2, 14, 2, 10, 10, 5, 2, 24, 4, 10, 5, 10, 2, 16, 5, 16, 5, 5, 2, 28, 2, 5, 10, 16, 5, 14, 2, 10, 5, 14, 2, 32, 2, 5, 10, 10, 5, 14, 2, 24, 9, 5, 2, 28, 5, 5, 5, 16, 2, 28, 5], 'A007875': [1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 2, 1, 1, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 4, 1, 1, 2, 2, 2, 2, 1, 2, 2, 2, 1, 4, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 4, 1, 2, 2, 1, 2, 4, 1, 2, 2, 4, 1, 2, 1, 2, 2, 2, 2, 4, 1, 2, 1, 2, 1, 4, 2, 2, 2, 2, 1, 4], 'A219009': [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0, -1, 1, 1, 1, -1, 0, -1, 0, 1, 1, -1, 0, 0, 1, 0, 0, -1, -1, -1, -1, 1, 1, 1, 0, -1, 1, 1, 0, -1, -1, -1, 0, 0, 1, -1, -1, 0, 0, 1, 0, -1, 0, 1, 0, 1, 1, -1, 0, -1, 1, 0, 0, 1, -1, -1, 0, 1, -1, -1, 0, -1, 1, 0, 0, 1, -1, -1, -1, 1, 1, -1, 0, 1, 1, 1, 0, -1, 0, 1, 0, 1, 1, 1, 1, -1], 'A355837': [1, -2, -2, 0, -2, 4, -2, 2, 0, 4, -2, 0, -2, 4, 4, 0, -2, 0, -2, 0, 4, 4, -2, -4, 0, 4, 2, 0, -2, -8, -2, -2, 4, 4, 4, 0, -2, 4, 4, -4, -2, -8, -2, 0, 0, 4, -2, 0, 0, 0, 4, 0, -2, -4, 4, -4, 4, 4, -2, 0, -2, 4, 0, 0, 4, -8, -2, 0, 4, -8, -2, 0, -2, 4, 0, 0, 4, -8, -2, 0, 0, 4, -2, 0, 4, 4, 4, -4, -2, 0, 4, 0, 4, 4, 4, 4], 'A369255': [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0]}
    arrays = {'A038548': g, 'A140773': h, 'A007875': a,
              'A219009': b, 'A355837': c, 'A369255': [x % 2 for x in h]}
    for seq, expected in source_prefixes.items():
        assert arrays[seq][1:len(expected) + 1] == expected, ('source prefix', seq)
        print('FIDELITY', seq, 'published terms matched:', len(expected))

    ab = convolution(a, b)
    hc = convolution(h, c)
    for n in range(1, LIMIT + 1):
        es = exponents(n, spf)
        triple_count = prod(comb(e + 2, 2) for e in es)
        square_divisors = prod(e // 2 + 1 for e in es)
        assert 2 * h[n] == triple_count + square_divisors, ('count bridge', n)
        assert h[n] % 2 == ab[n] == hc[n] == classification(es), ('target', n)
    print('SOURCE_REPORTED_RANGE_REPRODUCED: both original convolutions match H(n)%2 for 1 <= n <=', LIMIT)
    print('EXPECTED_SETTLEMENT_CHECK: exponent classification and count bridge match on the same range')

    count = 0
    for length in range(7):
        for es in combinations_with_replacement(range(1, 17), length):
            two_h = prod(comb(e + 2, 2) for e in es) + prod(e // 2 + 1 for e in es)
            assert two_h % 2 == 0, ('integrality', es)
            s = prod(e // 2 + 1 for e in es)
            t = prod(e + 1 if e % 2 == 0 else e + 2 for e in es)
            assert two_h == s * (t + 1), ('factorization', es)
            assert (two_h // 2) % 2 == classification(es), ('signature parity', es)
            count += 1
    print('EXPONENT_MULTISET_CHECK: passed', count, 'multisets; 0..6 exponents, each 1..16')
    for n in [1, 2, 3, 4, 6, 16, 36, 81, 96, 256, 900]:
        print('EXAMPLE', n, 'H=', h[n], 'A*B=', ab[n], 'H*C=', hc[n])
    assert h[2] % 2 == h[3] % 2 == 0 and h[6] % 2 == 1
    print('NONMULTIPLICATIVITY_CHECK: a(2)=a(3)=0, a(6)=1')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

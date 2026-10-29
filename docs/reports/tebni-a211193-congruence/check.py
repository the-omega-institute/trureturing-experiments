#!/usr/bin/env python3
"""Exact finite checks for OEIS A211193; no all-index proof is certified here."""
from fractions import Fraction as Q
from math import comb, factorial
from pathlib import Path
from hashlib import sha256

N = 200
SOURCE_TERMS = [1, 1, 3, 10, 45, 221, 1315, 8324, 60809, 464113,
                3993811, 35342814, 349085869, 3486862653, 38870528411,
                429139127416, 5345350992113, 63994963427393,
                887692696733827, 11284513262684914, 175285847038616301,
                2298693217837384957, 40805829165456572691]


def exp_zero(f):
    """Ordinary rational coefficients of exp(f), f[0] = 0, via y' = f'y."""
    assert f[0] == 0
    y = [Q(1)]
    for n in range(1, len(f)):
        y.append(sum((k * f[k] * y[n-k] for k in range(1, n+1)), Q(0)) / n)
    return y


def ordinary_route():
    log = [Q(0)] + [Q((-1)**(n+1), n) for n in range(1, N+1)]
    exponent = [log[0]] + [log[n] + log[n-1] for n in range(1, N+1)]
    b = exp_zero(exponent)
    a = exp_zero([Q(0)] + b[1:])
    return ([factorial(n) * b[n] for n in range(N+1)],
            [factorial(n) * a[n] for n in range(N+1)])


def differential_route():
    # B' = (1 + log(1+x))B and A' = B'A, using normalized integer coefficients.
    t = [1] + [(-1)**(j+1) * factorial(j-1) for j in range(1, N)]
    b, a = [1], [1]
    for n in range(N):
        b.append(sum(comb(n, j) * t[j] * b[n-j] for j in range(n+1)))
        a.append(sum(comb(n, j) * b[j+1] * a[n-j] for j in range(n+1)))
    return b, a


def stirling_route():
    s = [[0]*(N+1) for _ in range(N+1)]
    s[0][0] = 1
    for m in range(N):
        for k in range(m+2):
            s[m+1][k] = (s[m][k-1] if k else 0) - m*s[m][k]
    c = [1] + [sum(factorial(n)//factorial(n-k) * s[n-k][k]
                   for k in range(1, n//2+1)) for n in range(1, N+1)]
    b = [1] + [c[n] + n*c[n-1] for n in range(1, N+1)]
    assert c[1] == 0
    assert all(c[n] % n == 0 for n in range(1, N+1))
    return b


def bell_comtet_route():
    # A008296's four-term recurrence, with the (0,0) entry adjoined.
    t = [[0]*(N+1) for _ in range(N+1)]
    t[0][0] = 1
    bell = [1]
    a = [1]
    for n in range(1, N+1):
        bell.append(sum(comb(n-1, j)*bell[j] for j in range(n)))
        for k in range(1, n+1):
            t[n][k] = (t[n-1][k-1] + (k-n+1)*t[n-1][k]
                       + ((n-1)*t[n-2][k-1] if n >= 2 else 0))
        a.append(sum(bell[k]*t[n][k] for k in range(1, n+1)))
    return a


def transfer_check(b, a, limit=40):
    # E_m(H) = b_(m+1)/(m+1), H = (B-1-x)/x. Negative weights are allowed.
    h = [0] + [b[m+1]//(m+1) for m in range(1, limit+1)]
    assert any(v < 0 for v in h)
    d = [[0]*(limit+1) for _ in range(limit+1)]
    d[0][0] = 1
    for m in range(1, limit+1):
        for k in range(1, m+1):
            d[m][k] = sum(comb(m-1, j-1)*h[j]*d[m-j][k-1]
                          for j in range(1, m+1))
    for n in range(1, limit+1):
        total = 0
        for k in range(1, n+1):
            m = n-k
            e = sum(comb(m, j)*d[j][k] for j in range(m+1))
            falling = factorial(n)//factorial(n-k)
            assert falling % n == 0
            total += falling*e
        assert a[n]-1 == total, n


def main():
    bq, aq = ordinary_route()
    b, a = differential_route()
    assert bq == b and aq == a
    assert all(v.denominator == 1 for v in bq+aq)
    print(f'Exact rational series and independent integer ODE recurrences agree: n=0..{N}.')
    assert a == bell_comtet_route()
    assert b == stirling_route()
    print(f'Signed Stirling B formula and OEIS Bell/Comtet A formula agree: n=0..{N}.')
    assert a[:len(SOURCE_TERMS)] == SOURCE_TERMS
    fixture = Path(__file__).with_name('b211193.txt').read_bytes()
    table = {}
    for line in fixture.decode().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            n, v = map(int, line.split())
            table[n] = v
    assert set(table) == set(range(N+1))
    assert all(a[n] == table[n] for n in range(N+1))
    print(f'Fidelity: all {len(table)} OEIS b-file entries agree; sha256={sha256(fixture).hexdigest()}.')
    assert all(b[n] % n == 0 for n in range(2, N+1))
    assert all((a[n]-1) % n == 0 for n in range(2, N+1))
    print(f'Expected settlement finite check: n divides b_n and a_n-1 for n=2..{N}.')
    transfer_check(b, a)
    print('Divided-power transfer identity agrees: n=1..40, including negative H weights.')
    for n in [2, 3, 4, 6, 12, 31]:
        print(f'a({n})={a[n]}')
    assert all(a[n] < 0 for n in range(31, N+1, 2))
    print('Source negative odd-index observation agrees for odd n=31..199; no infinite sign claim tested.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

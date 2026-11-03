#!/usr/bin/env python3
"""Exact, independent finite checks for Hanna's OEIS A258655 sign conjecture."""

from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def divisors(n):
    require(n > 0, 'positive divisor domain')
    result = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            result.append(d)
            if d * d != n:
                result.append(n // d)
    return sorted(result)


def published_formula(N):
    ds = divisors(N)
    return (-sum(ds) + sum(d for d in ds if d % 4 == 2)
            + 2 * sum(d for d in ds if d % 8 in (1, 4, 7)))


def theta_coefficients(degree):
    f = [0] * (degree + 1)
    f[0] = 1
    for j in range(1, isqrt(degree) + 1):
        f[j * j] += 1
        if 2 * j * j <= degree:
            f[2 * j * j] += 1
    return f


def logarithmic_coefficients(f):
    # If B_N = N [q^N] log F, the coefficient identity q F' = F sum B_N q^N
    # gives this triangular recurrence directly from the NAME's sparse F.
    require(f[0] == 1, 'unit constant coefficient')
    support = [j for j in range(1, len(f)) if f[j]]
    b = [0] * len(f)
    for N in range(1, len(f)):
        v = N * f[N]
        for j in support:
            if j >= N:
                break
            v -= f[j] * b[N - j]
        b[N] = v
    return b


def logarithm_by_definition(f):
    # A separate exact rational implementation of log(1+H) = sum (-1)^(k+1) H^k/k.
    degree = len(f) - 1
    h = f[:]
    h[0] = 0
    power = [1] + [0] * degree
    answer = [Fraction(0) for _ in f]
    for k in range(1, degree + 1):
        power = [sum(power[i] * h[N - i] for i in range(N + 1))
                 for N in range(degree + 1)]
        for N in range(degree + 1):
            answer[N] += Fraction((-1) ** (k + 1) * power[N], k)
    return answer


def chi(d):
    require(d % 2 == 1, 'odd character domain')
    return 1 if d % 8 in (1, 7) else -1


def prime_factors(n):
    out = []
    p = 2
    while p * p <= n:
        r = 0
        while n % p == 0:
            n //= p
            r += 1
        if r:
            out.append((p, r))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def odd_part_sums(m):
    ds = divisors(m * m)
    return sum(ds), sum(chi(d) * d for d in ds)


def square_factorization(n):
    e, m = 0, n
    while m % 2 == 0:
        e, m = e + 1, m // 2
    S, T = odd_part_sums(m)
    product = 1
    for p, r in prime_factors(m):
        product *= sum((chi(p) * p) ** j for j in range(2 * r + 1))
    require(T == product and 0 < T <= S, f'odd factorization at {m}')
    return T if e == 0 else T + (12 - 2 ** (2 * e + 1)) * S


def main():
    fixture = Path(__file__).with_name('b258655.txt').read_bytes()
    expected_hash = 'f9638076cc38cc9467a77ce6bbd5b7e19db330dd4e6e7a6afc67566560fa0c89'
    require(sha256(fixture).hexdigest() == expected_hash, 'source fixture checksum')
    rows = [tuple(map(int, line.split())) for line in fixture.decode().splitlines()
            if line.strip() and not line.startswith('#')]
    require([n for n, _ in rows] == list(range(1, 1001)), 'source table domain')
    for n, expected in rows:
        actual = published_formula(n * n)
        require(actual == expected, f'source table at {n}')
        require(actual == square_factorization(n), f'square factorization at {n}')
        require(actual < 0 if n % 4 == 0 else actual > 0, f'sign at {n}')
    print('FIDELITY: all 1000 published A258655 b-file values match the divisor formula.')
    print('FACTOR_AND_SIGN: all n = 1..1000 match the proposed factorization and strict signs.')

    degree = 10000
    b = logarithmic_coefficients(theta_coefficients(degree))
    for N in range(1, degree + 1):
        require(b[N] == published_formula(N), f'NAME-to-formula at degree {N}')
    small_log = logarithm_by_definition(theta_coefficients(64))
    for N in range(1, 65):
        require(small_log[N] == Fraction(b[N], N), f'formal-log definition at {N}')
    for n in range(1, 101):
        require(b[n * n] == rows[n - 1][1], f'NAME square-index fidelity at {n}')
        require(b[n * n] < 0 if n % 4 == 0 else b[n * n] > 0,
                f'NAME square-index sign at {n}')
    print('NAME: exact qF\' / F recurrence through degree 10000 agrees with the divisor formula.')
    print('LOG_DEFINITION: rational truncated log(1+H) agrees through degree 64.')
    print('NAME_SQUARE_SIGNS: all n = 1..100 match the source table and strict signs.')

    units = (1, 3, 5, 7)
    require(all(chi(x * y) == chi(x) * chi(y) for x in units for y in units),
            'mod-8 character multiplication table')
    for m in range(1, 102, 2):
        for e in range(9):
            n = (2 ** e) * m
            require(published_formula(n * n) == square_factorization(n),
                    f'two-adic edge case at {(e, m)}')
    print('BOUNDARIES: all 16 character products; e = 0..8, odd m = 1..101 pass.')
    for n in (1, 2, 3, 4, 5, 8, 12, 16, 20, 32, 99, 100):
        print(f'a({n}) = {published_formula(n*n)}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

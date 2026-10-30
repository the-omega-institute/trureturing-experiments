#!/usr/bin/env python3
"""Exact, finite checks for Peter Bala's OEIS A281428 conjecture 1.

Only Python's standard library is needed. The optional b-file is read from
--b-file or downloaded from OEIS. Finite checks do not prove eventual divisibility.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, gcd
from pathlib import Path
from urllib.request import urlopen


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def ode_coefficients(degree):
    """EGF coefficients of C'=S*C^4, S'=C^5, C(0)=1, S(0)=0."""
    c, s = [1], [0]
    powers = {k: [] for k in range(2, 7)}
    for j in range(degree + 1):
        for k in range(2, 7):
            previous = c if k == 2 else powers[k - 1]
            powers[k].append(sum(comb(j, i) * c[i] * previous[j - i]
                                 for i in range(j + 1)))
        if j < degree:
            c.append(sum(comb(j, i) * s[i] * powers[4][j - i]
                         for i in range(j + 1)))
            s.append(powers[5][j])
    return c, s, powers


def root_coefficients(nmax):
    """Independent ordinary-series route: C^2=sum binom(3n,n)x^(2n)/3^n."""
    q = [Fraction(1)]
    for n in range(1, nmax + 1):
        square = Fraction(comb(3 * n, n), 3 ** n)
        q.append((square - sum(q[i] * q[n - i] for i in range(1, n))) / 2)
    return q


def supported_denominator(q):
    d = q.denominator
    for p in (2, 3):
        while d % p == 0:
            d //= p
    return d == 1


def exponential_egf(inner):
    e = [1]
    for n in range(len(inner) - 1):
        e.append(sum(comb(n, i) * e[i] * inner[n - i + 1]
                     for i in range(n + 1)))
    return e


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--b-file', type=Path)
    args = parser.parse_args()
    url = 'https://oeis.org/A281428/b281428.txt'
    raw = args.b_file.read_bytes() if args.b_file else urlopen(url, timeout=60).read()
    table = {}
    for line in raw.decode().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            n, value = map(int, line.split())
            table[n] = value
    nmax = 200
    require(sorted(table) == list(range(nmax + 1)), 'expected b-file indices 0..200')
    degree = 2 * nmax
    c, s, powers = ode_coefficients(degree)
    a = c[::2]
    require(all(a[n] == table[n] for n in range(nmax + 1)), 'b-file mismatch')
    print(f'SOURCE {url}')
    print(f'BFILE_SHA256 {sha256(raw).hexdigest()}')
    print(f'FIDELITY_BFILE matched {nmax + 1} values, indices 0..{nmax}')
    print(f'PREFIX {a[:11]}')

    residues = [value % 385 for value in a]
    require(residues[:14] == [1, 1, 17, 95, 315, 105] + [0] * 8,
            'source mod-385 example mismatch')
    require(all(r == 0 for r in residues[6:]), 'source mod-385 onset fails')
    print(f'FIDELITY_MOD385 {residues[:14]}; zero for all checked n=6..{nmax}')

    for n in range(2, nmax + 1):
        rhs = ((288 * n * n - 576 * n + 317) * a[n - 1]
               - 9 * (2 * n - 3) ** 2 * (6 * n - 11) * (6 * n - 7) * a[n - 2])
        require(16 * a[n] == rhs, f'source scalar recurrence fails at {n}')
    print(f'FIDELITY_SCALAR_RECURRENCE n=2..{nmax}')

    require(all(c[j] == 0 for j in range(1, degree + 1, 2)), 'C not even')
    require(all(s[j] == 0 for j in range(0, degree + 1, 2)), 'S not odd')
    for j in range(degree + 1):
        delta = 1 if j == 0 else 0
        require(4 * powers[6][j] - 3 * powers[2][j] - delta
                - (9 * j * (j - 1) * powers[6][j - 2] if j >= 2 else 0) == 0,
                f'source algebraic equation fails at degree {j}')
    print(f'IDENTITIES parity and (4-9*x^2)*C^6-3*C^2-1 through degree {degree}')

    integral = [0] + powers[4][:-1]
    ep = exponential_egf(integral)
    em = exponential_egf([-v for v in integral])
    require(all(ep[j] + em[j] == 2 * c[j] for j in range(degree + 1)),
            'literal cosh(integral C^4) fixed point fails')
    require(all(ep[j] - em[j] == 2 * s[j] for j in range(degree + 1)),
            'literal sinh(integral C^4) identity fails')
    print(f'LITERAL_FIXED_POINT C=cosh(integral_0^x C^4), S=sinh(integral_0^x C^4), degree 0..{degree}')

    q = root_coefficients(nmax)
    for n in range(nmax + 1):
        require(q[n] * factorial(2 * n) == a[n], f'independent root mismatch {n}')
        require(supported_denominator(q[n]), f'denominator has a prime outside 2,3 at {n}')
    print(f'INDEPENDENT_ROOT_AND_DENOMINATORS {nmax + 1} even coefficients; only primes 2,3')

    moduli = comparisons = 0
    for modulus in range(1, degree + 1):
        if gcd(modulus, 6) != 1:
            continue
        onset = (modulus + 1) // 2
        require(onset <= nmax, 'vacuous onset check')
        moduli += 1
        for n in range(onset, nmax + 1):
            require(a[n] % modulus == 0, f'divisibility fails: N={modulus}, n={n}')
            comparisons += 1
    print(f'EXPECTED_SETTLEMENT {moduli} moduli 1<=N<={degree}, gcd(N,6)=1; '
          f'{comparisons} comparisons at ceil(N/2)<=n<={nmax}')
    print('LIMIT finite exact checks only; universal claim and negation are not proved in Lean')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

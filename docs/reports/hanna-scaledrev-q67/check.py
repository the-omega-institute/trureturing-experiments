#!/usr/bin/env python3
"""Exact finite checks of Hanna's q=6,7 scaled-reversion conjectures.

No conjectural formula generates the coefficients.  At degree n the literal
composition A(x-x*A(q*x)/q)=x has diagonal coefficient a[n], with every
other term depending on a[1],...,a[n-1].  Computing in Z/mZ is exact: this
triangular recurrence uses only addition and multiplication, never division.
"""
from fractions import Fraction
import hashlib
from urllib.request import Request, urlopen


def coefficients(q, count, modulus=None):
    """Solve the literal source equation through degree count, without division."""
    def norm(z):
        return z if modulus is None else z % modulus
    a = [0] * (count + 1)
    g = [0] * (count + 1)
    powers = [[0] * (count + 1) for _ in range(count + 1)]
    a[1] = g[1] = powers[1][1] = 1
    for n in range(2, count + 1):
        g[n] = norm(-pow(q, n - 2) * a[n - 1])
        powers[1][n] = g[n]
        for k in range(2, n + 1):
            powers[k][n] = norm(sum(g[j] * powers[k - 1][n - j]
                                     for j in range(1, n - k + 2)))
        a[n] = norm(-sum(a[k] * powers[k][n] for k in range(1, n)))
        assert norm(sum(a[k] * powers[k][n] for k in range(1, n + 1))) == 0
    return a


def mul(f, g, count):
    return [sum(f[j] * g[n - j] for j in range(n + 1))
            for n in range(count + 1)]


def literal_rational_check(q, a, count=12):
    """Separate Horner composition, using the actual rational division by q."""
    g = [Fraction(0)] * (count + 1)
    g[1] = Fraction(1)
    for n in range(2, count + 1):
        g[n] = -Fraction(a[n - 1] * q ** (n - 1), q)
    value = [Fraction(0)] * (count + 1)
    for k in range(count, -1, -1):
        value = mul(value, g, count)
        value[0] += a[k]
    assert value == [Fraction(int(n == 1)) for n in range(count + 1)]


def btable(seq):
    url = f'https://oeis.org/{seq}/b{seq[1:]}.txt'
    with urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=60) as r:
        raw = r.read()
    rows = {}
    for line in raw.decode().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            n, a = line.split()
            rows[int(n)] = int(a)
    assert set(range(1, 81)).issubset(rows)
    print(f'{seq} source={url} sha256={hashlib.sha256(raw).hexdigest()}', flush=True)
    return rows


def quadratic_inverse_check(residues, count):
    # Independently form B(X+X^2) by the binomial theorem modulo 2.
    from math import comb
    for n in range(count + 1):
        c = sum((residues[k] % 2) * comb(k, n - k)
                for k in range((n + 1) // 2, n + 1)) % 2
        assert c == int(n == 1), (n, c)
    # The opposite composition (X+X^2)(B)=B+B^2 is checked separately.
    for n in range(count + 1):
        c = (residues[n] + (residues[n // 2] if n % 2 == 0 else 0)) % 2
        assert c == int(n == 1), (n, c)


def main():
    for q, seq, modulus in [(6, 'A393858', 14), (7, 'A393859', 8)]:
        published = btable(seq)
        integers = coefficients(q, 80)
        assert all(integers[n] == published[n] for n in range(1, 81))
        literal_rational_check(q, integers)
        print(f'q={q}: first 80 exact integers match the published table; '
              'literal rational composition through degree 12 equals X', flush=True)
        print(f'q={q}: a(1..8)={integers[1:9]}', flush=True)
        residues = coefficients(q, 500, modulus)
        assert all(residues[n] == integers[n] % modulus for n in range(81))
        if q == 6:
            assert all((residues[n] % 2 == 1) == (n & (n - 1) == 0)
                       for n in range(1, 501))
            assert all(residues[n] % 7 == 1 for n in range(1, 501))
            quadratic_inverse_check(residues, 500)
            print('q=6: parity iff power of two holds for every 1<=n<=500; '
                  'all 500 coefficients equal 1 mod 7', flush=True)
            print('q=6: odd indices=' + str([n for n in range(1, 501)
                                             if residues[n] % 2]), flush=True)
            print('q=6: both B(X+X^2)=X and B+B^2=X hold through degree 500 over F2', flush=True)
        else:
            assert all(residues[n] == 1 for n in range(1, 501))
            print('q=7: all 500 coefficients equal 1 mod 8', flush=True)
    print('Arithmetic: integers, Fraction and exact quotient rings only; '
          'finite checks do not prove the unbounded conjectures.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

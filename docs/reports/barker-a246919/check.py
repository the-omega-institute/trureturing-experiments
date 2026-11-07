#!/usr/bin/env python3
"""Exact checks for Colin Barker's OEIS A246919 generating-function conjecture.

The source's initial values are credited to Colin Barker and OEIS A246919,
revision 27. Arithmetic uses only Python integers and SymPy exact expressions.
No finite check here proves the unbounded conjecture.
"""
import argparse
from hashlib import sha256
from math import isqrt
from pathlib import Path

import sympy as sp

SOURCE_PREFIX = [
    0, 0, 7, 0, 19, 14, 37, 13, 61, 38, 91, 28, 127, 74, 169, 49,
    217, 122, 271, 76, 331, 182, 397, 109, 469, 254, 547, 148, 631,
    338, 721, 193, 817, 434, 919, 244, 1027, 542, 1141, 301, 1261,
    662, 1387, 364, 1519, 794, 1657, 433, 1801, 938, 1951, 508, 2107,
]
NUMERATOR = [0, 0, 0, 7, 0, 19, 14, 16, 13, 4, -4, 1, -11, 1, 2, 0, 4]


def predicted(n):
    """Proposed maximum, with the three endpoint-only exceptions removed."""
    if n in (0, 1, 2, 4):
        return 0
    if n % 2:
        k = (n - 1) // 2
        return 3 * k * k + 3 * k + 1
    if n % 4 == 2:
        k = (n - 2) // 4
        return 6 * k * k + 6 * k + 2
    k = n // 4
    return 3 * k * k + 1


def factor_solutions(n):
    """Enumerate all eligible (x,y), including both orientations of each pair.

    U=2y-2x+n, V=2y+2x-n; UV=3n^2, U+V=4y,
    V-U+2n=4x. Positive y forces U,V positive.
    """
    product = 3 * n * n
    result = set()
    for d in range(1, isqrt(product) + 1):
        if product % d:
            continue
        e = product // d
        for u, v in {(d, e), (e, d)}:
            if (u + v) % 4 or (v - u + 2 * n) % 4:
                continue
            x, y = (v - u + 2 * n) // 4, (u + v) // 4
            assert y > 0
            assert x * x - y * y - n * x + n * n == 0
            if x not in (0, n):
                result.add((x, y))
    return result


def source_pari_maximum(n):
    """Literal integer translation of the source's PARI longest(n) program."""
    m = 12 * n * n
    values = []
    for f in range(1, isqrt(m) + 1):
        if m % f:
            continue
        g = m // f
        if (f + g) % 2 == 0:
            x = (f + g) // 2
            if x % 4 == 0:
                values.append(x // 4)
    if len(values) == 1:
        return 0
    return next(y for y in values if y != n)


def direct_solutions(n):
    """Independent scan of x in [-3n^2,3n^2], checking the original equation.

    This interval contains all solutions: the factorization gives
    |2x-n|=|V-U|/2 <= (3n^2-1)/2. Every square test uses isqrt.
    """
    result = set()
    for x in range(-3 * n * n, 3 * n * n + 1):
        if x in (0, n):
            continue
        square = x * x - n * x + n * n
        y = isqrt(square)
        if y > 0 and y * y == square:
            result.add((x, y))
    return result


def rational_series_coefficients(limit):
    """Exact coefficient division by D=(1-X^4)^3, constant coefficient one."""
    coefficients = []
    for n in range(limit + 1):
        value = NUMERATOR[n] if n < len(NUMERATOR) else 0
        for shift, multiplier in [(4, 3), (8, -3), (12, 1)]:
            if n >= shift:
                value += multiplier * coefficients[n - shift]
        coefficients.append(value)
    return coefficients


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bfile', type=Path, help='Optional downloaded OEIS b246919.txt')
    args = parser.parse_args()
    print('Arithmetic: Python integers and SymPy exact rational algebra')
    print('Source: https://oeis.org/A246919 ; revision 27, 2026-06-28')
    for n, source_value in enumerate(SOURCE_PREFIX, start=1):
        assert source_pari_maximum(n) == source_value, n
        assert predicted(n) == source_value, n
    print(f'FIDELITY: source PARI program reproduces all {len(SOURCE_PREFIX)} printed terms')

    for n in range(1, 301):
        solutions = factor_solutions(n)
        actual = max((y for _, y in solutions), default=0)
        assert actual == predicted(n) == source_pari_maximum(n), (n, actual)
    print('MAXIMUM: complete admissible factor-pair enumeration agrees for n=1..300')
    for n in range(1, 81):
        assert direct_solutions(n) == factor_solutions(n), n
    print('CROSS_CHECK: direct original-equation scan agrees as sets for n=1..80')

    limit = 10000
    coefficients = rational_series_coefficients(limit)
    assert coefficients == [predicted(n) for n in range(limit + 1)]
    print('GENERATING_FUNCTION: exact coefficients agree for n=0..10000')
    if args.bfile:
        data = args.bfile.read_bytes()
        rows = [tuple(map(int, line.split())) for line in data.decode().splitlines()
                if line.strip() and not line.lstrip().startswith('#')]
        assert [n for n, _ in rows] == list(range(1, 10001))
        assert all(predicted(n) == value == coefficients[n] for n, value in rows)
        print('FIDELITY: all 10000 source b-file terms agree')
        print('BFILE_SHA256:', sha256(data).hexdigest())

    x, t = sp.symbols('x t')
    s0, s1, s2 = 1 / (1 - t), t / (1 - t) ** 2, t * (1 + t) / (1 - t) ** 3
    odd = x * (3 * s2 + 3 * s1 + s0).subs(t, x ** 2) - x
    two_mod_four = x ** 2 * (6 * s2 + 6 * s1 + 2 * s0).subs(t, x ** 4) - 2 * x ** 2
    zero_mod_four = (3 * s2 + s0).subs(t, x ** 4) - 1 - 4 * x ** 4
    numerator = x ** 3 * (7 + 19 * x ** 2 + 14 * x ** 3 + 16 * x ** 4 +
                          13 * x ** 5 + 4 * x ** 6 - 4 * x ** 7 + x ** 8 -
                          11 * x ** 9 + x ** 10 + 2 * x ** 11 + 4 * x ** 13)
    denominator = (1 - x) ** 3 * (1 + x) ** 3 * (1 + x ** 2) ** 3
    assert sp.expand(denominator - (1 - x ** 4) ** 3) == 0
    assert sp.cancel(odd + two_mod_four + zero_mod_four - numerator / denominator) == 0
    print('SYMBOLIC_IDENTITY: residue-class sum minus printed rational function = 0')
    assert sp.Poly(sp.expand(numerator), x).nth(16) == 4

    actual16 = max(y for _, y in factor_solutions(16))
    recurrence16 = 3 * predicted(12) - 3 * predicted(8) + predicted(4)
    assert actual16 == 49 and recurrence16 == 45
    assert (55, 49) in factor_solutions(16)
    print('SEPARATE_SOURCE_CORRECTION: printed recurrence at n=16 gives 45; actual maximum=49')
    print('SCOPE: finite checks and exact symbolic identity; no Lean proof of the conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

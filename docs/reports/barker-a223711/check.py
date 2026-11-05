#!/usr/bin/env python3
"""Exact fidelity and finite checks for Hardin--Barker OEIS A223711.

Only Python's standard library is required. These checks are not a proof of the
unbounded counting claim. Source: https://oeis.org/A223711, revision 9.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
from pathlib import Path
import hashlib

DISPLAYED = (
    4, 16, 53, 155, 418, 1066, 2611, 6205, 14408, 32852, 73825,
    163951, 360574, 786574, 1704095, 3670193, 7864516, 16777432,
    35651821, 75497731, 159383834, 335544626, 704643403, 1476395365,
    3087008128, 6442451356, 13421773241, 27917287895,
)


def unimodal(s):
    """Literal unmarked weak peak; include the empty auxiliary word."""
    return not s or any(
        all(s[i] <= s[i + 1] for i in range(k))
        and all(s[i] >= s[i + 1] for i in range(k, len(s) - 1))
        for k in range(len(s))
    )


def source_formula(n):
    return Fraction(-2 * (2 ** (3 + n) - 9)
                    + (5 + 2 ** (3 + n)) * n + n * n, 2)


def expanded_formula(n):
    return Fraction(n * n + 5 * n + 18 + (8 * n - 16) * 2 ** n, 2)


def normal_forms(n):
    """Enumerate length tuples, without choosing or marking a summit."""
    yield (0,) * n
    for a in range(n + 1):
        for b in range(1, n - a + 1):
            c = n - a - b
            yield (0,) * a + (1,) * b + (0,) * c
    for a in range(n + 1):
        for b in range(n - a + 1):
            for c in range(1, n - a - b + 1):
                for d in range(n - a - b - c + 1):
                    e = n - a - b - c - d
                    yield (0,) * a + (1,) * b + (2,) * c + (1,) * d + (0,) * e


def trim(p):
    p = list(map(Fraction, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(*polys):
    out = [Fraction(0)] * max(map(len, polys))
    for p in polys:
        for i, v in enumerate(p):
            out[i] += v
    return trim(out)


def mul(*polys):
    out = [Fraction(1)]
    for p in polys:
        new = [Fraction(0)] * (len(out) + len(p) - 1)
        for i, u in enumerate(out):
            for j, v in enumerate(p):
                new[i + j] += u * v
        out = trim(new)
    return out


def scaled(c, p):
    return trim([c * v for v in p])


def main():
    print('Source: OEIS A223711 revision 9; Hardin (2013), Barker conjectures (2018).')
    print('Arithmetic: Python integers and fractions.Fraction; no floating point.')
    # Independent literal arrays: two labeled columns, no symmetry quotient.
    for n in range(9):
        fibers = Counter()
        accepted = 0
        for entries in product(range(2), repeat=2 * n):
            rows = tuple(entries[2 * i] + entries[2 * i + 1] for i in range(n))
            cols = tuple(sum(entries[2 * i + j] for i in range(n)) for j in range(2))
            fibers[rows] += 1
            assert unimodal(cols)
            accepted += unimodal(rows) and unimodal(cols)
        assert all(size == 2 ** word.count(1) for word, size in fibers.items())
        forms = list(normal_forms(n))
        assert len(forms) == len(set(forms)), ('duplicate normal form', n)
        literal_words = {word for word in product(range(3), repeat=n) if unimodal(word)}
        assert set(forms) == literal_words, ('normal form coverage', n)
        weighted = sum(2 ** word.count(1) for word in forms)
        assert accepted == weighted == source_formula(n) == expanded_formula(n)
        if n:
            assert accepted == DISPLAYED[n - 1]
        print(f'Literal arrays n={n}: total={2 ** (2 * n)}, admissible={accepted}, '
              f'unimodal_sum_words={len(forms)}; fibers and unique forms OK')
    # Cross-multiply the proposed generating functions in Q[x], exactly.
    u, v, x = [1, -1], [1, -2], [0, 1]
    den = mul(u, u, u, v, v)
    numerator = mul(x, [2, -3, 2], [2, -3, 2])
    aux_numerator = add(mul(u, u, v, v), mul([0, 2], u, v), x)
    assert den == list(map(Fraction, [1, -7, 19, -25, 16, -4]))
    assert numerator == list(map(Fraction, [0, 4, -12, 17, -12, 4]))
    assert add(aux_numerator, scaled(-1, den)) == numerator
    partial_fraction_numerator = add(
        mul(v, v), mul(u, v, v), scaled(7, mul(u, u, v, v)),
        scaled(4, mul(u, u, u)), scaled(-12, mul(u, u, u, v)),
    )
    assert aux_numerator == partial_fraction_numerator
    print('Exact polynomial cross-products: peak-class GF, source GF, partial fractions OK')
    # Recover coefficients solely from the printed rational numerator/denominator.
    coefficients = []
    for n in range(211):
        rhs = numerator[n] if n < len(numerator) else 0
        coeff = rhs - sum(den[k] * coefficients[n - k]
                          for k in range(1, min(n, len(den) - 1) + 1))
        assert coeff.denominator == 1
        coefficients.append(coeff.numerator)
        aux_coeff = comb(n + 2, 2) + (n + 1) + 7 + 4 * (n + 1) * 2 ** n - 12 * 2 ** n
        assert aux_coeff == source_formula(n) == expanded_formula(n)
        assert coeff == (0 if n == 0 else source_formula(n))
    assert tuple(coefficients[1:29]) == DISPLAYED
    data = Path(__file__).with_name('b223711.txt').read_bytes()
    pairs = [tuple(map(int, line.split())) for line in data.decode().splitlines()
             if line.strip() and not line.lstrip().startswith('#')]
    assert [n for n, _ in pairs] == list(range(1, 211))
    for n, a in pairs:
        assert a == source_formula(n) == coefficients[n]
        if n >= 6:
            assert a == (7 * coefficients[n - 1] - 19 * coefficients[n - 2]
                         + 25 * coefficients[n - 3] - 16 * coefficients[n - 4]
                         + 4 * coefficients[n - 5])
    print('Fidelity: all 28 displayed terms and all 210 Hardin b-file entries match')
    print('Expected-settlement finite check: source formula = GF coefficients for n=1..210')
    print('Auxiliary n=0 count is 1; source GF coefficient at n=0 is 0')
    print('Empirical source order-five recurrence: n=6..210 OK')
    print('b223711.txt sha256=' + hashlib.sha256(data).hexdigest())
    print('Boundary: no Lean proof or all-n enumeration proof is asserted by this script')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact finite checks for Paul D. Hanna's OEIS A381361.

The source data retain their CC-BY-SA-4.0 license; see README.md.
No floating point or alternate elliptic/algebraic source formula is used.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path

LIMIT = 500
DEGREE = 24
SOURCE_PREFIX = [
    1, 1, 3, 7, 41, 225, 1435, 11815, 108945, 1062145,
    12367475, 154967175, 2052583225, 30729193825,
    489419016075, 8162349262375, 149689173742625,
    2891193923066625, 58124967786683875, 1262035025301370375,
    28655826402568055625, 674057118324247590625,
    16913418325411880586875,
]


def integer_recurrence(limit):
    """a = EGF(A), c = EGF(A^2), b = EGF(A^(-2))."""
    a = [1] + [0] * (limit + 1)
    c, b = [], []
    for n in range(limit + 1):
        c.append(sum(comb(n, k) * a[k] * a[n - k]
                     for k in range(n + 1)))
        b.append(1 if n == 0 else
                 -sum(comb(n, k) * c[k] * b[n - k]
                      for k in range(1, n + 1)))
        a[n + 1] = sum(comb(n, k) * a[k] * abs(b[n - k])
                       for k in range(n + 1))
    return a[:limit + 1], c, b


def mul(f, g):
    return [sum((f[k] * g[n - k] for k in range(n + 1)), Q(0))
            for n in range(DEGREE + 1)]


def inverse(f):
    assert f[0] == 1
    g = [Q(1)]
    for n in range(1, DEGREE + 1):
        g.append(-sum((f[k] * g[n - k] for k in range(1, n + 1)), Q(0)))
    return g


def primitive(f):
    return [Q(0)] + [f[n - 1] / n for n in range(1, DEGREE + 1)]


def exp_zero_constant(f):
    assert f[0] == 0
    g = [Q(1)]
    for n in range(1, DEGREE + 1):
        g.append(sum((k * f[k] * g[n - k] for k in range(1, n + 1)), Q(0)) / n)
    return g


def source_step(f):
    return exp_zero_constant(primitive([abs(x) for x in inverse(mul(f, f))]))


def main():
    path = Path(__file__).with_name('source-b381361.txt')
    raw = path.read_bytes()
    rows = [tuple(map(int, line.split())) for line in raw.decode().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert [n for n, _ in rows] == list(range(LIMIT + 1))
    published = [v for _, v in rows]
    a, c, b = integer_recurrence(LIMIT)
    assert a[:len(SOURCE_PREFIX)] == SOURCE_PREFIX
    assert a == published
    print('Source: Paul D. Hanna, OEIS A381361, revision 25 (2026-04-23).')
    print('b-file sha256:', sha256(raw).hexdigest())
    print('Fidelity: all 501 b-file terms at indices 0..500 match the integer recurrence.')
    print('a[0..10] =', a[:11])
    assert c[:5] == [1, 2, 8, 32, 192]
    assert b[:5] == [1, -2, 0, 16, -64]
    print('c[0..4] =', c[:5])
    print('b[0..4] =', b[:5])

    # Independent ordinary rational-coefficient Picard iteration of the literal source map.
    f = [Q(1)] + [Q(0)] * DEGREE
    for _ in range(DEGREE + 1):
        f = source_step(f)
    assert source_step(f) == f
    assert [factorial(n) * f[n] for n in range(DEGREE + 1)] == a[:DEGREE + 1]
    recip = inverse(mul(f, f))
    assert [factorial(n) * recip[n] for n in range(DEGREE + 1)] == b[:DEGREE + 1]
    print('Literal formal equation: rational Picard fixed point agrees through degree 24.')

    carries = [(i, j) for i in range(5) for j in range(5) if i + j >= 5]
    assert all(comb(i + j, i) % 5 == 0 for i, j in carries)
    assert all(x % 5 == 0 for seq in (a, b, c) for x in seq[5:])
    assert all(abs(x) % 5 == 0 for x in b[5:])
    assert a[4] % 5 == 1 and a[5] == 225
    print('Carry helper:', len(carries), 'pairs i,j < 5 with i+j >= 5 checked exactly.')
    print('Expected settlement: a, b and c divisible by 5 at all 496 indices 5..500.')
    print('Sharpness: a[4] = 41 is not divisible by 5; a[5] = 225.')
    print('These are finite checks; the infinite statement still requires its induction proof.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

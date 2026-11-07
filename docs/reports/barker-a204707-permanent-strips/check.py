#!/usr/bin/env python3
"""Exact checks for Hardin/Barker OEIS A204707; finite evidence, not an all-n proof."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import hashlib

HERE = Path(__file__).resolve().parent
ROWS = tuple(product((0, 1), repeat=3))

def perms(r, s):
    return tuple(r[j] * s[j + 1] + r[j + 1] * s[j] for j in range(2))

def literal_source(rows):
    values = [p for r, s in zip(rows, rows[1:]) for p in perms(r, s)]
    return bool(values) and values[0] > 0 and all(p == values[0] for p in values)

def nonsingular(rows):
    return all(r[j] * s[j + 1] - r[j + 1] * s[j] != 0
               for r, s in zip(rows, rows[1:]) for j in range(2))

def bfile(name):
    raw = (HERE / name).read_bytes()
    values = dict(tuple(map(int, line.split())) for line in raw.decode().splitlines()
                  if line.strip() and not line.lstrip().startswith('#'))
    print(f'SOURCE {name} terms={len(values)} sha256={hashlib.sha256(raw).hexdigest()}')
    return values

# Each field pair (a,b) represents exactly a+b*sqrt(17).
def add(x, y):
    return (x[0] + y[0], x[1] + y[1])

def scale(c, x):
    return (c * x[0], c * x[1])

def mul(x, y):
    return (x[0] * y[0] + 17 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def power(x, n):
    y = (F(1), F(0))
    while n:
        if n % 2:
            y = mul(y, x)
        x = mul(x, x)
        n //= 2
    return y

def radical(n):
    z = add(mul(power((F(1), F(-1)), n), (F(-19), F(5))),
            mul(power((F(1), F(1)), n), (F(19), F(5))))
    # Division by sqrt(17): (a+b*sqrt(17))/sqrt(17)=b+(a/17)*sqrt(17).
    return add((F(1), F(0)), scale(F(1, 2 ** (n + 1)), (z[1], z[0] / 17)))

def polyadd(a, b):
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0))
            for i in range(max(len(a), len(b)))]

def polymul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c

assert all(isinstance(x, int) for r in ROWS for x in r)
expected = {
    (0, 0, 0): (), (0, 0, 1): (), (1, 0, 0): (),
    (0, 1, 0): ((1, 0, 1), (1, 1, 1)),
    (0, 1, 1): ((1, 0, 1), (1, 1, 0)),
    (1, 1, 0): ((0, 1, 1), (1, 0, 1)),
    (1, 1, 1): ((0, 1, 0), (1, 0, 1)),
    (1, 0, 1): ((0, 1, 0), (0, 1, 1), (1, 1, 0), (1, 1, 1)),
}
active_a = {(0, 1, 0), (0, 1, 1), (1, 1, 0), (1, 1, 1)}
active_b = {(1, 0, 1)}
neighbors = {r: tuple(s for s in ROWS if perms(r, s) == (1, 1)) for r in ROWS}
assert neighbors == expected
for r in ROWS:
    for s in ROWS:
        for j in range(2):
            det = r[j] * s[j + 1] - r[j + 1] * s[j]
            assert (det != 0) == (perms(r, s)[j] == 1)
for r in active_a:
    assert sum(s in active_a for s in neighbors[r]) == 1
    assert sum(s in active_b for s in neighbors[r]) == 1
assert set(neighbors[(1, 0, 1)]) == active_a
print('ROW_PAIRS 64 checked; all neighbor lists and quotient degrees match')

source_a = bfile('b204707.txt')
source_b = bfile('b183682.txt')
assert set(source_a) == set(range(1, 211))
assert set(source_b) == set(range(1, 201))

# Directly enumerate all labelled arrays, independently of the row-walk recurrence.
for n in range(1, 6):
    ca = cb = 0
    for rows in product(ROWS, repeat=n + 1):
        valid = literal_source(rows)
        ns = nonsingular(rows)
        ones = all(r == (1, 1, 1) for r in rows)
        assert valid == (ones or ns)
        assert not (ones and ns)
        ca += valid
        cb += ns
    assert (ca, cb) == (source_a[n], source_b[n])
    print(f'BRUTE n={n} labelled_arrays={8 ** (n + 1)} a={ca} b={cb} decomposition=OK')

# Nonconsecutive columns change the first term, so they are excluded from the reading.
arbitrary = sum(len({r[j] * s[k] + r[k] * s[j] for j, k in combinations(range(3), 2)}) == 1
                and r[0] * s[1] + r[1] * s[0] > 0 for r in ROWS for s in ROWS)
assert arbitrary == 7 and source_a[1] == 13
print('FIDELITY offset=1; consecutive a(1)=13; arbitrary-column a(1)=7')

# The seven five-row examples printed in the source entry.
example_lines = ['101010111011110011011', '111111010101011110101', '101101111110101011111', '111010101011110110101', '101111011110011101011']
examples = [tuple(tuple(int(x) for x in line[3*k:3*k+3]) for line in example_lines)
            for k in range(7)]
assert all(literal_source(rows) for rows in examples)
print('FIDELITY all 7 source n=4 example arrays satisfy the literal predicate')

# Two independent eight-state transfer propagations, one for each possible positive permanent.
counts_by_p = {p: {r: 1 for r in ROWS} for p in (1, 2)}
a = {0: 0}
b = {0: 5}  # Active rows only: auxiliary walk count, not the source n=0 count.
u = v = 1
for n in range(1, 211):
    for p in (1, 2):
        prev = counts_by_p[p]
        counts_by_p[p] = {s: sum(prev[r] for r in ROWS if perms(r, s) == (p, p)) for s in ROWS}
    b[n] = sum(counts_by_p[1].values())
    p2 = sum(counts_by_p[2].values())
    assert p2 == 1
    a[n] = b[n] + p2
    u, v = u + v, 4 * u
    assert b[n] == 4 * u + v
    assert a[n] == source_a[n]
    if n in source_b:
        assert b[n] == source_b[n]
    assert radical(n) == (F(a[n]), F(0))
    if n >= 2:
        assert b[n] == b[n - 1] + 4 * b[n - 2]
    if n >= 4:
        assert a[n] == 2 * a[n - 1] + 3 * a[n - 2] - 4 * a[n - 3]
print('TRANSFER all 210 A204707 terms and all 200 A183682 terms match exactly')
print('RADICAL exact Q(sqrt(17)) values match for every 1<=n<=210')

# GF coefficient checks and the exact polynomial identity for adding x/(1-x).
d = [1, -2, -3, 4]
num = [0, 13, 7, -24]
for n in range(211):
    assert sum(d[j] * a[n - j] for j in range(min(n, 3) + 1)) == (num[n] if n < 4 else 0)
assert polyadd(polymul([0, 1], [1, -1, -4]), polymul([0, 12, 20], [1, -1])) == num
qmat = ((1, 1), (4, 0))
q2 = tuple(tuple(sum(qmat[i][k] * qmat[k][j] for k in range(2)) for j in range(2)) for i in range(2))
assert q2 == tuple(tuple(qmat[i][j] + (4 if i == j else 0) for j in range(2)) for i in range(2))
print('GF exact polynomial identity and coefficients through x^210 verified; Q^2=Q+4I')
print('FINITE_CHECKS_ONLY: no all-height theorem or Lean proof is claimed')
print('ALL_CHECKS_OK')

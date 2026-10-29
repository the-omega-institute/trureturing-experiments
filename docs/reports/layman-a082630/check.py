#!/usr/bin/env python3
"""Exact finite audits for Layman's A082630; not an infinite proof.

Source: https://oeis.org/A082630, revision 68, NAME, DATA, EXAMPLE,
and Zhuorui He's Python generator. Source data: https://oeis.org/LICENSE.
The successor-leaf/floor construction also appears in
https://github.com/NiuVector/A085260-proof/blob/bed0c63a10c1250e7cbc86e7c7bfdb24be361d80/proof_note_en.md.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

SOURCE_ROWS = [
    [1, 1],
    [1, 2, 1],
    [1, 2, 3, 1],
    [1, 2, 5, 3, 4, 1],
    [1, 2, 5, 8, 3, 7, 4, 5, 1],
]
SOURCE_PREFIX = [
    1, 2, 5, 8, 19, 30, 71, 112, 265, 418, 989, 1560, 3691, 5822,
    13775, 21728, 51409, 81090, 191861, 302632, 716035, 1129438,
    2672279, 4215120, 9973081, 15731042, 37220045, 58709048,
    138907099, 219105150, 518408351, 817711552, 1934726305,
    3051741058, 7220496869,
]


def simultaneous(old):
    """Inspect only old gaps; inserted values cannot create same-round work."""
    new = []
    for x, y in zip(old, old[1:]):
        new.append(x)
        if 7 * y <= 11 * x:
            new.append(x + y)
    return new + old[-1:]


def certified_prefix(row):
    """Every internal gap returned is terminal, so this prefix persists."""
    for i, (x, y) in enumerate(zip(row, row[1:])):
        if 7 * y <= 11 * x:
            return row[:i + 1]
    return row[:]


def source_generator(count):
    """Literal arithmetic of Zhuorui He's generator in OEIS revision 68."""
    out = []
    x = y = 1
    while len(out) < count:
        out.append(x)
        while 7 * (x + y) <= 11 * x:
            y += x
        x += y
    return out


def terminal_tree(count):
    """Left descend to a terminal gap; next use its immediate right sibling."""
    out = []
    x = y = 1
    while len(out) < count:
        while 7 * y <= 11 * x:
            y += x
        assert 0 < x < y
        out.append(x)
        x, y = y, y - x
    return out


def floor_successor(count):
    out = [1, 2]
    while len(out) < count:
        a, b = out[-2:]
        q = 1 + (11 * b + 7 * a) // (7 * b)
        out.append(q * b - a)
    return out[:count]


def recurrence(count):
    out = [1, 2, 5, 8]
    while len(out) < count:
        out.append(4 * out[-2] - out[-4])
    return out[:count]


def matmul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bfile', type=Path, help='optional downloaded OEIS 1000-term table')
    args = parser.parse_args()
    n = 1000
    a = source_generator(n)
    assert a == terminal_tree(n) == floor_successor(n) == recurrence(n)
    assert a[:len(SOURCE_PREFIX)] == SOURCE_PREFIX
    print(f'FIDELITY OEIS DATA: all {len(SOURCE_PREFIX)} displayed entries agree')

    row = [1, 1]
    previous = []
    for t in range(20):
        if t < len(SOURCE_ROWS):
            assert row == SOURCE_ROWS[t], (t, row)
        if t == 5:
            assert row[:6] == [1, 2, 5, 8, 11, 3]
        if t == 6:
            assert row[:6] == [1, 2, 5, 8, 19, 11]
        prefix = certified_prefix(row)
        assert prefix == a[:len(prefix)]
        assert prefix[:len(previous)] == previous
        assert all(11 * x < 7 * y for x, y in zip(prefix, prefix[1:]))
        print(f'SIMULTANEOUS round={t} row_length={len(row)} stable_prefix_length={len(prefix)}')
        previous = prefix
        row = simultaneous(row)
    print('FIDELITY OEIS EXAMPLE: full rounds 0..4 and stated prefixes at rounds 5,6 agree')

    for i in range(n - 2):
        r = F(a[i], a[i + 1])
        lo, hi, q = (F(1, 2), F(7, 11), 3) if i % 2 == 0 else (F(2, 5), F(3, 7), 2)
        assert lo <= r < hi
        assert q == 1 + (F(11, 7) + r).__floor__()
        assert a[i + 2] + a[i] == q * a[i + 1]
    assert all(a[i + 4] + a[i] == 4 * a[i + 2] for i in range(n - 4))
    print(f'EXPECTED SETTLEMENT: {n-2} alternating 3,2 coefficients and {n-4} lag-2 recurrence checks agree')

    # Exact endpoint checks for the proposed interval induction.
    assert 2 <= F(11, 7) + F(1, 2)
    assert F(11, 7) + F(7, 11) < 3
    assert 1 <= F(11, 7) + F(2, 5)
    assert F(11, 7) + F(3, 7) == 2  # excluded endpoint
    assert 1 / (3 - F(1, 2)) == F(2, 5)
    assert 1 / (3 - F(7, 11)) == F(11, 26) < F(3, 7)
    assert F(1, 2) <= 1 / (2 - F(2, 5)) == F(5, 8)
    assert 1 / (2 - F(3, 7)) == F(7, 11)  # excluded endpoint
    b = [[-1, 3], [-2, 5]]
    b2 = matmul(b, b)
    assert b2 == [[-5, 12], [-8, 19]]
    assert [[b2[i][j] + int(i == j) for j in range(2)] for i in range(2)] == [
        [4 * b[i][j] for j in range(2)] for i in range(2)]
    print('EXACT HELPERS: rational interval endpoints and B^2 + I = 4*B agree')

    if args.bfile:
        data = args.bfile.read_bytes()
        pairs = [tuple(map(int, line.split())) for line in data.decode().splitlines()
                 if line.strip() and not line.lstrip().startswith('#')]
        assert [i for i, _ in pairs] == list(range(1, 1001))
        assert [v for _, v in pairs] == a
        print(f'OEIS BFILE: all 1000 entries agree; sha256={sha256(data).hexdigest()}')
    print('FINITE AUDIT ONLY: infinite stabilization and all-index induction remain to be formalized')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

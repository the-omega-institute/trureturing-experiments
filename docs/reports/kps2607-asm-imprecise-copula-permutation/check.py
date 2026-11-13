#!/usr/bin/env python3
"""Exact pre-probe checks for the KPS permutation-copula question.

Only Python's standard library is used. No Lean settlement is attempted.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import factorial


def corners(a):
    n = len(a)
    q = [[0] * (n + 1) for _ in range(n + 1)]
    for r in range(1, n + 1):
        for s in range(1, n + 1):
            q[r][s] = a[r - 1][s - 1] + q[r - 1][s] + q[r][s - 1] - q[r - 1][s - 1]
    return q


def second_differences(q):
    n = len(q) - 1
    return [[q[r][s] + q[r - 1][s - 1] - q[r][s - 1] - q[r - 1][s]
             for s in range(1, n + 1)] for r in range(1, n + 1)]


def check_asm(a):
    for row in a + [list(col) for col in zip(*a)]:
        assert all(x in (-1, 0, 1) for x in row)
        assert sum(row) == 1
        nonzero = [x for x in row if x]
        assert nonzero[0] == nonzero[-1] == 1
        assert all(x == -y for x, y in zip(nonzero, nonzero[1:]))
        partial = 0
        for x in row:
            partial += x
            assert partial in (0, 1)


def fidelity():
    # arXiv:2607.06178v2, Examples 2.1 and 2.2.
    a = [[0, 1, 0, 0], [1, -1, 1, 0], [0, 1, -1, 1], [0, 0, 1, 0]]
    expected = [[0, 1, 1, 1], [1, 1, 2, 2], [1, 2, 2, 3], [1, 2, 3, 4]]
    q = corners(a)
    assert [row[1:] for row in q[1:]] == expected
    assert second_differences(q) == a
    check_asm(a)
    scaled = [[Fraction(x, 4) for x in row] for row in expected]
    assert scaled == [[Fraction(0), Fraction(1, 4), Fraction(1, 4), Fraction(1, 4)],
                      [Fraction(1, 4), Fraction(1, 4), Fraction(1, 2), Fraction(1, 2)],
                      [Fraction(1, 4), Fraction(1, 2), Fraction(1, 2), Fraction(3, 4)],
                      [Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)]]
    # Check the source's boundary-rectangle quasi-copula definition, without
    # replacing it by a test for a proper copula (the example is not one).
    for i in range(5):
        for j in range(i, 5):
            for k in range(5):
                for l in range(k, 5):
                    if i == 0 or j == 4 or k == 0 or l == 4:
                        assert q[i][k] + q[j][l] >= q[i][l] + q[j][k]
    print('FIDELITY Examples 2.1/2.2: scaling, corner sums, inverse, ASM, boundary rectangles OK')


def bracket(n):
    return ([[r * s // n for s in range(n + 1)] for r in range(n + 1)],
            [[(r * s + n - 1) // n for s in range(n + 1)] for r in range(n + 1)])


def family_premises():
    rectangles = 0
    for n in range(1, 41):
        p, q = bracket(n)
        for f in (p, q):
            assert all(f[r][0] == f[0][r] == 0 for r in range(n + 1))
            assert all(f[r][n] == f[n][r] == r for r in range(n + 1))
            assert all(0 <= x <= n for row in f for x in row)
            assert all(f[r + 1][s] - f[r][s] in (0, 1)
                       for r in range(n) for s in range(n + 1))
            assert all(f[r][s + 1] - f[r][s] in (0, 1)
                       for r in range(n + 1) for s in range(n))
            a = second_differences(f)
            check_asm(a)
            assert corners(a) == f
        assert all(n * p[r][s] <= r * s <= n * q[r][s]
                   for r in range(n + 1) for s in range(n + 1))
        # n*C(r,s) = r*s: its boundaries, range and rectangle inequalities
        # are checked as integers. n is strictly positive.
        assert all(0 <= r * s <= n * n for r in range(n + 1) for s in range(n + 1))
        assert all(r * n == n * r for r in range(n + 1))
        for i in range(n + 1):
            for j in range(i, n + 1):
                for k in range(n + 1):
                    for l in range(k, n + 1):
                        assert i * k + j * l - i * l - j * k == (j - i) * (l - k) >= 0
                        assert q[i][k] + p[j][l] - p[i][l] - p[j][k] >= 0
                        assert p[i][k] + q[j][l] - p[i][l] - p[j][k] >= 0
                        assert q[i][k] + q[j][l] - q[i][l] - p[j][k] >= 0
                        assert q[i][k] + q[j][l] - p[i][l] - q[j][k] >= 0
                        rectangles += 1
    print(f'PREMISES n=1..40: {rectangles} rectangles (including degenerate), four IC inequalities, two ASMs, C=rs/n OK')


def permutation_corners(pi):
    n = len(pi)
    a = [[int(pi[i] == j) for j in range(n)] for i in range(n)]
    return corners(a)


def exhaustive_small_brackets():
    for n in range(1, 9):
        p, q = bracket(n)
        seen = witnesses = 0
        for pi in permutations(range(n)):
            r = permutation_corners(pi)
            seen += 1
            if all(p[i][j] <= r[i][j] <= q[i][j]
                   for i in range(n + 1) for j in range(n + 1)):
                witnesses += 1
        assert seen == factorial(n) and witnesses > 0
        print(f'EXHAUSTIVE bracket n={n}: permutations={seen}, witnesses={witnesses}; no small refutation claimed')


def haar(n, start, length):
    return {r: 1 if r < start + length // 2 else -1 for r in range(start, start + length)}


def dot(a, b):
    return sum(v * b.get(i, 0) for i, v in a.items())


def block_basis(n, m):
    blocks = []
    for k in range(2, m):
        ell = m + 1 - k
        length, width = n // (2 ** k), n // (2 ** ell)
        assert length % 2 == width % 2 == 0 and length * width == n // 2
        shape = []
        for a in range(0, n, length):
            hi = haar(n, a, length)
            assert sum(hi.values()) == 0
            assert sum(r * v for r, v in hi.items()) == -(length ** 2) // 4
            assert sum(v for r, v in hi.items() if r >= a) == 0
            assert sum(v for r, v in hi.items() if r >= a + length) == 0
            for b in range(0, n, width):
                hj = haar(n, b, width)
                h = {r * n + s: x * y for r, x in hi.items() for s, y in hj.items()}
                assert dot(h, h) == n // 2
                shape.append((k, a, b, length, width, h))
        assert len(shape) == 2 * n
        blocks.extend(shape)
    for x, y in combinations(blocks, 2):
        assert dot(x[-1], y[-1]) == 0
    return blocks


def check_haar_permutation(pi, m, blocks):
    n = len(pi)
    r = permutation_corners(pi)
    e = {i * n + j: n * r[i][j] - i * j for i in range(n) for j in range(n)}
    hsum = {}
    t = 0
    per_shape = {k: 0 for k in range(2, m)}
    for k, a, b, length, width, h in blocks:
        interior_points = sum(a < i + 1 < a + length and b < pi[i] + 1 < b + width
                              for i in range(n))
        if interior_points:
            continue
        per_shape[k] += 1
        t += 1
        # Verify the cancellation of the counting part independently.
        assert sum(r[x // n][x % n] * v for x, v in h.items()) == 0
        assert dot(e, h) * 64 == -(n ** 2)
        for x, v in h.items():
            hsum[x] = hsum.get(x, 0) + v
    assert all(count >= n for count in per_shape.values())
    assert t >= n * (m - 2)
    assert dot(e, hsum) * 64 == -(n ** 2) * t
    assert dot(hsum, hsum) * 2 == n * t
    energy = dot(e, e)
    squares = sum((32 * v + n * hsum.get(x, 0)) ** 2 for x, v in e.items())
    assert 2048 * energy - n ** 3 * t == 2 * squares >= 0


def multiscale_checks():
    for m in (3, 4, 5, 6):
        n = 2 ** m
        blocks = block_basis(n, m)
        if n == 8:
            sample = permutations(range(n))
            expected_count = factorial(n)
        else:
            # Deterministic examples, with reversals, boundary-rich points and
            # bit-reversal examples in addition to affine permutations.
            values = {tuple(range(n)), tuple(reversed(range(n))),
                      tuple(int(f'{i:0{m}b}'[::-1], 2) for i in range(n))}
            for a in (1, 3, 5, 7):
                for b in (0, 1, n // 2, n - 1):
                    values.add(tuple((a * i + b) % n for i in range(n)))
            sample = sorted(values)
            expected_count = len(sample)
        checked = 0
        for pi in sample:
            assert sorted(pi) == list(range(n))
            check_haar_permutation(pi, m, blocks)
            checked += 1
        assert checked == expected_count
        print(f'HAAR n={n}: all {len(blocks)} basis blocks pairwise orthogonal; {checked} permutations satisfy empty-block coefficients, counts and square identity')


def threshold():
    for m in (2051, 2052, 4096):
        n = 2 ** m
        gap = (m - 2) * n ** 4 - 2048 * n ** 4
        assert gap == (m - 2050) * n ** 4 > 0
        print(f'THRESHOLD m={m}: n bit_length={n.bit_length()}, normalized contradiction gap={m - 2050} > 0')
    print('LIMIT: S_(2^2051) was not enumerated; the universal family requires the stated mathematical argument and subsequent Lean verification')


def main():
    fidelity()
    family_premises()
    exhaustive_small_brackets()
    multiscale_checks()
    threshold()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

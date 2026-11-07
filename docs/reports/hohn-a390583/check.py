#!/usr/bin/env python3
"""Exact checks for Hohn's OEIS A390583 repeated-largest-radius conjecture."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BOUND = 500


def descartes(row):
    curvatures = [Fraction(1, radius) for radius in row]
    return sum(curvatures) ** 2 == 2 * sum(x ** 2 for x in curvatures)


def cube_root(n):
    """Return the exact nonnegative cube root, or None. No real arithmetic."""
    lo, hi = 0, n + 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** 3 <= n:
            lo = mid
        else:
            hi = mid
    return lo if lo ** 3 == n else None


def cube_witness(difference):
    assert difference > 0
    m = cube_root(difference)
    if m is not None:
        return (1, m)
    if difference % 2 == 0:
        m = cube_root(difference // 2)
        if m is not None:
            return (2, m)
    return None


def normalization(row):
    a, b, c, d = row
    assert 0 < a < b <= c == d and gcd(*row) == 1 and descartes(row)
    assert c * (b - a) ** 2 == 4 * a * b * (a + b)
    g = gcd(a, b)
    u, v = a // g, b // g
    t = v - u
    A = 4 * u * v * (u + v)
    D = gcd(t ** 2, A)
    expected_D = 1 if t % 2 else (4 if t % 4 else 8)
    assert gcd(u, v) == gcd(t, u * v) == gcd(g, c) == 1
    assert D == expected_D
    assert t ** 2 % D == A % D == 0
    assert g == t ** 2 // D and c == A // D
    assert D * (b - a) == t ** 3
    witness = cube_witness(b - a)
    assert witness is not None
    factor, m = witness
    assert m >= 1 and b - a == factor * m ** 3
    return D, t


def main():
    data = (ROOT / 'b390583.txt').read_bytes()
    entries = []
    for line in data.decode().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            n, value = map(int, line.split())
            assert n == len(entries) + 1
            entries.append(value)
    assert len(entries) == 4000
    rows = [tuple(entries[i:i + 4]) for i in range(0, len(entries), 4)]
    expected_first_ten = [
        (1, 2, 24, 24), (1, 3, 12, 12), (2, 3, 120, 120),
        (2, 10, 15, 15), (3, 4, 336, 336), (3, 5, 120, 120),
        (3, 8, 22, 264), (4, 5, 720, 720), (4, 11, 33, 132),
        (5, 6, 1320, 1320),
    ]
    assert rows[:10] == expected_first_ten
    assert rows == sorted(set(rows))
    for row in rows:
        assert 0 < min(row) and list(row) == sorted(row)
        assert gcd(*row) == 1 and descartes(row)
    print('SOURCE_DATA_SHA256', sha256(data).hexdigest())
    print('FIDELITY first_ten_rows=10 all_source_rows=1000 exact_equation_and_gcd=OK')
    repeated = [row for row in rows if row[2] == row[3]]
    source_cases = [normalization(row) for row in repeated]
    source_counts = Counter(D for D, t in source_cases)
    assert len(repeated) == 498 and source_counts == {1: 278, 4: 137, 8: 83}
    print('SOURCE_REPEATED', len(repeated), 'D_COUNTS', dict(sorted(source_counts.items())),
          'MAX_T', max(t for D, t in source_cases))

    # All positive a <= b <= BOUND: the cleared equation fixes c uniquely
    # when a < b. The a=b case has positive right side and zero left side.
    checked = []
    for a in range(1, BOUND + 1):
        for b in range(a + 1, BOUND + 1):
            numerator, denominator = 4 * a * b * (a + b), (b - a) ** 2
            c, remainder = divmod(numerator, denominator)
            if remainder == 0 and b <= c and gcd(a, b, c) == 1:
                row = (a, b, c, c)
                checked.append(normalization(row))
    assert checked
    print('SMALL_CASES all_positive_a_le_b_le', BOUND, 'primitive_repeated_rows', len(checked),
          'D_COUNTS', dict(sorted(Counter(D for D, t in checked).items())))

    # Check the proposed gcd classification independently on reduced u,v,
    # including pairs whose repeated radius would not be largest.
    pair_counts = Counter()
    for u in range(1, BOUND + 1):
        for v in range(u + 1, BOUND + 1):
            if gcd(u, v) != 1:
                continue
            t, A = v - u, 4 * u * v * (u + v)
            D = gcd(t ** 2, A)
            assert D == (1 if t % 2 else (4 if t % 4 else 8))
            assert gcd(t ** 2 // D, A // D) == 1
            pair_counts[D] += 1
    print('REDUCED_PAIRS all_coprime_1_le_u_lt_v_le', BOUND,
          'count', sum(pair_counts.values()), 'D_COUNTS', dict(sorted(pair_counts.items())))

    # A source family formula is an auxiliary error, not the target conjecture.
    # A000447(k)=k*(4*k*k-1)/3; at k=3 the literal A390583 family gives c=420.
    k = 3
    auxiliary_radius = 12 * (k * (4 * k * k - 1) // 3)
    assert auxiliary_radius == 420
    assert not descartes((k, k + 2, auxiliary_radius, auxiliary_radius))
    assert (3, 5, 120, 120) in rows
    print('SOURCE_AUXILIARY_DISCREPANCY literal_k3_family=(3,5,420,420) fails_equation;',
          'source_table_row=(3,5,120,120) passes; target_conjecture_unaffected')

    nonprimitive = (3, 6, 72, 72)
    assert descartes(nonprimitive) and gcd(*nonprimitive) == 3
    assert cube_witness(3) is None
    print('NEGATIVE_CONTROL common_gcd_3_row=(3,6,72,72) difference=3 not_cube_or_twice_cube')
    print('SCOPE finite exact checks; no universal proof or Lean settlement is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact checks for the three Chen--Chen OEIS parity-positivity conjectures.

Source: https://oeis.org/A392664/internal (m=-3), A392665 (m=-4),
A392666 (m=-5), each entry's %F recurrence and %e displayed rows.
Uses only Python 3 standard-library integer arithmetic.
Finite checks support the proposed induction; they do not prove the conjectures.
"""

from functools import lru_cache

# Mathematical data transcribed from the three entries' displayed examples.
SOURCE_ROWS = {-5: {1: [1],
      2: [1, 1],
      3: [1, -8, 1],
      4: [1, 23, 23, 1],
      5: [1, -106, -414, -106, 1],
      6: [1, 405, 5210, 5210, 405, 1],
      7: [1, -1644, -54585, -145880, -54585, -1644, 1],
      8: [1, 6547, 530721, 3079435, 3079435, 530721, 6547, 1],
      9: [1, -26222, -4966352, -55849394, -117018530, -55849394, -4966352, -26222, 1],
      10: [1, 104849, 45588716, 925915724, 3563737526, 3563737526, 925915724, 45588716, 104849, 1]},
 -4: {1: [1],
      2: [1, 1],
      3: [1, -6, 1],
      4: [1, 11, 11, 1],
      5: [1, -44, -154, -44, 1],
      6: [1, 117, 1562, 1562, 117, 1],
      7: [1, -370, -12689, -34364, -12689, -370, 1],
      8: [1, 1087, 95853, 568339, 568339, 95853, 1087, 1],
      9: [1, -3288, -695972, -8072936, -17050170, -8072936, -695972, -3288, 1],
      10: [1, 9833, 4960580, 104809652, 409138334, 409138334, 104809652, 4960580, 9833, 1]},
 -3: {1: [1],
      2: [1, 1],
      3: [1, -4, 1],
      4: [1, 3, 3, 1],
      5: [1, -14, -30, -14, 1],
      6: [1, 17, 262, 262, 17, 1],
      7: [1, -48, -1497, -4192, -1497, -48, 1],
      8: [1, 79, 8157, 50003, 50003, 8157, 79, 1],
      9: [1, -178, -42128, -514222, -1100066, -514222, -42128, -178, 1],
      10: [1, 333, 214200, 4829952, 19299834, 19299834, 4829952, 214200, 333, 1]}}


def triangle(m, last_n):
    """One-based row and column indexing, endpoints fixed before interior update."""
    rows = [[], [1]]
    for n in range(2, last_n + 1):
        row = [1] * n
        prev = rows[n - 1]
        for k in range(2, n):
            row[k - 1] = ((m * (n - k) + 1) * prev[k - 2]
                          + (m * (k - 1) + 1) * prev[k - 1])
        rows.append(row)
    return rows


@lru_cache(None)
def entry(m, n, k):
    """Independent top-down evaluation with zero exterior and endpoint priority."""
    if not 1 <= k <= n:
        return 0
    if k == 1 or k == n:
        return 1
    return ((1 + m * n - m * k) * entry(m, n - 1, k - 1)
            + (1 + m * k - m) * entry(m, n - 1, k))


def inspect(m, rows):
    r = -m
    even_entries = odd_interiors = second_column_steps = 0
    row_sum = 2
    for n in range(1, len(rows)):
        row = rows[n]
        assert len(row) == n and row[0] == row[-1] == 1, (m, n, "boundary")
        assert row == row[::-1], (m, n, "reflection")
        if n % 2 == 0:
            assert all(v > 0 for v in row), (m, n, "even positivity")
            even_entries += n
            assert sum((-1) ** j * v for j, v in enumerate(row)) == 0
        elif n >= 3:
            assert all(v < 0 for v in row[1:-1]), (m, n, "odd interior sign")
            odd_interiors += n - 2
        if n == 1:
            assert sum(row) == 1
        elif n == 2:
            assert sum(row) == row_sum
        else:
            # Source %C/%F row-sum identity, checked without division.
            row_sum *= 2 + m * (n - 2)
            assert sum(row) == row_sum, (m, n, "row sum")
        if n >= 2 and n + 2 < len(rows):
            x_n = row[1]
            forcing = r * (r - 2) * n + 2 - r * r
            assert rows[n + 2][1] == (r - 1) ** 2 * x_n + forcing
            if n >= 4:
                assert forcing >= 3 * r * r - 8 * r + 2
                assert 3 * r * r - 8 * r + 2 == 5 + 10 * (r - 3) + 3 * (r - 3) ** 2
                assert forcing >= 5
            second_column_steps += 1
    assert rows[2][1] == 1
    assert rows[4][1] == 2 * r * r - 6 * r + 3 == 2 * r * (r - 3) + 3
    return even_entries, odd_interiors, second_column_steps


def main():
    fidelity_entries = 0
    for m, displayed in SOURCE_ROWS.items():
        rows = triangle(m, 400)
        for n, expected in displayed.items():
            assert rows[n] == expected, (m, n, "source fidelity")
            for k, value in enumerate(expected, 1):
                assert entry(m, n, k) == value, (m, n, k, "independent evaluation")
                fidelity_entries += 1
        counts = inspect(m, rows)
        print(f"OEIS m={m}: displayed rows 1..10 reproduced; rows 1..400 checked")
        print(f"  positive even entries={counts[0]}; negative odd interiors={counts[1]}; "
              f"two-step column identities={counts[2]}")
    print(f"FIDELITY_OK: {fidelity_entries} displayed entries, two implementations")
    total_even = 0
    for r in range(3, 33):
        counts = inspect(-r, triangle(-r, 80))
        total_even += counts[0]
    print(f"GENERAL_ROUTE_OK: integer r=3..32, rows 1..80, {total_even} positive even entries")
    print("Reflection, exact row sums, even alternating sums, second-column identities and bounds: OK")
    print("Scope: finite exact checks only; the unbounded assertion requires the proposed induction.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite checks for Spezia's OEIS A392131/A392132 Pell circulants.

Uses only Python's standard library. Finite checks do not prove a limit.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import comb, prod

SOURCE_DET = [1, 0, 1, 9, 640, 252100, 597791376, 8252767236301, 665993168101048320, 313398861092499186897912, 859858277306561684298748471552, 13751722491910335588624048111158497953, 1281922118247965546737056182821273600000000000, 696508759202841219994984822691854206568914419368757980]
SOURCE_PER = [1, 0, 1, 9, 772, 317660, 755353604, 10245351511241, 806176992721352256, 369503021786486764387032, 989096947940706813698038153600, 15475597969523884739005106748390084993, 1415492646595459083986339197781974314319509760, 756726561532204931161254789744755112027328520513039972]
C = F(6, 5)
Q = F(12, 5)
R = F(11, 12)


def pell(limit):
    values = [0, 1]
    for _ in range(2, limit + 1):
        values.append(2 * values[-1] + values[-2])
    return values


def matrix(n, values):
    # Source Mathematica starts with RotateRight[l], then rotates each row right.
    if n == 0:
        return []
    row = values[n - 1:n] + values[:n - 1]
    rows = []
    for _ in range(n):
        rows.append(row)
        row = row[-1:] + row[:-1]
    # Independent modular-coordinate implementation, matching Mathlib.circulant.
    assert rows == [[values[n - 1 - ((i - j) % n)] for j in range(n)]
                    for i in range(n)]
    return rows


def permanent_dp(a):
    n = len(a)
    dp = [0] * (1 << n)
    dp[0] = 1
    for mask in range((1 << n) - 1):
        row = mask.bit_count()
        if not dp[mask]:
            continue
        for col in range(n):
            if not mask & (1 << col):
                dp[mask | (1 << col)] += dp[mask] * a[row][col]
    return dp[-1]


def permanent_ryser(a):
    n = len(a)
    result = 0
    for mask in range(1 << n):
        factor = prod(sum(a[i][j] for j in range(n) if mask & (1 << j))
                      for i in range(n))
        result += (-1) ** (n - mask.bit_count()) * factor
    return result


def determinant_bareiss(a):
    n = len(a)
    if n == 0:
        return 1
    b = [row[:] for row in a]
    previous, sign = 1, 1
    for k in range(n - 1):
        if b[k][k] == 0:
            pivot_row = next((i for i in range(k + 1, n) if b[i][k]), None)
            if pivot_row is None:
                return 0
            b[k], b[pivot_row] = b[pivot_row], b[k]
            sign = -sign
        pivot = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = b[i][j] * pivot - b[i][k] * b[k][j]
                assert numerator % previous == 0
                b[i][j] = numerator // previous
            b[i][k] = 0
        previous = pivot
    return sign * b[-1][-1]


def leibniz(a):
    n = len(a)
    unsigned, signed = 0, 0
    for sigma in permutations(range(n)):
        term = prod(a[i][sigma[i]] for i in range(n))
        inversions = sum(sigma[i] > sigma[j]
                         for i in range(n) for j in range(i + 1, n))
        unsigned += term
        signed += (-1) ** inversions * term
    return unsigned, signed


def cycles(n):
    # The least vertex starts the cycle: exactly one cyclic-rotation representative.
    for k in range(2, n + 1):
        for support in combinations(range(n), k):
            for tail in permutations(support[1:]):
                yield (support[0],) + tail


def cycle_bound(n):
    return F(6 * n, 11) * R ** n / (1 - R ** n)


def check_cycles(n, a, scale, permanent):
    weights = []
    rooted_by_total = {}
    for cycle in cycles(n):
        edges = list(zip(cycle, cycle[1:] + cycle[:1]))
        displacements = [(i - j) % n for i, j in edges]
        assert all(1 <= d <= n - 1 for d in displacements)
        total = sum(displacements)
        assert total > 0 and total % n == 0
        weight = prod(F(a[i][j], scale) for i, j in edges)
        assert 0 <= weight <= C ** len(cycle) / Q ** total
        weights.append(weight)
        rooted_by_total[total] = (rooted_by_total.get(total, F(0))
                                  + len(cycle) * C ** len(cycle))
    activity = sum(weights, F(0))
    assert activity <= cycle_bound(n)
    # The composition majorant includes nonsimple walks and every possible root.
    for total, rooted_weight in rooted_by_total.items():
        assert rooted_weight <= n * C * (1 + C) ** (total - 1)
    if n <= 6:
        all_subsets = prod(1 + w for w in weights)
        assert F(permanent, scale ** n) <= all_subsets
    # A rational substitute for exp(activity) when activity < 1.
    if activity < 1:
        assert F(permanent, scale ** n) <= 1 / (1 - activity)
    print(f"CYCLES n={n} count={len(weights)} activity={activity} bound={cycle_bound(n)}")


def main():
    values = pell(200)
    assert values[:7] == [0, 1, 2, 5, 12, 29, 70]
    assert matrix(4, values) == [[5, 0, 1, 2], [2, 5, 0, 1],
                                 [1, 2, 5, 0], [0, 1, 2, 5]]
    assert SOURCE_DET[4] == 640 and SOURCE_PER[4] == 772
    print("FIDELITY OEIS A392131 revision 16 determinant; A392132 revision 10 permanent")
    print("FIDELITY source example n=4: determinant=640 permanent=772")
    for j in range(2, 200):
        ratio = F(values[j + 1], values[j])
        assert Q <= ratio <= F(5, 2)
    for N in range(2, 201):
        for m in range(N):
            assert F(values[m], values[N]) <= C / Q ** (N - m)
    # Independently enumerate compositions for the exact weighted count.
    for m in range(1, 18):
        sum_by_parts = sum(F(comb(m - 1, k - 1)) * C ** k
                           for k in range(1, m + 1))
        by_cuts = F(0)
        for cuts in range(1 << (m - 1)):
            by_cuts += C ** (1 + cuts.bit_count())
        assert by_cuts == sum_by_parts == C * (1 + C) ** (m - 1)
    print("EXACT_LEMMAS Pell entry majorant: 2<=N<=200; compositions: 1<=m<=17")
    for n in range(17):
        a = matrix(n, values)
        per = permanent_dp(a)
        det = determinant_bareiss(a)
        if n < len(SOURCE_DET):
            assert det == SOURCE_DET[n]
            assert per == SOURCE_PER[n]
        if n <= 9:
            assert per == permanent_ryser(a)
        if n <= 8:
            assert (per, det) == leibniz(a)
        if n >= 2:
            assert det > 0
            scale = values[n - 1]
            assert per >= scale ** n
            assert abs(det - scale ** n) <= per - scale ** n
            ratio = F(per, det)
            print(f"MATRIX n={n} det={det} per={per} ratio={ratio}")
            if n >= 3:
                row_sum = sum((F(a[0][j], scale) for j in range(1, n)), F(0))
                assert row_sum < F(6, 7)
                for i in range(n):
                    for j in range(n):
                        if i != j:
                            d = (i - j) % n
                            assert F(a[i][j], scale) <= C / Q ** d
                if n <= 8:
                    check_cycles(n, a, scale, per)
    assert R == (1 + C) / Q and R < 1
    for n in [50, 100, 200, 500, 1000]:
        bound = cycle_bound(n)
        assert 0 < bound < F(1, 2)
        # Prose route gives |per/det - 1| <= 2B/(1-2B) on this tail.
        error_bound = 2 * bound / (1 - 2 * bound)
        assert error_bound > 0
        # Exact certificate for increasingly small tail estimates.
        if n >= 100:
            assert error_bound < F(1, 50)
        if n >= 200:
            assert error_bound < F(1, 10000)
        if n >= 500:
            assert error_bound < F(1, 10 ** 12)
        print(f"TAIL n={n} B={bound} ratio_error_bound={error_bound}")
    print("SCOPE exact source terms 0..13; matrix comparisons 2..16; cycles 3..8")
    print("SCOPE finite tests support the outline; the infinite limit is not certified here")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()

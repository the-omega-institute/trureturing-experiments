#!/usr/bin/env python3
"""Exact checks of OEIS A005121 revision 88 and Mohapatra's tail triangle.

The program uses only the Python standard library and integer arithmetic.
These finite checks do not prove the all-index conjecture.
"""
from functools import lru_cache
from math import comb, factorial

LIMIT = 50
OEIS = (
    1, 1, 4, 32, 436, 9012, 262760, 10270696, 518277560,
    32795928016, 2542945605432, 237106822506952, 26173354092593696,
    3375693096567983232, 502995942483693043200, 85750135569136650473360,
    16583651916595710735271248, 3611157196483089769387182064,
    879518067472225603327860638128,
)


def stirling_table(limit):
    table = [[0] * (limit + 1) for _ in range(limit + 1)]
    table[0][0] = 1
    for n in range(1, limit + 1):
        for k in range(1, n + 1):
            table[n][k] = k * table[n - 1][k] + table[n - 1][k - 1]
    return table


def lengyel(table, limit):
    a = [0] * (limit + 1)
    a[1] = 1
    for n in range(2, limit + 1):
        a[n] = sum(table[n][k] * a[k] for k in range(1, n))
    return a


def mohapatra(limit):
    rows = [[], [0, 1]]
    for n in range(2, limit + 1):
        row = [0] * (n + 1)
        for k in range(1, n):
            row[k] = k * sum(rows[n - 1][j] for j in range(k, n))
        row[n] = sum(row[1:n])
        rows.append(row)
    return rows


@lru_cache(None)
def partitions(n):
    """Generate each partition once by adjoining the largest element."""
    if n == 0:
        return ((),)
    answer = []
    for p in partitions(n - 1):
        answer.append(p + ((n - 1,),))
        for i in range(len(p)):
            answer.append(p[:i] + (p[i] + (n - 1,),) + p[i + 1:])
    assert len(answer) == len(set(answer))
    return tuple(answer)


def refines(p, q):
    return all(any(set(b).issubset(c) for c in q) for b in p)


def partition_chains(n):
    ps = partitions(n)

    @lru_cache(None)
    def count(p):
        if len(p) == 1:
            return 1
        return sum(count(q) for q in ps if q != p and refines(p, q))

    return count(tuple((i,) for i in range(n)))


def multiply(p, q, degree):
    return [sum(p[i] * q[n - i] for i in range(n + 1))
            for n in range(degree + 1)]


def check_truncations(rows, table):
    for degree in range(1, 13):
        columns = [[0] * (degree + 1) for _ in range(degree + 2)]
        for k in range(1, degree + 1):
            for n in range(k, degree + 1):
                columns[k][n] = rows[n][k]
        tails = [[sum(columns[j][n] for j in range(k, degree + 1))
                  for n in range(degree + 1)]
                 for k in range(degree + 2)]
        for k in range(1, degree + 1):
            for n in range(degree + 1):
                lhs = tails[k][n] - (k * tails[k][n - 1] if n else 0)
                rhs = tails[k + 1][n] + (rows[k][k] if n == k else 0)
                assert lhs == rhs, ('tail equation', degree, k, n)
        inverse = [1] + [0] * degree
        reconstructed = [0] * (degree + 1)
        for k in range(1, degree + 1):
            inverse = multiply(inverse, [k ** r for r in range(degree + 1)], degree)
            for n in range(k, degree + 1):
                assert inverse[n - k] == table[n][k], ('Stirling column', degree, k, n)
                reconstructed[n] += rows[k][k] * inverse[n - k]
        assert reconstructed == tails[1], ('tail elimination', degree)


def main():
    table = stirling_table(LIMIT)
    a = lengyel(table, LIMIT)
    rows = mohapatra(LIMIT)
    assert tuple(a[1:len(OEIS) + 1]) == OEIS
    assert tuple(rows[n][n] for n in range(1, len(OEIS) + 1)) == OEIS
    # An independent exact inclusion-exclusion evaluation of Stirling numbers.
    for n in range(1, 13):
        for k in range(1, n + 1):
            numerator = sum((-1) ** (k - j) * comb(k, j) * j ** n
                            for j in range(k + 1))
            quotient, remainder = divmod(numerator, factorial(k))
            assert remainder == 0 and quotient == table[n][k]
    for n in range(1, 7):
        assert partition_chains(n) == a[n]
        for k in range(1, n + 1):
            assert sum(len(p) == k for p in partitions(n)) == table[n][k]
    assert partition_chains(3) == 4  # OEIS's explicitly displayed example.
    for n in range(1, LIMIT + 1):
        assert rows[n][n] == a[n], ('diagonal', n)
        total = sum(rows[n][1:])
        assert total == sum(table[n][k] * rows[k][k] for k in range(1, n + 1))
        assert total == (rows[n][n] if n == 1 else 2 * rows[n][n])
    check_truncations(rows, table)
    print('Source: OEIS A005121, revision 88, timestamp 2026-11-05T17:38:02-05:00')
    print('Fidelity: all 19 published DATA terms reproduced by the Lengyel recurrence')
    print('Fidelity: the 4 strict partition chains for n=3 reproduced; chains checked through n=6')
    print('Independent Stirling inclusion-exclusion check: 1 <= k <= n <= 12')
    print('Required OEIS/triangle comparison:')
    for n in range(1, 13):
        print(f'n={n}: OEIS={OEIS[n - 1]} Lengyel={a[n]} R(n,n)={rows[n][n]}')
    print('First five triangle rows:')
    for n in range(1, 6):
        print(f'n={n}: {rows[n][1:]}')
    print('Diagonal equality and row/Stirling identity: n=1,...,50')
    print('Truncated tail equations, elimination, and Stirling columns: N=1,...,12')
    print('All computations use exact integers; finite checks are not an all-index proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

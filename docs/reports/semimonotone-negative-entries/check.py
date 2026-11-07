#!/usr/bin/env python3
"""Exact finite checks of arXiv:2603.00639v1, Remark 4.1 (first bullet).

No numerical LP or floating-point operations are used. A failure on a support S
is equivalent, by positive scaling, to feasibility of x_i >= 1 and Bx <= -1
(ordinary) or Bx <= 0 (strict). Integer Fourier--Motzkin elimination decides
these closed rational systems; back substitution supplies a Fraction witness.
The matrix predicate checks every nonempty support, not just x = (1,...,1).
"""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from itertools import combinations, product
from math import gcd, lcm


def normalize(rows):
    """Rows (a_1,...,a_d,b) encode a.x <= b; positive rescaling only."""
    result = set()
    for row in rows:
        den = lcm(*(Q(a).denominator for a in row))
        ints = tuple(int(Q(a) * den) for a in row)
        common = gcd(*ints)
        if common:
            ints = tuple(a // common for a in ints)
        if not any(ints[:-1]):
            if ints[-1] < 0:
                return None
        else:
            result.add(ints)
    return tuple(sorted(result))


def feasible(rows, dimension):
    rows = normalize(rows)
    if rows is None:
        return None
    if dimension == 0:
        return ()
    positive = [r for r in rows if r[-2] > 0]
    negative = [r for r in rows if r[-2] < 0]
    projected = [r[:-2] + (r[-1],) for r in rows if r[-2] == 0]
    for p in positive:
        for q in negative:
            projected.append(tuple(-q[-2] * a + p[-2] * b
                                   for a, b in zip(p[:-2] + (p[-1],),
                                                   q[:-2] + (q[-1],))))
    prefix = feasible(projected, dimension - 1)
    if prefix is None:
        return None
    lower, upper = [], []
    for r in rows:
        residual = r[-1] - sum(a * x for a, x in zip(r[:-2], prefix))
        if r[-2]:
            (upper if r[-2] > 0 else lower).append(residual / Q(r[-2]))
    lo = max(lower) if lower else None
    hi = min(upper) if upper else None
    if lo is not None and hi is not None:
        assert lo <= hi
    value = lo if lo is not None else (hi if hi is not None else Q(0))
    answer = prefix + (value,)
    assert all(sum(a * x for a, x in zip(r[:-1], answer)) <= r[-1] for r in rows)
    return answer


def principal(a, subset):
    return tuple(tuple(a[i][j] for j in subset) for i in subset)


@cache
def full_support_failure(a, strict):
    n = len(a)
    rows = [tuple(-int(i == j) for j in range(n)) + (-1,) for i in range(n)]
    rows += [tuple(row) + (0 if strict else -1,) for row in a]
    witness = feasible(rows, n)
    if witness is not None:
        assert all(x > 0 for x in witness)
        image = [sum(c * x for c, x in zip(row, witness)) for row in a]
        assert all(y <= 0 if strict else y < 0 for y in image)
    return witness is not None


@cache
def has_property(a, strict):
    n = len(a)
    # Every proper support is contained in a one-index deletion.
    if any(not has_property(principal(a, tuple(j for j in range(n) if j != i)), strict)
           for i in range(n)):
        return False
    return n == 0 or not full_support_failure(a, strict)


def exact_orders(a, strict):
    n = len(a)
    profile = [(len(s), has_property(principal(a, s), strict))
               for m in range(1, n + 1) for s in combinations(range(n), m)]
    return [k for k in range(1, n)
            if all(good for m, good in profile if m == n - k)
            and all(not good for m, good in profile if m > n - k)]


def strongly_connected(a, vertices):
    for start in vertices:
        reached = {start}
        frontier = [start]
        while frontier:
            i = frontier.pop()
            for j in vertices:
                if i != j and a[i][j] < 0 and j not in reached:
                    reached.add(j)
                    frontier.append(j)
        if len(reached) != len(vertices):
            return False
    return True


def check_member(a, k, strict):
    n = len(a)
    assert all(a[i][i] > 0 if strict else a[i][i] >= 0 for i in range(n))
    assert all(sum(a[i][j] < 0 for j in range(n)) >= k for i in range(n))
    assert all(sum(a[i][j] < 0 for i in range(n)) >= k for j in range(n))
    # A separate exact test of the outline's stronger graph consequence.
    for m in range(k):
        for deleted in combinations(range(n), m):
            remaining = tuple(i for i in range(n) if i not in deleted)
            assert len(remaining) >= 2 and strongly_connected(a, remaining)


def inverse_and_determinant(a):
    n = len(a)
    rows = [[Q(x) for x in a[i]] + [Q(i == j) for j in range(n)] for i in range(n)]
    det = Q(1)
    for j in range(n):
        pivot = next(i for i in range(j, n) if rows[i][j])
        if pivot != j:
            rows[j], rows[pivot] = rows[pivot], rows[j]
            det = -det
        v = rows[j][j]
        det *= v
        rows[j] = [x / v for x in rows[j]]
        for i in range(n):
            if i != j:
                v = rows[i][j]
                rows[i] = [x - v * y for x, y in zip(rows[i], rows[j])]
    return tuple(tuple(row[n:]) for row in rows), det


def fidelity():
    # Example 3.2 (label example:4.2 in main.tex), including the displayed inverse.
    a = ((1, Q(1, 2), -1, -1), (Q(1, 2), 1, -1, -1),
         (-1, -1, 1, 0), (-1, -1, 0, 1))
    counts = [sum(has_property(principal(a, s), False)
                  for s in combinations(range(4), m)) for m in range(1, 5)]
    assert counts == [4, 6, 0, 0]
    assert exact_orders(a, False) == [2]
    inverse, det = inverse_and_determinant(a)
    expected = tuple(tuple(Q(x, 5) for x in row)
                     for row in ((4, -6, -2, -2), (-6, 4, -2, -2),
                                 (-2, -2, 1, -4), (-2, -2, -4, 1)))
    assert inverse == expected and det == Q(-5, 4)
    assert all(sum(a[i][k] * inverse[k][j] for k in range(4)) == int(i == j)
               for i in range(4) for j in range(4))
    check_member(a, 2, False)
    print('FIDELITY Example 3.2: good principal blocks by size =', counts,
          '; exact ordinary orders = [2]; det =', det,
          '; displayed inverse reproduced')


def controls():
    assert feasible([(2, 1), (-3, -1)], 1) == (Q(1, 3),)
    assert feasible([(2, 0), (-1, -1)], 1) is None
    assert feasible([(-1,)], 0) is None
    boundary = ((1, -1), (-1, 1))
    assert has_property(boundary, False) and not has_property(boundary, True)
    # Independent closed-form decision for all 81 matrices of order two.
    for entries in product((-1, 0, 1), repeat=4):
        a, b, c, d = entries
        matrix = ((a, b), (c, d))
        determinant = a * d - b * c
        ordinary = a >= 0 and d >= 0 and (b >= 0 or c >= 0 or determinant >= 0)
        strict = a > 0 and d > 0 and (b >= 0 or c >= 0 or determinant > 0)
        assert has_property(matrix, False) == ordinary
        assert has_property(matrix, True) == strict
    print('CONTROLS: exact feasible/infeasible systems; strict zero boundary; all 81 order-two decisions agree with the determinant criterion')


def census(n, matrices, label):
    counts = {False: Counter(), True: Counter()}
    total = 0
    for a in matrices:
        total += 1
        for strict in (False, True):
            for k in exact_orders(a, strict):
                check_member(a, k, strict)
                counts[strict][k] += 1
    print(label, 'matrices =', total)
    for strict in (False, True):
        print(' ', 'strict' if strict else 'ordinary',
              'exact-order membership counts =', {k: counts[strict][k] for k in range(1, n)})
    assert all(sum(counts[strict].values()) > 0 for strict in (False, True))


def main():
    controls()
    fidelity()
    matrices3 = (tuple(tuple(entries[3*i:3*i+3]) for i in range(3))
                 for entries in product((-1, 0, 1), repeat=9))
    census(3, matrices3, 'EXHAUSTIVE n=3, every entry in {-1,0,1}:')
    positions = [(i, j) for i in range(4) for j in range(4) if i != j]
    def matrices4():
        for entries in product((-1, 0), repeat=12):
            values = dict(zip(positions, entries))
            yield tuple(tuple(1 if i == j else values[i, j] for j in range(4)) for i in range(4))
    census(4, matrices4(), 'EXHAUSTIVE n=4, diagonal 1 and every off-diagonal entry in {-1,0}:')
    print('Every qualifying matrix passed the row/column bound and all deletions of fewer than k vertices.')
    print('Finite tests only: no assertion about untested real matrices or unbounded dimensions is certified by this program.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

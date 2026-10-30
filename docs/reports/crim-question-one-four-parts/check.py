#!/usr/bin/env python3
"""Exact, independent checks for arXiv:2606.16828v1, Question 1.

No Lean proof or unbounded computational certification is claimed.
The evaluator uses the source's conjugate/delete/conjugate definition.
The candidate formula and strategy never enter the Grundy recursion.
"""

from collections import Counter
from functools import cache
from itertools import combinations_with_replacement
from math import comb


WIDTH = 30


def conjugate(position):
    return tuple(sum(part >= j for part in position)
                 for j in range(1, max(position, default=0) + 1))


def row_delete(position, i):
    return position[:i] + position[i + 1:]


def source_moves(position):
    transposed = conjugate(position)
    return {row_delete(position, i) for i in range(len(position))} | {
        conjugate(row_delete(transposed, j)) for j in range(len(transposed))}


def column_height_delete(position, t):
    assert 1 <= t <= len(position)
    below = position[t] if t < len(position) else 0
    assert position[t - 1] > below, (position, t)
    return tuple(x - 1 if i < t else x for i, x in enumerate(position)
                 if i >= t or x > 1)


def height_moves(position):
    return {row_delete(position, i) for i in range(len(position))} | {
        column_height_delete(position, t)
        for t in range(1, len(position) + 1)
        if position[t - 1] > (position[t] if t < len(position) else 0)}


@cache
def grundy(position):
    followers = source_moves(position)
    assert all(sum(q) < sum(position) for q in followers)
    used = {grundy(q) for q in followers}
    answer = 0
    while answer in used:
        answer += 1
    return answer


def four_part_criterion(a, b, c, d):
    return a % 2 == 0 and (
        (d >= 2 and not (b % 2 == 1 and c % 2 == 1))
        or (d == 1 and c == 1 and (b % 2 == 0 or b == 1))
        or (d == 1 and c >= 2 and b % 2 == 1)
        or (d == 1 and c == 2 and (b >= 3 or a == 2)))


def pstar(position):
    match position:
        case ():
            return True
        case (a,):
            return False
        case (a, b):
            return a % 2 == 0
        case (a, b, c):
            return a > 1 and a % 2 == 1 and b % 2 == 1 and (c > 1 or b == 1)
        case (a, b, c, d):
            return four_part_criterion(a, b, c, d)
        case _:
            raise AssertionError(position)


def four_part_strategy(position):
    a, b, c, d = position
    assert not pstar(position)
    if a % 2 == 0:
        if d > 1:
            return "01 even a, d>1: R1", row_delete(position, 0)
        if c == 1:
            return "02 even a, d=c=1: R1", row_delete(position, 0)
        return "03 even a, remaining: C4", column_height_delete(position, 4)
    if a == 1:
        return "04 all ones: C4", column_height_delete(position, 4)
    if c == 1 and b == 1:
        return "05 odd a>1, c=b=1: C1", column_height_delete(position, 1)
    if c == 1:
        return "06 odd a, c=1<b: R2", row_delete(position, 1)
    if b % 2 == 1:
        return "07 odd a,b, c>1: R4", row_delete(position, 3)
    if d > 1:
        if c % 2 == 1:
            return "08 odd a, even b, d>1, odd c: R2", row_delete(position, 1)
        return "09 odd a, even b,c, d>1: C1", column_height_delete(position, 1)
    if c % 2 == 1:
        return "10 odd a, even b, d=1, odd c: C2", column_height_delete(position, 2)
    if c >= 4:
        return "11 odd a, even b,c>=4, d=1: C3", column_height_delete(position, 3)
    if b == 2:
        return "12 odd a, b=c=2, d=1: C3", column_height_delete(position, 3)
    return "13 odd a, even b>c=2, d=1: C1", column_height_delete(position, 1)


def main():
    print("Source: Basic-Gottlieb-Krnc, arXiv:2606.16828v1, Question 1")
    print("Arithmetic: Python integers only; no random sampling")
    # Fidelity to Proposition 10: all positive two-part inputs in this range.
    fidelity_count = 0
    for a in range(1, WIDTH + 1):
        for b in range(1, a + 1):
            assert (grundy((a, b)) == 0) == (a % 2 == 0), (a, b)
            fidelity_count += 1
    print(f"FIDELITY Proposition 10: {fidelity_count} two-part inputs, width <= {WIDTH}: OK")
    assert grundy(()) == 0 and grundy((1, 1, 1)) == 1
    print("BOUNDARY empty SG=0; [1,1,1] SG=1 (Theorem 11 needs a>1): OK")
    counts, losing, strategies = Counter(), Counter(), Counter()
    for rows in range(5):
        for ascending in combinations_with_replacement(range(1, WIDTH + 1), rows):
            p = tuple(reversed(ascending))
            assert conjugate(conjugate(p)) == p, p
            followers = source_moves(p)
            assert followers == height_moves(p), p
            actual = grundy(p) == 0
            proposed = pstar(p)
            assert actual == proposed, (p, grundy(p), proposed)
            assert (not any(pstar(q) for q in followers)) == proposed, p
            if rows == 4 and not proposed:
                case, q = four_part_strategy(p)
                assert q in followers and pstar(q) and grundy(q) == 0, (case, p, q)
                strategies[case] += 1
            counts[rows] += 1
            losing[rows] += actual
    assert counts[4] == comb(WIDTH + 3, 4) == 40920
    assert sum(counts.values()) == comb(WIDTH + 4, 4) == 46376
    assert len(strategies) == 13
    for rows in range(5):
        print(f"CLASSIFICATION rows={rows}: checked={counts[rows]}, P={losing[rows]}: OK")
    print("Source moves = row/column-height moves; strict cell descent; P*/N* transitions: OK")
    for case, count in sorted(strategies.items()):
        print(f"STRATEGY {case}: checked={count}: OK")
    for p, expected in [((4, 3, 3, 2), False), ((6, 4, 3, 2), True)]:
        assert (grundy(p) == 0) == expected
        print(f"EXAMPLE {list(p)}: SG={grundy(p)}: {'P' if expected else 'N'}")
    print("SCOPE: width <= 30, at most four positive parts; no unbounded proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

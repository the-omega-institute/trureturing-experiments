#!/usr/bin/env python3
"""Exact finite fidelity and identity checks for Hanna's A381356 observation.

Sources: OEIS A381356 revision 11, A381587 revision 36, A306346 revision 36.
All arithmetic is on integers. Finite checks do not prove the infinite claim.
Run with Python 3; no packages or network access are needed.
"""
from itertools import groupby, product

A381356_TERMS = [1, 3, 1, 3, 1, 3, 1, 1, 1, 5, 1, 1, 1, 1, 1, 7, 1, 3, 1, 1, 1, 1, 1, 7, 1, 1, 1, 5, 1, 3, 1, 1, 1, 1, 1, 7, 1, 1, 1, 3, 1, 5, 1, 1, 1, 5, 1, 3, 1, 1, 1, 1, 1, 7, 1, 1, 1, 3, 1, 1, 1, 3, 1, 5, 1, 1, 1, 3, 1, 5, 1, 1, 1, 5, 1, 3, 1, 1, 1, 1, 1, 7, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 3, 1, 5]
A381587_TERMS = [1, 1, 2, 1, 2, 1, 1, 1, 2, 1, 3, 1, 1, 1, 2, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2, 1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2, 1, 3, 1, 3, 1, 3, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3]
A306346_TERMS = [0, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 3, 1, 7, 1, 1, 1, 1, 1, 5, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 3, 1, 5, 1, 1, 1, 7, 1, 1, 1, 1, 1, 3, 1, 7, 1]
PRINTED_W_ROWS = [[1], [1], [2], [1, 2], [1, 1, 1, 2], [1, 3, 1, 1, 1, 2], [1, 3, 1, 1, 1, 3, 1, 1, 1, 2], [1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2], [1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2], [1, 3, 1, 3, 1, 3, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2], [1, 3, 1, 3, 1, 3, 1, 1, 1, 5, 1, 1, 1, 1, 1, 7, 1, 3, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 2]]
PRINTED_S_ROWS = [[0], [0, 1], [0, 1, 1, 1], [0, 1, 1, 1, 3, 1], [0, 1, 1, 1, 3, 1, 1, 1, 3, 1], [0, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 3, 1], [0, 1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 1, 3, 1, 1, 1, 1, 1, 3, 1, 3, 1, 3, 1]]

def runs(word):
    """Run lengths via maximal equal-value groups."""
    return [sum(1 for _ in group) for _, group in groupby(word)]


def runs_by_endpoints(word):
    """Independent implementation via positions where adjacent values differ."""
    if not word:
        return []
    endpoints = [0]
    endpoints += [i for i in range(1, len(word)) if word[i] != word[i - 1]]
    endpoints.append(len(word))
    return [right - left for left, right in zip(endpoints, endpoints[1:])]


def append_orbit(seed, count):
    rows = [[seed]]
    for _ in range(count - 1):
        word = rows[-1]
        rows.append(word + runs(word[::-1]))
    return rows


def flattened_triangle(count):
    rows = [[1]]
    flattened = [1]
    for _ in range(count - 1):
        row = runs(flattened[::-1])
        rows.append(row)
        flattened += row
    return rows


def row_only_triangle(count):
    rows = [[1], [1], [2]]
    for _ in range(3, count):
        word = rows[-1]
        rows.append(runs(word[::-1]) + word)
    return rows[:count]


def main():
    words = [list(w) for size in range(9) for w in product(range(3), repeat=size)]
    for word in words:
        lengths = runs(word)
        assert lengths == runs_by_endpoints(word)
        assert sum(lengths) == len(word)
        assert all(length > 0 for length in lengths)
        assert runs(word[::-1]) == lengths[::-1]
    small = [word for word in words if len(word) <= 4]
    boundary_checks = 0
    for x in small:
        for y in small:
            if x and y and x[-1] != y[0]:
                assert runs(x + y) == runs(x) + runs(y)
                boundary_checks += 1
    print(f"RUNS: {len(words)} words over {{0,1,2}} of lengths 0..8; independent endpoints agree")
    print(f"UNEQUAL_BOUNDARY: {boundary_checks} eligible pairs of nonempty words of lengths <= 4")

    w = flattened_triangle(27)
    alternate = row_only_triangle(27)
    assert w == alternate
    assert w[:len(PRINTED_W_ROWS)] == PRINTED_W_ROWS
    s = append_orbit(0, 24)
    u = append_orbit(2, 24)
    assert s[:len(PRINTED_S_ROWS)] == PRINTED_S_ROWS
    flattened20 = [entry for row in w[:20] for entry in row]
    assert len(flattened20) == 9259
    assert flattened20[:len(A381587_TERMS)] == A381587_TERMS
    assert s[-1][:len(A306346_TERMS)] == A306346_TERMS
    assert runs(s[-1])[:len(A381356_TERMS)] == A381356_TERMS
    print(f"FIDELITY: all {len(PRINTED_W_ROWS)} complete displayed W rows and {len(PRINTED_S_ROWS)} displayed S rows agree")
    print(f"FIDELITY: A381587 flattened first 20 rows has {len(flattened20)} terms as its source states")
    print(f"FIDELITY: printed prefixes A381356={len(A381356_TERMS)}, A381587={len(A381587_TERMS)}, A306346={len(A306346_TERMS)} agree")

    for i, (sn, un) in enumerate(zip(s, u), start=1):
        assert un == [2] + sn[1:]
        assert runs(un) == runs(sn)
        assert w[i + 1] == un[::-1]  # W_(i+2) = reverse U_i
        assert w[i + 2] == runs(sn) + un[::-1]  # W_(i+3)
        assert len(sn) >= i
        assert sn[0] == 0 and all(value > 0 for value in sn[1:])
        if i >= 2:
            assert sn[-1] == 1
    for row in w[2:]:
        assert row[-1] == 2
        assert len(row) == 1 or row[-2] != 2
    for row in w[3:]:
        assert row[0] == 1
    print("IDENTITY: flattened and row-only definitions agree through W_27")
    print("IDENTITY: seed replacement, reversal conjugacy, and RUNS-prefix formula hold for n=1..24")
    print(f"S_23_LENGTH: {len(s[22])}")

    for i, sn in enumerate(s):
        # The last finite run may extend later; every preceding run is complete.
        stable = runs(sn)[:-1]
        for later in s[i:]:
            assert runs(later)[:len(stable)] == stable
        for later_row in w[i + 3:]:
            assert later_row[:len(stable)] == stable
    stable25 = runs(s[21])[:-1]
    assert len(stable25) == 10797
    assert w[24][:len(stable25)] == stable25
    print("FINITE_STABILIZATION: every completed run at each tested stage persists in all later tested words and rows")
    print(f"W_25_CERTIFIED_FINITE_PREFIX: {len(stable25)} entries, matching the source's stated prefix length")
    print("SCOPE: these checks are finite; infinite run finiteness and row convergence remain proof obligations")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

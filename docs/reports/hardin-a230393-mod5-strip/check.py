#!/usr/bin/env python3
"""Exact finite checks of the literal OEIS A230393 arrays and proposed transfer.

Run with Python >= 3.12 and SymPy 1.14.0. No network access is used.
Finite agreement is evidence for a proposed proof, not an all-height proof.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import hashlib
import sys
import sympy as sp

assert sys.version_info >= (3, 12)
ROOT = Path(__file__).resolve().parent
ROWS = list(product(range(5), repeat=2))
STATES = list(product(range(5), (False, True), (False, True)))
SOURCE_DATA = [
    0, 0, 0, 0, 2, 8, 32, 120, 430, 1508, 5216, 17872, 60866, 206448,
    698312, 2357640, 7949742, 26782492, 90175912, 303496048, 1021160706,
    3435198296, 11554530832, 38860960600, 130691626574, 439504634068,
    1477972726608, 4970048569744, 16712785339938, 56199558987264,
]
T = sp.Matrix([
    [2, 1, 0, 0, 0],
    [2, 0, 1, 0, 0],
    [2, 0, 0, 1, 0],
    [1, 0, 0, 0, 1],
    [8, 2, 2, 2, 2],
])
# b-a modulo 5; pending flag for left cell; pending flag for right cell.
CLASSES = [(1, False, True), (2, True, False), (0, False, True),
           (2, False, True), (1, False, False)]


def row_ok(above, row, below):
    """Each cell sees its successor horizontally, above, or below."""
    return all(
        row[1-j] == (row[j]+1) % 5
        or (above is not None and above[j] == (row[j]+1) % 5)
        or (below is not None and below[j] == (row[j]+1) % 5)
        for j in range(2)
    )


def literal_grid_ok(grid):
    """Independent coordinate-based interpretation of the source definition."""
    n = len(grid)
    return bool(n) and grid[0][0] == 0 and all(
        any(0 <= i+di < n and 0 <= j+dj < 2
            and grid[i+di][j+dj] == (grid[i][j]+1) % 5
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)))
        for i in range(n) for j in range(2)
    )


def brute_count(n):
    answer = 0
    witnesses = []
    for tail in product(range(5), repeat=2*n-1):
        cells = (0,) + tail
        grid = list(zip(cells[::2], cells[1::2]))
        if literal_grid_ok(grid):
            answer += 1
            witnesses.append(grid)
    return answer, witnesses


def pair_counts(max_n):
    """Literal row-pair DP, with no normalized states or pending flags."""
    counts = {1: sum(row_ok(None, (0, b), None) for b in range(5))}
    pairs = Counter({((0, b), below): 1
                     for b in range(5) for below in ROWS
                     if row_ok(None, (0, b), below)})
    for n in range(2, max_n+1):
        counts[n] = sum(k for (above, row), k in pairs.items()
                        if row_ok(above, row, None))
        nxt = Counter()
        for (above, row), k in pairs.items():
            for below in ROWS:
                if row_ok(above, row, below):
                    nxt[(row, below)] += k
        pairs = nxt
    return counts


def transitions(state):
    b, p, q = state
    for c, d in ROWS:
        if p and c != 1:
            continue
        if q and d != (b+1) % 5:
            continue
        yield ((d-c) % 5,
               0 != (c+1) % 5 and d != (c+1) % 5,
               b != (d+1) % 5 and c != (d+1) % 5)


def reflection(state):
    b, p, q = state
    return ((-b) % 5, q, p)


def closure(seed, expand):
    current = set(seed)
    while True:
        nxt = current | set(expand(current))
        if nxt == current:
            return current
        current = nxt


def main():
    bfile = ROOT / 'b230393.txt'
    data = {int(a): int(b) for line in bfile.read_text().splitlines()
            if line.strip() and not line.startswith('#')
            for a, b in [line.split()]}
    assert sorted(data) == list(range(1, 211))
    assert [data[n] for n in range(1, 31)] == SOURCE_DATA
    print('SOURCE: OEIS A230393 revision 12; offset 1,5; Barker formula Mar 17 2018')
    print('FIDELITY: 30 displayed terms equal the b-file prefix')
    print('b230393.txt sha256:', hashlib.sha256(bfile.read_bytes()).hexdigest())
    brute = {}
    for n in range(1, 6):
        brute[n], witnesses = brute_count(n)
        assert brute[n] == data[n]
        if n == 5:
            print('height-5 literal arrays:', witnesses)
    print('EXHAUSTIVE: all normalized arrays at heights 1..5:', list(brute.values()))
    pairs = pair_counts(14)
    assert all(pairs[n] == data[n] for n in pairs)
    print('LITERAL ROW-PAIR DP: heights 1..14 agree with the OEIS b-file')

    trans = {s: Counter(transitions(s)) for s in STATES}
    initial = Counter((b, b != 1, 0 != (b+1) % 5) for b in range(5))
    accepting = {s for s in STATES if not s[1] and not s[2]}
    reachable = closure(initial, lambda r: (t for s in r for t in trans[s]))
    productive = closure(accepting, lambda r: (s for s in STATES if set(trans[s]) & r))
    live = reachable & productive
    assert (len(reachable), len(productive), len(live)) == (13, 17, 10)
    members = [{cl, reflection(cl)} for cl in CLASSES]
    assert set.union(*members) == live and all(len(m) == 2 for m in members)
    assert [sum(initial[s] for s in m) for m in members] == [2, 0, 0, 0, 0]
    assert [[s in accepting for s in m] for m in members] == [[False]*2]*4+[[True]*2]
    for i, m in enumerate(members):
        for s in m:
            row = [sum(trans[s][t] for t in dest) for dest in members]
            assert row == list(T.row(i)), (s, row)
    print('STATE CHECK: 20 states; 13 reachable; 17 productive; 10 live')
    print('CLASS ORDER (delta,left-pending,right-pending):', CLASSES)
    print('QUOTIENT:', T.tolist())
    print('BOUNDARIES: initial=[2,0,0,0,0]; terminal=[0,0,0,0,1]')

    v = sp.Matrix([[2, 0, 0, 0, 0]])
    full = initial
    for n in range(1, 211):
        quotient_count = int(v[4])
        full_count = sum(k for s, k in full.items() if s in accepting)
        assert quotient_count == full_count == data[n], n
        nxt = Counter()
        for s, k in full.items():
            for t, multiplicity in trans[s].items():
                nxt[t] += k*multiplicity
        full = nxt
        v = v*T
    print('ALL 210 TERMS: unpruned 20-state DP = 5-class transfer = OEIS b-file')

    x = sp.Symbol('x')
    M = sp.eye(5)-x*T
    denominator = 1-4*x+4*x**3+9*x**4+6*x**5
    assert sp.expand(M.det()) == denominator
    assert sp.expand(M.adjugate()[0, 4]) == x**4
    assert T**5-4*T**4+4*T**2+9*T+6*sp.eye(5) == sp.zeros(5)
    print('det(I-xT):', sp.expand(M.det()))
    print('adj(I-xT)[1,5]:', sp.expand(M.adjugate()[0, 4]))
    print('MATRIX IDENTITY: T^5 - 4T^4 + 4T^2 + 9T + 6I = 0')
    # Put the terminal indicator into the all-ones convention of scalarRelationShift.
    P = sp.eye(5)
    for i in range(4):
        P[i, 4] = -1
    P_inv = P.inv()
    B = P_inv*T*P
    assert P*sp.ones(5, 1) == sp.eye(5)[:, 4]
    assert 2*sp.eye(5)[0, :]*P == sp.Matrix([[2, 0, 0, 0, -2]])
    assert B**5-4*B**4+4*B**2+9*B+6*sp.eye(5) == sp.zeros(5)
    print('SIGNED BASIS TRANSFER B=P^-1 T P:', B.tolist())
    print('SIGNED BASIS: P*1=e5; initial=2(e1-e5); verified exact inverse and power identity')

    coefficients = {0: 0, **data}
    for n in range(211):
        c = sum(k*coefficients.get(n-j, 0) for j, k in [(0, 1),(1, -4),(3, 4),(4, 9),(5, 6)] if n >= j)
        assert c == (2 if n == 5 else 0), n
    for n in range(6, 211):
        assert data[n] == 4*data[n-1]-4*data[n-3]-9*data[n-4]-6*data[n-5]
    print('GF CONVOLUTION: coefficients 0..210 equal 2*x^5')
    print('RECURRENCE: all indices 6..210 checked with exact integers')
    print('LIMITATION: finite exact evidence; no Lean proof of the all-height claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

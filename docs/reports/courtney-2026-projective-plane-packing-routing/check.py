#!/usr/bin/env python3
"""Exact finite tests for Courtney's whole-line packing model (Python 3 stdlib)."""
from collections import Counter, deque
from itertools import combinations, permutations, product
from math import factorial
from random import Random


def plane(q):
    """PG(2,q) for the prime q used below, with canonical homogeneous coordinates."""
    assert q >= 2 and all(q % d for d in range(2, q))
    vectors = []
    for v in product(range(q), repeat=3):
        if any(v):
            first = next(x for x in v if x)
            if first == 1:
                vectors.append(v)
    points = tuple(vectors)
    lines = tuple(frozenset(i for i, v in enumerate(points)
                            if sum(a * b for a, b in zip(v, w)) % q == 0)
                  for w in vectors)
    n = q*q + q + 1
    assert len(points) == len(lines) == n
    assert len(set(lines)) == n
    assert all(len(l) == q+1 for l in lines)
    assert all(sum(i in l for l in lines) == q+1 for i in range(n))
    assert all(len(a & b) == 1 for a, b in combinations(lines, 2))
    assert all(sum(i in l and j in l for l in lines) == 1
               for i, j in combinations(range(n), 2))
    return lines


def fidelity(q, lines):
    """Source Section 2 adjacency, Proposition 2.2, and Definition 2.1, exactly."""
    n = len(lines)
    A = [[0 if i == j else sum(i in l and j in l for l in lines)
          for j in range(n)] for i in range(n)]
    assert A == [[int(i != j) for j in range(n)] for i in range(n)]
    assert all(sum(row) == (q+1)*q for row in A)
    # The all-one vector and e_i-e_(n-1) form an exact eigenbasis.
    for i in range(n-1):
        v = [int(j == i) - int(j == n-1) for j in range(n)]
        assert [sum(A[j][k]*v[k] for k in range(n)) for j in range(n)] == [-x for x in v]
    # |-1-(r-2)| <= 2 sqrt((d-1)(r-1)); both sides nonnegative, square exactly.
    assert q*q <= 4*q*q
    print(f'FIDELITY q={q} N={n} d=r={q+1} nu=1 A=J-I eigenvalues={q*(q+1)},-1 multiplicities=1,{n-1}')


def line_generators(n, lines):
    result = set()
    for line in lines:
        l = sorted(line)
        for p in permutations(l):
            g = list(range(n))
            for x, y in zip(l, p):
                g[x] = y
            if g != list(range(n)):
                result.add(tuple(g))
    return sorted(result)


def fano_bfs(lines):
    n = len(lines)
    gens = line_generators(n, lines)
    start = bytes(range(n))
    dist = {start: 0}
    queue = deque([start])
    while queue:
        s = queue.popleft()
        for g in gens:
            t = bytes(s[i] for i in g)
            if t not in dist:
                dist[t] = dist[s]+1
                queue.append(t)
    assert len(dist) == factorial(n)
    diameter = max(dist.values())
    assert diameter <= 4*2+1 <= 5*2
    assert diameter >= 3
    print(f'FANO_BFS states={len(dist)} generators={len(gens)} transitions={len(dist)*len(gens)} diameter={diameter} histogram={dict(sorted(Counter(dist.values()).items()))}')
    return dist


def route(q, lines, initial):
    """Boundary correction, regular bipartite edge coloring, serialized R-C-R."""
    n = len(lines)
    state = list(initial)
    assert sorted(state) == list(range(n))
    B = lines[0]
    moves = 0

    def perform(line, assignment):
        nonlocal moves
        assert line in lines
        assert set(assignment) <= line
        before = state[:]
        after = state[:]
        for p, t in assignment.items():
            after[p] = t
        assert sorted(after[p] for p in line) == sorted(before[p] for p in line)
        assert all(before[p] == after[p] for p in range(n) if p not in line)
        if before != after:
            moves += 1
            state[:] = after

    fixed = set()
    for b in sorted(B):
        if state[b] != b:
            x = state.index(b)
            line, = [l for l in lines if b in l and x in l]
            perform(line, {b: state[x], x: state[b]})
        fixed.add(b)
        assert all(state[c] == c for c in fixed)
    boundary_moves = moves
    assert boundary_moves <= q+1
    assert all(state[c] == c for c in B)
    assert set(state[p] for p in range(n) if p not in B) == set(range(n))-B

    a, b = sorted(B)[:2]
    rows = [l for l in lines if a in l and l != B]
    cols = [l for l in lines if b in l and l != B]
    assert len(rows) == len(cols) == q
    grid = []
    for row in rows:
        row_cells = []
        for col in cols:
            cell, = row & col
            assert cell not in B
            row_cells.append(cell)
        grid.append(row_cells)
    assert set(sum(grid, [])) == set(range(n))-B
    assert len(set(sum(grid, []))) == q*q
    destrow = {p: i for i in range(q) for p in grid[i]}
    edges = {(i, k): [] for i in range(q) for k in range(q)}
    for i in range(q):
        for p in grid[i]:
            edges[i, destrow[state[p]]].append(state[p])
    colors = {}
    for color in range(q):
        degree = q-color
        assert all(sum(len(edges[i, k]) for k in range(q)) == degree for i in range(q))
        assert all(sum(len(edges[i, k]) for i in range(q)) == degree for k in range(q))
        # Check Hall's condition on every left set, not just the matching result.
        for mask in range(1 << q):
            S = [i for i in range(q) if (mask >> i) & 1]
            neighbors = {k for i in S for k in range(q) if edges[i, k]}
            assert len(S) <= len(neighbors)
        owner = {}

        def augment(i, seen):
            for k in range(q):
                if edges[i, k] and k not in seen:
                    seen.add(k)
                    if k not in owner or augment(owner[k], seen):
                        owner[k] = i
                        return True
            return False

        for i in range(q):
            assert augment(i, set())
        assert len(owner) == q
        for k, i in owner.items():
            colors[edges[i, k].pop()] = color
    assert all(not e for e in edges.values())
    for i, row in enumerate(rows):
        assignment = {grid[i][colors[state[p]]]: state[p] for p in grid[i]}
        assert len(assignment) == q
        perform(row, assignment)
        assert all(state[c] == c for c in B)
    for j, col in enumerate(cols):
        tokens = [state[grid[i][j]] for i in range(q)]
        assignment = {grid[destrow[t]][j]: t for t in tokens}
        assert len(assignment) == q
        perform(col, assignment)
        assert all(state[c] == c for c in B)
    for i, row in enumerate(rows):
        perform(row, {p: p for p in grid[i]})
        assert all(state[c] == c for c in B)
    assert state == list(range(n))
    assert moves-boundary_moves <= 3*q
    assert moves <= 4*q+1 <= 5*q
    return moves


def main():
    rng = Random(260502498)
    for q in [2, 3, 5, 7]:
        lines = plane(q)
        fidelity(q, lines)
        n = len(lines)
        if q == 2:
            cases = fano_bfs(lines)
            costs = [route(q, lines, s) for s in cases]
            print(f'CONSTRUCTION q=2 exhaustive_cases={len(costs)} max_moves={max(costs)} bound={4*q+1}')
        else:
            cases = {tuple(range(n)), tuple(reversed(range(n)))}
            cases.update(tuple((i+k) % n for i in range(n)) for k in range(n))
            for i, j in combinations(range(n), 2):
                s = list(range(n)); s[i], s[j] = s[j], s[i]; cases.add(tuple(s))
            for _ in range(1000):
                s = list(range(n)); rng.shuffle(s); cases.add(tuple(s))
            costs = [route(q, lines, s) for s in sorted(cases)]
            print(f'CONSTRUCTION q={q} tested_cases={len(cases)} max_moves={max(costs)} bound={4*q+1} exhaustive=False')
        if q == 3:
            generators = len(line_generators(n, lines))
            assert generators == 299
            states = factorial(n)
            print(f'PG23_BFS_NOT_RUN states={states} generators={generators} full_transitions={states*generators} one_byte_distance_storage={states} eight_byte_queue_storage={8*states}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

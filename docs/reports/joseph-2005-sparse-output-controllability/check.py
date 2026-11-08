#!/usr/bin/env python3
"""Exact finite checks for Joseph's sparse output controllability question.

Uses only Python's standard library. All linear algebra is over Fraction.
Finite checks support, and do not prove, the universal rank-cut theorem.
"""
from fractions import Fraction as Q
from itertools import combinations, product
import json
from math import comb
from random import Random


def mat(rows):
    return tuple(tuple(Q(x) for x in row) for row in rows)


def identity(n):
    return mat([[int(i == j) for j in range(n)] for i in range(n)])


def mm(a, b):
    return tuple(tuple(sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
                       for j in range(len(b[0]))) for i in range(len(a)))


def mv(a, v):
    return tuple(sum((x * y for x, y in zip(row, v)), Q(0)) for row in a)


def columns(a):
    return [tuple(row[j] for row in a) for j in range(len(a[0]))]


def rank(cols, n):
    """Row elimination of a matrix supplied as a list of columns."""
    a = [[Q(v[i]) for v in cols] for i in range(n)]
    r = 0
    for j in range(len(cols)):
        pivot = next((i for i in range(r, n) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        z = a[r][j]
        a[r] = [x / z for x in a[r]]
        for i in range(r + 1, n):
            z = a[i][j]
            if z:
                a[i] = [x - z * y for x, y in zip(a[i], a[r])]
        r += 1
        if r == n:
            break
    return r


def det(rows):
    """Independent Leibniz recursion, used to cross-check small ranks."""
    if not rows:
        return Q(1)
    if len(rows) == 1:
        return rows[0][0]
    return sum(((-1) ** j * rows[0][j] *
                det([row[:j] + row[j+1:] for row in rows[1:]])
                for j in range(len(rows))), Q(0))


def rank_by_minors(cols, n):
    for k in range(min(n, len(cols)), 0, -1):
        for rs in combinations(range(n), k):
            for cs in combinations(range(len(cols)), k):
                if det([[cols[j][i] for j in cs] for i in rs]):
                    return k
    return 0


def solve(cols, target):
    n = len(target)
    assert len(cols) == n
    a = [[cols[j][i] for j in range(n)] + [target[i]] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if a[i][j])
        a[j], a[pivot] = a[pivot], a[j]
        z = a[j][j]
        a[j] = [x / z for x in a[j]]
        for i in range(n):
            if i != j:
                z = a[i][j]
                a[i] = [x - z * y for x, y in zip(a[i], a[j])]
    return [a[i][-1] for i in range(n)]


def analyse(A, B, C, s, exhaustive=True):
    d, m, p = len(A), len(B[0]), len(C)
    assert d > 0 and m > 0 and p > 0 and 1 <= s <= m
    powers = [identity(d)]
    for _ in range(max(2*d, d*(p+1))):
        powers.append(mm(A, powers[-1]))
    blocks = [columns(mm(mm(C, power), B)) for power in powers]
    tail = sum(blocks[d:2*d], [])
    r = rank(tail, p)
    H = d * (r + 1)
    cutvals = []
    for mask in range(1 << d):
        J = tuple(j for j in range(d) if mask >> j & 1)
        val = rank(tail + sum([blocks[j] for j in J], []), p) + s * (d-len(J))
        cutvals.append((val, J))
    mincut, Jmin = min(cutvals)
    supports = list(combinations(range(m), s))
    # Exhaustively enumerate every maximal allowed early support tuple.
    best, best_schedule = -1, None
    for sched in product(supports, repeat=d):
        early = [blocks[j][a] for j, S in enumerate(sched) for a in S]
        val = rank(tail + early, p)
        if val > best:
            best, best_schedule = val, sched
    assert mincut == best, (d, m, p, s, mincut, best)
    assert rank(tail, p) == r
    # Every later d-block spans the same output tail; all later columns lie in it.
    for k in range(d, H+1):
        later = sum(blocks[k:k+d], [])
        if k+d <= len(blocks):
            assert rank(later, p) == r
            assert rank(tail + later, p) == r
    # Directly compare the uncompressed intersection cuts and common independent sets.
    matroid_subsets = 0
    if d*m <= 8:
        ground = [(j, a) for j in range(d) for a in range(m)]
        uncompressed, common_max = d*m+p, 0
        for mask in range(1 << (d*m)):
            X = [e for i, e in enumerate(ground) if mask >> i & 1]
            lr = rank(tail + [blocks[j][a] for j, a in X], p) - r
            counts = [sum(j == k for j, a in X) for k in range(d)]
            f = lr + sum(min(s, m-k) for k in counts)
            uncompressed = min(uncompressed, f)
            if all(k <= s for k in counts) and lr == len(X):
                common_max = max(common_max, len(X))
            matroid_subsets += 1
        assert uncompressed == common_max == mincut-r
    schedule = []
    if mincut >= p:
        # Greedy independent columns modulo U from the exhaustive early optimum.
        early = []
        for j, S in enumerate(best_schedule):
            for a in S:
                v = blocks[j][a]
                if rank(tail+early+[v], p) > rank(tail+early, p):
                    early.append(v)
                    schedule.append((j, a))
        chosen_tail = []
        for t in range(r):
            found = False
            for j in range(d*(t+1), d*(t+2)):
                for a in range(m):
                    v = blocks[j][a]
                    if rank(chosen_tail+[v], p) > len(chosen_tail):
                        chosen_tail.append(v)
                        schedule.append((j, a))
                        found = True
                        break
                if found:
                    break
            assert found
        assert len({j for j, a in schedule if j >= d}) == r
        assert all(sum(j == k for j, a in schedule) <= s for k in range(H))
        selected = [blocks[j][a] for j, a in schedule]
        assert len(selected) == rank(selected, p) == p
        # Solve one arbitrary nonzero-state, nonzero-output instance, then run the
        # physical recurrence; this checks the time reversal and input coordinates.
        x0 = tuple(Q(i+1, i+2) for i in range(d))
        target = tuple(Q(2*i-1, i+1) for i in range(p))
        drift = mv(C, mv(powers[H], x0))
        coeffs = solve(selected, [x-y for x, y in zip(target, drift)])
        inputs = [[Q(0) for _ in range(m)] for _ in range(H)]
        for (j, a), z in zip(schedule, coeffs):
            inputs[H-1-j][a] = z
        x = x0
        for u in inputs:
            assert sum(z != 0 for z in u) <= s
            x = tuple(a+b for a, b in zip(mv(A, x), mv(B, u)))
        assert mv(C, x) == target
    direct_schedules = 0
    direct_max = None
    if exhaustive and comb(m, s)**H <= 4096:
        # This computation uses the actual whole-horizon columns, no quotient/cut.
        direct_max = 0
        for sched in product(supports, repeat=H):
            val = rank([blocks[j][a] for j, S in enumerate(sched) for a in S], p)
            direct_max = max(direct_max, val)
            direct_schedules += 1
        assert (direct_max == p) == (mincut >= p)
    return dict(d=d, m=m, p=p, s=s, tail_rank=r, horizon=H,
                mincut=mincut, minimizer=list(Jmin), feasible=mincut >= p,
                constructed_schedule=schedule, direct_max=direct_max,
                direct_schedules=direct_schedules, matroid_subsets=matroid_subsets)


def source_examples():
    A1 = mat([[0,1,0,0,0],[0,0,1,0,0],[0,0,0,0,0],
              [0,0,0,0,1],[0,0,0,0,0]])
    B1 = mat([[1,1],[0,0],[1,0],[0,0],[0,1]])
    C1 = mat([[1,0,0,0,0],[0,1,0,0,0],[0,0,0,1,0]])
    A2 = mat([[0,1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]])
    B2 = mat([[1,1],[1,0],[0,0],[0,1]])
    C2 = mat([[1,0,0,0],[0,0,1,0]])
    out = []
    for number, A, B, C, expected in [(1,A1,B1,C1,False),(2,A2,B2,C2,True)]:
        d, p = len(A), len(C)
        blocks, AB = [], B
        for _ in range(2*d+1):
            blocks.append(columns(mm(C, AB)))
            AB = mm(A, AB)
        ranks = [rank(sum(blocks[i:i+d], []), p) for i in range(d+1)]
        drops = [ranks[i]-ranks[i+1] for i in range(d)]
        avg = max(Q(sum(drops[:i+1]), i+1) for i in range(d))
        assert ranks[0] == p and max(drops) == 2 and avg == 1
        assert rank_by_minors(sum(blocks[:d], []), p) == ranks[0]
        res = analyse(A, B, C, 1)
        assert res['feasible'] == expected and res['mincut'] == 2
        res.update(example=number, ordinary_output_rank=ranks[0],
                   rank_drops=drops, max_drop=2, max_prefix_average=str(avg))
        if number == 2:
            only_second = [blocks[0][1], blocks[1][1]]
            assert rank(only_second, p) == rank_by_minors(only_second, p) == 2
        out.append(res)
    return out


def main():
    print('Arithmetic: fractions.Fraction; no floating point; seed=377')
    for x in source_examples():
        print('SOURCE_EXAMPLE ' + json.dumps(x, sort_keys=True))
    rng = Random(377)
    rows = []
    for d in range(1, 5):
        for trial in range(30):
            m, p = rng.randint(1, 3), rng.randint(1, min(4, d+1))
            A = [[rng.randint(-1, 1) for _ in range(d)] for _ in range(d)]
            if trial % 3 == 0:
                A = [[A[i][j] if i < j else 0 for j in range(d)] for i in range(d)]
            elif trial % 3 == 1:
                A[-1] = [0]*d
            B = [[rng.randint(-1, 1) for _ in range(m)] for _ in range(d)]
            C = [[rng.randint(-1, 1) for _ in range(d)] for _ in range(p)]
            s = rng.randint(1, m)
            rows.append(analyse(mat(A), mat(B), mat(C), s))
    print('RANDOM_SUMMARY ' + json.dumps(dict(
        systems=len(rows), feasible=sum(x['feasible'] for x in rows),
        infeasible=sum(not x['feasible'] for x in rows),
        full_horizon_exhaustive_systems=sum(x['direct_max'] is not None for x in rows),
        full_horizon_schedules=sum(x['direct_schedules'] for x in rows),
        all_ground_subsets_systems=sum(x['matroid_subsets'] > 0 for x in rows),
        all_ground_subsets=sum(x['matroid_subsets'] for x in rows)), sort_keys=True))
    # A non-prefix timing obstruction: k inputs matter only at delay one.
    for k in range(2, 5):
        d = 2*k
        A = [[0]*d for _ in range(d)]
        B = [[0]*k for _ in range(d)]
        C = [[0]*d for _ in range(k)]
        for i in range(k):
            A[k+i][i] = B[i][i] = C[i][k+i] = 1
        # Only delays 0 and 1 can contribute. A^2 = 0, so an exhaustive
        # two-delay check is also an exhaustive all-horizon check for this family.
        useful = [columns(mm(mat(C), mat(B))), columns(mm(mm(mat(C), mat(A)), mat(B)))]
        assert all(not z for row in mm(mat(A), mat(A)) for z in row)
        thresholds = []
        for s in range(1, k+1):
            supports = list(combinations(range(k), s))
            best = max(rank([useful[j][a] for j, S in enumerate(Ss) for a in S], k)
                       for Ss in product(supports, repeat=2))
            assert best == s
            thresholds.append(best == k)
        assert thresholds == [False]*(k-1)+[True]
        print('OBSTRUCTION_FAMILY ' + json.dumps(dict(k=k, minimum_sparsity=k,
              original_necessary_integer_bound=(k+1)//2, exact_max_ranks=list(range(1,k+1)))))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

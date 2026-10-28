#!/usr/bin/env python3
"""Exact finite checks for Csernák--Soukup, arXiv:2610.09650, Hypothesis 3.2(1).

Only Python's standard library and integer arithmetic are used. Recursive search
is cross-checked by a separately implemented bottom-up enumeration of all states.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import product
from math import comb

EDGES = ((0, 1), (0, 4), (0, 5), (0, 6), (1, 2), (1, 4),
         (2, 4), (2, 5), (3, 4), (4, 6), (5, 6))
INITIAL = (0, 0, 0, 9, 0, 3, 0)


def compositions(total, order):
    if order == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in compositions(total - first, order - 1):
                yield (first,) + tail


def stacked(c):
    return sum(x > 0 for x in c) == 1


def almost(c):
    return any(all(c[u] <= 1 for u in range(len(c)) if u != h)
               for h in range(len(c)))


def almost_at(total, order):
    configs = set()
    for h in range(order):
        for bits in product((0, 1), repeat=order - 1):
            pile = total - sum(bits)
            if pile >= 0:
                configs.add(bits[:h] + (pile,) + bits[h:])
    return configs


def recursive_checker(edges):
    directed = edges + tuple((v, u) for u, v in edges)

    def children(c):
        for u, v in directed:
            if c[u] >= 2:
                d = list(c)
                d[u] -= 2
                d[v] += 1
                assert sum(d) + 1 == sum(c)
                yield tuple(d)

    @lru_cache(None)
    def can_stack(c):
        return stacked(c) or any(can_stack(d) for d in children(c))

    return children, can_stack


def bottom_up(edges, order, max_total):
    # The transitions here use adjacency lists and coordinate expressions,
    # independently of the recursive checker's list mutation and directed list.
    neighbors = [set() for _ in range(order)]
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    good = set()
    bad_counts = []
    for total in range(max_total + 1):
        count = 0
        for c in compositions(total, order):
            ok = stacked(c)
            if not ok:
                for u in range(order):
                    if c[u] < 2:
                        continue
                    for v in neighbors[u]:
                        d = tuple(c[i] - (2 if i == u else 0)
                                  + (1 if i == v else 0) for i in range(order))
                        if d in good:
                            ok = True
                            break
                    if ok:
                        break
            if ok:
                good.add(c)
            else:
                count += 1
        bad_counts.append(count)
    return good, bad_counts


def fidelity():
    # Theorem 4.10 at n=2: every almost-stacked size-7 configuration on C_4 stacks.
    cycle4 = ((0, 1), (1, 2), (2, 3), (3, 0))
    _, check4 = recursive_checker(cycle4)
    a4 = almost_at(7, 4)
    assert len(a4) == 32
    assert a4 == {c for c in compositions(7, 4) if almost(c)}
    assert all(check4(c) for c in a4)
    print('FIDELITY Theorem 4.10: C4, total=7, almost_stacked=32, successes=32')
    # Theorem 4.4: stack(P_3)=2^3-1=7; enumerate each smaller total as well.
    path3 = ((0, 1), (1, 2))
    _, check3 = recursive_checker(path3)
    counts = [sum(not check3(c) for c in compositions(t, 3)) for t in range(2, 8)]
    assert all(n > 0 for n in counts[:-1]) and counts[-1] == 0
    print('FIDELITY Theorem 4.4: P3 bad_counts_at_2_to_7=' + str(counts))


def main():
    fidelity()
    assert len(set(EDGES)) == 11 and all(0 <= u < v < 7 for u, v in EDGES)
    seen = {0}
    while True:
        grown = seen | {v for u, v in EDGES if u in seen} | {u for u, v in EDGES if v in seen}
        if grown == seen:
            break
        seen = grown
    assert seen == set(range(7))
    assert {(0, 1), (0, 4), (1, 4)} <= set(EDGES)
    print('GRAPH order=7 edges=11 connected=True triangle=(0,1,4)')
    children, can_stack = recursive_checker(EDGES)
    good, unrestricted_bad = bottom_up(EDGES, 7, 15)
    failures = []
    for total in range(2, 16):
        configs = almost_at(total, 7)
        all_configs = list(compositions(total, 7))
        assert len(all_configs) == comb(total + 6, 6)
        assert configs == {c for c in all_configs if almost(c)}
        assert all(can_stack(c) == (c in good) for c in configs)
        failures.append(sum(not can_stack(c) for c in configs))
    assert failures[6:] == [74, 37, 17, 6, 0, 1, 1, 0]
    a12 = almost_at(12, 7)
    assert len(a12) == 448
    assert all(sum(c) == 12 and sum(x > 1 for x in c) == 1 for c in a12)
    assert all(can_stack(c) for c in a12)
    # Independently check concrete positive paths, with at most 11 moves each.
    lengths = Counter()
    for start in sorted(a12):
        c = start
        moves = 0
        while not stacked(c):
            d = next(d for d in children(c) if can_stack(d))
            delta = tuple(y - x for x, y in zip(c, d))
            sources = [i for i, x in enumerate(delta) if x == -2]
            targets = [i for i, x in enumerate(delta) if x == 1]
            assert len(sources) == len(targets) == 1
            u, v = sources[0], targets[0]
            assert (min(u, v), max(u, v)) in EDGES and c[u] >= 2
            assert all(delta[i] == 0 for i in range(7) if i not in (u, v))
            c = d
            moves += 1
            assert moves <= 11
        lengths[moves] += 1
    print('POSITIVE total=12 cases=448 checked_paths=448 lengths=' + str(dict(sorted(lengths.items()))))
    reachable = {INITIAL}
    queue = deque([INITIAL])
    while queue:
        c = queue.popleft()
        for d in children(c):
            if d not in reachable:
                reachable.add(d)
                queue.append(d)
    histogram = Counter(map(sum, reachable))
    assert len(reachable) == 225
    assert sum(1 for c in reachable for _ in children(c)) == 615
    assert all(d in reachable for c in reachable for d in children(c))
    assert all(sum(x > 0 for x in c) >= 2 for c in reachable)
    assert histogram == Counter({2: 13, 3: 23, 4: 35, 5: 54, 6: 36, 7: 25,
                                 8: 21, 9: 9, 10: 4, 11: 4, 12: 1})
    assert INITIAL not in good and not can_stack(INITIAL)
    print('NEGATIVE initial=' + str(INITIAL) + ' reachable=225 transition_occurrences=615 closed=True stacked_endpoints=0')
    print('NEGATIVE total_histogram=' + str(dict(sorted(histogram.items()))))
    for c in [(0, 0, 0, 11, 1, 1, 0), (0, 0, 0, 13, 0, 1, 0)]:
        assert almost(c) and not can_stack(c) and c not in good
    print('ALMOST_STACKED bad_counts_at_2_to_15=' + str(failures))
    print('UNRESTRICTED bad_counts_at_2_to_15=' + str(unrestricted_bad[2:]))
    assert all(n > 0 for n in unrestricted_bad[2:15]) and unrestricted_bad[15] == 0
    print('THRESHOLDS restricted_least=12 unrestricted_least=15')
    print('INDEPENDENT bottom_up_states=' + str(comb(22, 7)) + ' max_total=15 exact_integer_arithmetic=True')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

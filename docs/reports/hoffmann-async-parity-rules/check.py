#!/usr/bin/env python3
"""Exact checks for Hoffmann arXiv:2610.14026v1, Sections 2, 3.3, 3.6.

All arithmetic used as evidence is integer arithmetic. Graph reachability and
constructed paths are independent checks; neither is an all-size proof.
"""
from collections import Counter, deque
from itertools import product

BUNDLES = {
    'PM': frozenset((6, 9)), 'plmr': frozenset((1, 14)),
    'prml': frozenset((7, 8)), 'JLJR': frozenset((5, 10)),
    'R0R1': frozenset((4, 11)), 'L0L1': frozenset((2, 13)),
}
HC = BUNDLES['PM'] | BUNDLES['prml']
MIRROR_HC = BUNDLES['PM'] | BUNDLES['plmr']


def configurations(n):
    return product((0, 1), repeat=n)


def code(x, i):
    n = len(x)
    return sum(x[(i + j - 1) % n] << (3 - j) for j in range(4))


def step(x, i, active):
    if code(x, i) not in active:
        return x
    y = list(x)
    # Simultaneous complement of the two selected cells.
    for j in {i, (i + 1) % len(x)}:
        y[j] = 1 - x[j]
    return tuple(y)


def run_count(x):
    return sum(x[i] != x[(i + 1) % len(x)] for i in range(len(x)))


def rules():
    result = []
    for left, right in ((True, False), (False, True), (True, True)):
        for optional in product((False, True), repeat=3):
            selected = ['PM']
            selected += ['plmr'] if left else []
            selected += ['prml'] if right else []
            selected += [name for name, take in zip(('JLJR', 'R0R1', 'L0L1'), optional) if take]
            result.append((tuple(selected), frozenset().union(*(BUNDLES[b] for b in selected))))
    assert len(result) == len(set(a for _, a in result)) == 24
    return result


def orbit_representative(x):
    return min(int(''.join(map(str, x[i:] + x[:i])), 2) for i in range(len(x)))


def fidelity():
    # Appendix 2, Rule HC, matrix (I), n=5. Rows/columns are shift
    # representatives 1,7,11,31; entries are selection counts, divided by 5.
    reps = (1, 7, 11, 31)
    expected = ((4, 1, 0, 0), (1, 3, 0, 1), (1, 0, 4, 0), (0, 0, 0, 5))
    got = []
    for r in reps:
        x = tuple(map(int, format(r, '05b')))
        counts = Counter(orbit_representative(step(x, i, HC)) for i in range(5))
        assert set(counts) <= set(reps)
        got.append(tuple(counts[s] for s in reps))
    assert tuple(got) == expected, (got, expected)
    print('FIDELITY HC Appendix 2 n=5 transition-count matrix:', got)
    # Section 3.1 uses b for the number of runs of ONE color, not both.
    for word, expected_b in [('11111', 0), ('00100', 1), ('100111011', 2),
                             ('1010101', 3), ('0101010', 3)]:
        assert run_count(tuple(map(int, word))) == 2 * expected_b
    print('FIDELITY Section 3.1: all five block-count examples agree')
    # Section 3.6 HC's Boolean formula agrees with Table 1's active indices.
    for bits in product((0, 1), repeat=4):
        index = int(''.join(map(str, bits)), 2)
        l, z, r, _ = bits
        assert (index in HC) == (l != z and z == r)
    print('FIDELITY HC active indices:', sorted(HC))


def graph_check(n, active):
    states = list(configurations(n))
    reverse = {x: [] for x in states}
    zeros, ones = (0,) * n, (1,) * n
    transitions = 0
    for x in states:
        for i in range(n):
            y = step(x, i, active)
            assert sum(x) % 2 == sum(y) % 2
            if x in (zeros, ones):
                assert y == x
            reverse[y].append(x)
            transitions += 1
    reached = {zeros, ones}
    queue = deque(reached)
    while queue:
        for x in reverse[queue.popleft()]:
            if x not in reached:
                reached.add(x)
                queue.append(x)
    assert len(reached) == len(states), (n, active, len(reached))
    return len(states), transitions


def eliminate(x):
    """Explicit HC run erasure; returns the prescribed selected-site sequence."""
    n = len(x)
    initial_parity = sum(x) % 2
    sites = []
    while run_count(x):
        before = run_count(x)
        assert before > 0 and before % 2 == 0
        even_runs = []
        for start in range(n):
            if x[(start - 1) % n] == x[start]:
                continue
            length = 1
            while x[(start + length) % n] == x[start]:
                length += 1
            if length % 2 == 0:
                even_runs.append((start, length))
        assert even_runs
        start, length = even_runs[0]
        for j in range(length // 2):
            i = (start + 2 * j) % n
            assert x[(i - 1) % n] != x[i] == x[(i + 1) % n]
            assert code(x, i) in HC
            y = list(x)
            y[i] = x[(i - 1) % n]
            y[(i + 1) % n] = x[(i - 1) % n]
            y = tuple(y)
            assert step(x, i, HC) == y
            x = y
            sites.append(i)
        assert run_count(x) == before - 2
    assert x == (initial_parity,) * n
    assert len(sites) <= n * n
    return sites


def main():
    fidelity()
    family = rules()
    total_states = total_transitions = 0
    for n in (1, 3, 5, 7, 9, 11):
        for _, active in family:
            assert HC <= active or MIRROR_HC <= active
            states, transitions = graph_check(n, active)
            total_states += states
            total_transitions += transitions
        print(f'GRAPH n={n}: all 24 rules preserve parity, absorb constants, and reach a constant from every state')
    print(f'GRAPH_TOTAL rule/configuration pairs={total_states}; selections={total_transitions}')
    path_states = 0
    for n in range(1, 16, 2):
        maximum = 0
        for x in configurations(n):
            path = eliminate(x)
            maximum = max(maximum, len(path))
            # Reflection translates selected pair i,i+1 to n-2-i,n-1-i.
            reflected = tuple(reversed(x))
            for i in path:
                reflected = step(reflected, (n - 2 - i) % n, MIRROR_HC)
            assert reflected == (sum(x) % 2,) * n
            path_states += 1
        print(f'PATH n={n}: HC and reflected HC verified; maximum prescribed length={maximum}; bound={n*n}')
    print(f'PATH_TOTAL initial words={path_states}; both orientations checked')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

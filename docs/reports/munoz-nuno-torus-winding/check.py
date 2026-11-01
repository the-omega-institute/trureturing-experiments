#!/usr/bin/env python3
"""Exact checks for the threshold-three torus-winding question, arXiv:2610.37990v1.

Updates enumerate physical neighboring sites, without a residue-count shortcut.
The output contains finite evidence; the infinite family requires the argument
in the preregistration and subsequent formalization.
"""
import json
from collections import Counter
from fractions import Fraction

OFFSETS = tuple((a, b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))
F = (2, 1, 0, 2)
G = (2, 2, 0, 2)


def neighbors(size, i, j):
    assert size >= 3
    sites = {((i + a) % size, (j + b) % size) for a, b in OFFSETS}
    assert len(sites) == 8 and (i, j) not in sites
    return sites


def counts(state):
    size = len(state)
    assert all(len(row) == size for row in state)
    assert all(c in (0, 1, 2) for row in state for c in row)
    return [[sum(state[u][v] == (state[i][j] + 1) % 3
                 for u, v in neighbors(size, i, j))
             for j in range(size)] for i in range(size)]


def step(state, threshold=3):
    pc = counts(state)
    return [[(c + 1) % 3 if pc[i][j] >= threshold else c
             for j, c in enumerate(row)] for i, row in enumerate(state)]


def increment(a, b):
    return (b - a + 1) % 3 - 1


def loop_winding(loop):
    total = sum(increment(loop[k], loop[(k + 1) % len(loop)]) for k in range(len(loop)))
    assert total % 3 == 0
    return total // 3


def windings(state):
    size = len(state)
    return {'horizontal': [loop_winding(row) for row in state],
            'vertical': [loop_winding([state[i][j] for i in range(size)]) for j in range(size)]}


def square_sums(state):
    size = len(state)
    result = []
    for i in range(size):
        row = []
        for j in range(size):
            loop = [state[i][j], state[i][(j + 1) % size],
                    state[(i + 1) % size][(j + 1) % size], state[(i + 1) % size][j]]
            row.append(sum(increment(loop[k], loop[(k + 1) % 4]) for k in range(4)))
        result.append(row)
    return result


def population(state):
    counter = Counter(c for row in state for c in row)
    return [counter[c] for c in range(3)]


def diagonal(size, profile):
    assert size % 4 == 0
    return [[profile[(j - i) % 4] for j in range(size)] for i in range(size)]


def fidelity():
    # Section II: alpha=1/3 and alpha=3/8 give threshold three, exactly.
    for n in range(9):
        assert (Fraction(n, 8) >= Fraction(1, 3)) == (n >= 3)
        assert (Fraction(n, 8) >= Fraction(3, 8)) == (n >= 3)
    # Section V.4: axis-aligned bands translate by one site; exactly three
    # front columns advance, with exactly three predator neighbors each.
    records = []
    for size in range(3, 13):
        tested = 0
        for w0 in range(1, size - 1):
            for w1 in range(1, size - w0):
                w2 = size - w0 - w1
                row = [0] * w0 + [1] * w1 + [2] * w2
                state = [row[:] for _ in range(size)]
                expected = [row[1:] + row[:1] for _ in range(size)]
                pc = counts(state)
                actual = step(state)
                assert actual == expected
                assert sum(actual[i][j] != state[i][j]
                           for i in range(size) for j in range(size)) == 3 * size
                assert Counter(n for line in pc for n in line) == Counter({3: 3 * size, 0: size * (size - 3)})
                assert population(actual) == [size * w for w in (w0, w1, w2)]
                assert all(w == 1 for w in windings(state)['horizontal'])
                assert all(w == 0 for w in windings(state)['vertical'])
                assert all(x == 0 for line in square_sums(state) for x in line)
                tested += 1
        records.append({'L': size, 'positive_width_triples_checked': tested,
                        'changed_fraction': str(Fraction(3, size))})
    state = [[0, 0, 1, 1, 2, 2] for _ in range(6)]
    initial = [row[:] for row in state]
    for t in range(1, 7):
        state = step(state)
        assert (state == initial) == (t == 6)
    return {'reference_threshold_equivalence': True, 'axis_band_translation': records,
            'equal_width_L6_exact_state_period': 6}


def family_checks():
    records, witness = [], None
    for m in range(1, 17):
        size = 4 * m
        before, target = diagonal(size, F), diagonal(size, G)
        wb, wt = windings(before), windings(target)
        assert population(before) == [size * size // 4, size * size // 4, size * size // 2]
        assert wb == {'horizontal': [-m] * size, 'vertical': [m] * size}
        assert wt == {'horizontal': [0] * size, 'vertical': [0] * size}
        for state in (before, target):
            assert all(x == 0 for line in square_sums(state) for x in line)
        pc = counts(before)
        for i in range(size):
            for j in range(size):
                assert pc[i][j] == (4 if before[i][j] == 1 else 2)
                local = Counter(before[u][v] for u, v in neighbors(size, i, j))
                assert local == Counter({0: 2, 1: 2, 2: 4})
        for threshold in (3, 4):
            assert step(before, threshold) == target
            assert step(target, threshold) == target
        records.append({'m': m, 'L': size, 'thresholds': [3, 4],
                        'initial_population': population(before),
                        'initial_winding_pair': [-m, m], 'final_winding_pair': [0, 0],
                        'zero_square_sums_before_and_after': size * size,
                        'one_step_target_and_fixed_point': True})
        if size == 8:
            witness = {'L': size, 'initial_state': before, 'initial_predator_counts': pc,
                       'one_step_fixed_state': target, 'fixed_state_predator_counts': counts(target),
                       'initial_windings': wb, 'final_windings': wt,
                       'initial_square_sums': square_sums(before),
                       'final_square_sums': square_sums(target)}
    # Every possible four-residue color profile has zero elementary-square sums.
    profile_count = 0
    for code in range(3 ** 4):
        profile = tuple((code // (3 ** k)) % 3 for k in range(4))
        state = diagonal(4, profile)
        assert all(x == 0 for line in square_sums(state) for x in line)
        profile_count += 1
    return {'family_instances': records, 'profiles_with_zero_square_sums_checked': profile_count,
            'exact_L8_witness': witness}


def main():
    result = {'source': 'arXiv:2610.37990v1, Sections II and V.4',
              'arithmetic': 'integers and fractions.Fraction; no floating point',
              'fidelity': fidelity(), 'expected_refutation': family_checks()}
    print(json.dumps(result, indent=2, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

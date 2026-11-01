#!/usr/bin/env python3
"""Exact checks for arXiv:2610.30001v1, Conjecture 1 item 2.

Source fidelity: Section IV-C's degeneracy example and Table I, k=5,t=2.
Candidate: a 71-row exact-4 witness, padding, and translation by binary XOR.
Only Python integers, finite sets, and exact combinatorial counts are used.
This program is experimental evidence, not a Lean proof.
"""
from collections import Counter
from itertools import combinations
from functools import lru_cache
from math import comb


def word(s):
    return tuple(int(b) for b in s)


def distance(a, b):
    return sum(x != y for x, y in zip(a, b))


def column_counts(rows):
    return tuple(sum(row[j] for row in rows) for j in range(len(rows[0])))


def majority(rows):
    n, k = len(rows), len(rows[0])
    counts = column_counts(rows)
    return None if any(2 * c == n for c in counts) else tuple(int(2 * c > n) for c in counts)


def weighted_votes(rows, budgets, remaining):
    return tuple(sum((remaining - 2 * s) * (1 - 2 * row[j])
                     for row, s in zip(rows, budgets)) for j in range(len(rows[0])))


def decode(rows, t):
    """Algorithm 1: depth-first, increasing-column-order maximal-reliability ties.

    A zero vote, invalid signed residual budget, or failing final distance check
    rejects that child, restoring its parent's complete state before the next tie.
    No nonmaximal coordinate or alternative sign is explored. Success stops DFS.
    """
    k, n = len(rows[0]), len(rows)
    stats = Counter()
    first_steps = []

    def dfs(todo, assigned, budgets, remaining, history):
        stats['states'] += 1
        assert remaining == len(todo)
        assert set(todo).isdisjoint(assigned)
        assert len(todo) + len(assigned) == k
        assert budgets == tuple(t - sum(row[j] != b for j, b in assigned.items()) for row in rows)
        if not todo:
            estimate = tuple(assigned[j] for j in range(k))
            if all(distance(row, estimate) == t for row in rows):
                stats['accepted'] += 1
                return estimate
            stats['final_rejected'] += 1
            return None
        votes = weighted_votes(rows, budgets, remaining)
        highest = max(abs(votes[j]) for j in todo)
        tied = tuple(j for j in todo if abs(votes[j]) == highest)
        if len(first_steps) < 2:
            first_steps.append((tuple(todo), tuple(votes[j] for j in todo), tied))
        for j in tied:
            stats['choices'] += 1
            if votes[j] == 0:
                stats['zero_vote'] += 1
                continue
            bit = int(votes[j] < 0)
            next_budgets = tuple(s - int(row[j] != bit) for row, s in zip(rows, budgets))
            if any(s < 0 or s > remaining - 1 for s in next_budgets):
                stats['invalid_budget'] += 1
                continue
            result = dfs(tuple(l for l in todo if l != j), assigned | {j: bit},
                         next_budgets, remaining - 1, history + ((j, bit),))
            if result is not None:
                return result
        return None

    result = dfs(tuple(range(k)), {}, (t,) * n, k, ())
    return result, dict(sorted(stats.items())), first_steps


def independent_accepted_words(rows, t):
    """Independent state-DAG enumeration, recomputing all budgets from assignments.

    Unlike decode(), this visits invalid states, merges repeated assignments,
    explores every branch even after an acceptance, and stores accepted words.
    """
    stats = Counter()
    k = len(rows[0])

    @lru_cache(None)
    def visit(assignment):
        stats['states'] += 1
        todo = [j for j, b in enumerate(assignment) if b == -1]
        r = len(todo)
        budgets = [t - sum(b != -1 and row[j] != b for j, b in enumerate(assignment)) for row in rows]
        if any(s < 0 or r < s for s in budgets):
            stats['invalid_budget'] += 1
            return frozenset()
        if r == 0:
            return frozenset((assignment,)) if all(distance(row, assignment) == t for row in rows) else frozenset()
        votes = {j: sum((r - 2 * s) * (-1 if row[j] else 1) for row, s in zip(rows, budgets)) for j in todo}
        largest = max(map(abs, votes.values()))
        accepted = set()
        for j in todo:
            if abs(votes[j]) != largest or votes[j] == 0:
                continue
            child = list(assignment)
            child[j] = int(votes[j] < 0)
            accepted.update(visit(tuple(child)))
        return frozenset(accepted)

    accepted = visit((-1,) * k)
    return accepted, dict(sorted(stats.items()))


def validate_channel(rows, t, u):
    assert rows and len(set(rows)) == len(rows)
    assert all(len(row) == len(u) and all(b in (0, 1) for b in row) for row in rows)
    assert all(distance(row, u) == t for row in rows)
    assert 1 <= t and 2 * t < len(u)
    assert 3 <= len(rows) <= comb(len(u), t) and len(rows) % 2 == 1


def fidelity():
    rows = tuple(map(word, ('11000', '10100', '10010', '10001', '01100', '01010', '00110')))
    validate_channel(rows, 2, (0,) * 5)
    assert column_counts(rows) == (4, 3, 3, 3, 1)
    assert weighted_votes(rows, (2,) * 7, 5) == (-1, 1, 1, 1, 5)
    budgets = tuple(2 - row[4] for row in rows)
    assert budgets == (2, 2, 2, 1, 2, 2, 2)
    assert weighted_votes(rows, budgets, 4) == (-2, 2, 2, 2, -2)
    answer, stats, _ = decode(rows, 2)
    assert answer == (0,) * 5 and stats.get('invalid_budget', 0) > 0
    print('FIDELITY degeneracy example: counts=(4,3,3,3,1); first=4->0; tied second votes=(-2,2,2,2); backtracking returns 00000')
    sphere = tuple(tuple(int(j in support) for j in range(5)) for support in combinations(range(5), 2))
    for n, expected_total, expected_majority_failures in [(7, 120, 90), (9, 10, 0)]:
        total = weighted_failures = majority_failures = 0
        for subset in combinations(sphere, n):
            total += 1
            weighted_failures += decode(subset, 2)[0] != (0,) * 5
            majority_failures += majority(subset) != (0,) * 5
        assert (total, weighted_failures, majority_failures) == (expected_total, 0, expected_majority_failures)
        print(f'FIDELITY Table I k=5 t=2 N={n}: total={total}; weighted_failures={weighted_failures}; majority_failures={majority_failures}')


def construct():
    v = tuple(range(2, 9))
    triples = tuple(frozenset(v[(i + d) % 7] for d in range(3)) for i in range(7))
    a = tuple(frozenset((0, 1, *pair)) for pair in combinations(v, 2)) + tuple(s | {0} for s in triples)
    b = tuple(frozenset(s) for s in combinations(v, 4)) + tuple(s | {1} for s in triples) + (frozenset((1, 2, 3, 5)),)
    assert len(a) == 28 and len(b) == 43 and len(set(a + b)) == 71
    rows = tuple(tuple(int(j in s) for j in range(9)) for s in a + b)
    supplied = '''111100000 111010000 111001000 111000100 111000010 111000001 110110000 110101000 110100100 110100010 110100001 110011000 110010100 110010010 110010001 110001100 110001010 110001001 110000110 110000101 110000011 101110000 100111000 100011100 100001110 100000111 101000011 101100001 001111000 001110100 001110010 001110001 001101100 001101010 001101001 001100110 001100101 001100011 001011100 001011010 001011001 001010110 001010101 001010011 001001110 001001101 001001011 001000111 000111100 000111010 000111001 000110110 000110101 000110011 000101110 000101101 000101011 000100111 000011110 000011101 000011011 000010111 000001111 011110000 010111000 010011100 010001110 010000111 011000011 011100001 011101000'''
    assert rows == tuple(map(word, supplied.split()))
    return rows


def counterexample(rows):
    validate_channel(rows, 4, (0,) * 9)
    assert column_counts(rows[:28]) == (28, 21, 9, 9, 9, 9, 9, 9, 9)
    counts = column_counts(rows)
    assert counts == (28, 29, 33, 33, 32, 33, 32, 32, 32)
    assert all(2 * c < 71 for c in counts) and majority(rows) == (0,) * 9
    votes = weighted_votes(rows, (4,) * 71, 9)
    assert votes == (15, 13, 5, 5, 7, 5, 7, 7, 7)
    budgets = tuple(4 - row[0] for row in rows)
    weights = tuple(8 - 2 * s for s in budgets)
    assert weights == (2,) * 28 + (0,) * 43
    next_votes = weighted_votes(rows, budgets, 8)[1:]
    assert next_votes == (-28,) + (20,) * 7
    answer, stats, steps = decode(rows, 4)
    assert answer is None and steps[0][2] == (0,) and steps[1][2] == (1,)
    centers = tuple(mask for mask in range(1 << 9) if all(distance(row, tuple((mask >> j) & 1 for j in range(9))) == 4 for row in rows))
    assert centers == (0,)
    assert 2 * comb(7, 3) + 1 == 71
    nwmaj = (comb(8, 3) * 9) // (9 - 4) + 1
    assert nwmaj == 101 and 2 * comb(8, 3) + 1 == 113
    print(f'BASE k=9 t=4 N=71: column_errors={counts}')
    print(f'BASE first_votes={votes}; second_votes={next_votes}; A_weights=2 (28 rows); B_weights=0 (43 rows)')
    accepted, dag_stats = independent_accepted_words(rows, 4)
    assert not accepted and dag_stats == {'invalid_budget': 8, 'states': 25}
    print(f'BASE forced decisions: 0->0 then 1->1; algorithm1=FAILURE; DFS stats={stats}')
    print(f'BASE independent assignment DAG, including invalid terminals: stats={dag_stats}; accepted_words=0')
    print('BASE consistent centers among all 512 words: only 000000000; N_wc=71; N_wmaj=101; N_maj=113')


def family(rows):
    # Exact finite corroboration of the general padding argument, not an unbounded proof.
    for m in range(6):
        padded = tuple(row + (0,) * m for row in rows)
        k = 9 + m
        validate_channel(padded, 4, (0,) * k)
        assert majority(padded) == (0,) * k
        for decoded_padding in range(m):
            r = k - decoded_padding
            common_weight = r - 8
            assert common_weight > 0
            assert 71 * common_weight > 15 * common_weight
        result, stats, _ = decode(padded, 4)
        assert result is None
        centers = [mask for mask in range(1 << k) if all(sum(row[j] != ((mask >> j) & 1) for j in range(k)) == 4 for row in padded)]
        assert centers == [0]
        print(f'PADDING k={k}: FAILURE; unique_center=0; states={stats["states"]}; invalid_budget={stats.get("invalid_budget", 0)}')
    # Check all translations of the base instance, plus representative padded translations.
    for m in (0, 1, 3):
        padded = tuple(row + (0,) * m for row in rows)
        k = 9 + m
        masks = range(1 << k) if m == 0 else (0, 1, (1 << k) - 1, sum(1 << j for j in range(0, k, 2)))
        checked = 0
        for mask in masks:
            u = tuple((mask >> j) & 1 for j in range(k))
            translated = tuple(tuple(b ^ u[j] for j, b in enumerate(row)) for row in padded)
            validate_channel(translated, 4, u)
            assert majority(translated) == u
            assert decode(translated, 4)[0] is None
            checked += 1
        print(f'XOR k={k}: {checked} transmitted words checked; majority correct; weighted decoder FAILURE')


def main():
    fidelity()
    rows = construct()
    counterexample(rows)
    family(rows)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

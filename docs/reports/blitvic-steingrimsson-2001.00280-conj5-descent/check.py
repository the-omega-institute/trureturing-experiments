#!/usr/bin/env python3
"""Exact finite checks of the descent-convention ambiguity in Conjecture 5.

No general settlement is claimed. Only Python integer arithmetic is used.
"""

from collections import Counter
from itertools import combinations, permutations, product
from math import comb


def arrangements(n, k):
    """Definition 6: permutation with a [k]-color on each fixed point."""
    for perm in permutations(range(1, n + 1)):
        fixed = [i for i, value in enumerate(perm, 1) if i == value]
        for colors in product(range(1, k + 1), repeat=len(fixed)):
            yield perm, dict(zip(fixed, colors))


def forms(n, k, derangement=False):
    """Definition 7, enumerated independently of colored permutations."""
    alphabet = tuple(range(-k if derangement else 1 - k, 0))
    for negative_count in range(n + 1):
        m = n - negative_count
        for positive in permutations(range(1, m + 1)):
            if derangement and any(i == value for i, value in enumerate(positive, 1)):
                continue
            for positions in combinations(range(n), negative_count):
                for negative in product(alphabet, repeat=negative_count):
                    assigned = dict(zip(positions, negative))
                    positive_iter = iter(positive)
                    yield tuple(assigned[i] if i in assigned else next(positive_iter)
                                for i in range(n))


def word_form(perm, colors, k, derangement=False):
    raw = [-colors[i] if i in colors and (derangement or colors[i] != k) else value
           for i, value in enumerate(perm, 1)]
    ranks = {value: rank for rank, value in enumerate(sorted(x for x in raw if x > 0), 1)}
    return tuple(ranks[x] if x > 0 else x for x in raw)


def descents(word):
    return sum(left > right for left, right in zip(word, word[1:]))


def colored_descents(perm, colors, k, convention, terminal):
    def color(i):
        if i not in colors:
            return 0
        if convention == 'literal':
            return colors[i]
        if convention == 'zero_based':
            return colors[i] - 1
        assert convention == 'cyclic' and k > 0
        return colors[i] % k

    # Definition 4 compares color first, then value.
    pairs = tuple((color(i), value) for i, value in enumerate(perm, 1))
    if terminal:
        pairs += ((0, len(perm) + 1),)
    return descents(pairs)


def coefficients(distribution, n):
    return tuple(distribution[i] for i in range(n + 1))


def main():
    # Fidelity: Proposition 2(1), not the conjectured statistic comparison.
    counts = {}
    for k in range(4):
        for n in range(7):
            count = sum(1 for _ in arrangements(n, k))
            assert count == sum(1 for _ in forms(n, k, derangement=True))
            assert count == (1 if n == 0 else n * counts[n - 1, k] + (k - 1) ** n)
            counts[n, k] = count
        print('FIDELITY recurrence k=', k, 'A_k(0..6)=', tuple(counts[n, k] for n in range(7)))

    # Fidelity: the source's explicit 4-arrangement representation example.
    example = (6, 2, 1, 4, 5, 7, 3)
    coloring = {2: 4, 4: 1, 5: 4}
    assert word_form(example, coloring, 4, True) == (3, -4, 1, -1, -4, 4, 2)
    assert word_form(example, coloring, 4) == (5, 2, 1, -1, 4, 6, 3)
    print('FIDELITY source 4-arrangement example OK')

    boxes = 0
    for n in range(7):
        for k in range(1, 4):
            objects = tuple(arrangements(n, k))
            permutation_words = tuple(forms(n, k))
            derangement_words = tuple(forms(n, k, derangement=True))
            assert len(set(permutation_words)) == len(objects) == counts[n, k]
            assert len(set(derangement_words)) == len(objects)
            assert {word_form(p, f, k) for p, f in objects} == set(permutation_words)
            assert {word_form(p, f, k, True) for p, f in objects} == set(derangement_words)
            p_distribution = Counter(map(descents, permutation_words))
            d_distribution = Counter(map(descents, derangement_words))
            assert p_distribution == d_distribution
            literal = Counter(colored_descents(p, f, k, 'literal', True) for p, f in objects)
            shifted = Counter(colored_descents(p, f, k, 'zero_based', True) for p, f in objects)
            cyclic = Counter(colored_descents(p, f, k, 'cyclic', True) for p, f in objects)
            assert shifted == cyclic
            assert p_distribution[0] == comb(n + k - 1, k - 1)
            if n:
                assert literal[0] == 0
                assert shifted[0] == 1
            for convention in ('literal', 'zero_based', 'cyclic'):
                internal = Counter(colored_descents(p, f, k, convention, False)
                                   for p, f in objects)
                assert internal == p_distribution, (n, k, convention, internal, p_distribution)
            if n:
                boxes += 1
                print('SMALL_CASE n=', n, 'k=', k,
                      'word=', coefficients(p_distribution, n),
                      'literal_terminal=', coefficients(literal, n),
                      'zero_based_terminal=', coefficients(shifted, n),
                      'all_three_internal_conventions_match=True')
            if (n, k) == (2, 2):
                assert coefficients(literal, n) == (0, 4, 1)
                assert coefficients(shifted, n) == (1, 4, 0)
                assert coefficients(p_distribution, n) == (3, 2, 0)
                for p, f in objects:
                    print('WITNESS', p, tuple(f.get(i, 0) for i in range(1, n + 1)),
                          'literal_terminal=', colored_descents(p, f, k, 'literal', True),
                          'literal_internal=', colored_descents(p, f, k, 'literal', False))
                print('INDEPENDENT_PERMUTATION_FORMS', sorted(permutation_words))

    assert boxes == 18
    print('CONCLUSION: literal mismatch reproduced; terminal-convention robustness fails.')
    print('SCOPE: 0<=n<=6, 1<=k<=3; finite agreement does not prove an unbounded identity.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

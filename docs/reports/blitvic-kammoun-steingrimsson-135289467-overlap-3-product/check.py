#!/usr/bin/env python3
"""Exact finite checks of arXiv:2610.29846v1, Section 7 item 10.

All counts start from literal prescribed-window comparisons.  The conjectured
recurrence and the proposed poset decomposition are never inputs to counting.
Only Python 3 standard-library integer arithmetic is used.
"""
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, prod

PATTERN = (1, 3, 5, 2, 8, 9, 4, 6, 7)


def window_poset(pattern, overlap, occurrences):
    """Strict predecessor masks from all pairs in each prescribed window."""
    k = len(pattern)
    assert sorted(pattern) == list(range(1, k + 1))
    assert 1 <= overlap < k and occurrences >= 1
    step = k - overlap
    n = k + (occurrences - 1) * step
    pred = [0] * n
    for t in range(occurrences):
        for a in range(k):
            for b in range(k):
                if pattern[a] < pattern[b]:
                    pred[t * step + b] |= 1 << (t * step + a)
    return tuple(pred)


def closure(pred):
    p = list(pred)
    for k in range(len(p)):
        for y in range(len(p)):
            if p[y] & (1 << k):
                p[y] |= p[k]
    assert all(not (p[x] & (1 << x)) for x in range(len(p))), 'cycle'
    return tuple(p)


def extension_count(pred):
    """Ideal DP: append each enabled position; independent of any product law."""
    n = len(pred)
    full = (1 << n) - 1

    @lru_cache(None)
    def visit(ideal):
        if ideal == full:
            return 1
        total = 0
        remaining = full ^ ideal
        while remaining:
            bit = remaining & -remaining
            remaining -= bit
            v = bit.bit_length() - 1
            if pred[v] & ideal == pred[v]:
                total += visit(ideal | bit)
        return total

    result = visit(0)
    return result, visit.cache_info().currsize


def restricted(pred, vertices):
    p = closure(pred)
    return tuple(sum(1 << i for i, x in enumerate(vertices) if p[y] & (1 << x))
                 for y in vertices)


def packs(word, pattern, overlap, occurrences):
    k = len(pattern)
    assert len(word) == k + (occurrences - 1) * (k - overlap)
    return all((word[t * (k - overlap) + a] < word[t * (k - overlap) + b])
               == (pattern[a] < pattern[b])
               for t in range(occurrences) for a in range(k) for b in range(k))


def labelled_positions(j):
    labels = {}
    for t in range(j + 1):
        labels['a', t] = 6 * t
        labels['c', t] = 6 * t + 1
        labels['e', t] = 6 * t + 2
    for t in range(j):
        labels['b', t] = 6 * t + 3
        labels['u', t] = 6 * t + 4
        labels['v', t] = 6 * t + 5
    assert set(labels.values()) == set(range(6 * j + 3))
    return labels


def chain_presentation(j):
    labels = labelled_positions(j)
    p = [0] * len(labels)
    for t in range(j):
        chain = [('a', t), ('b', t), ('c', t), ('a', t + 1),
                 ('e', t), ('c', t + 1), ('e', t + 1), ('u', t), ('v', t)]
        for x, y in zip(chain, chain[1:]):
            p[labels[y]] |= 1 << labels[x]
    return tuple(p)


HEAD = (('a', 0), ('b', 0), ('c', 0), ('a', 1), ('e', 0), ('b', 1), ('c', 1))
REMOVED = (('a', 2), ('b', 2), ('e', 1), ('c', 2), ('u', 0), ('v', 0))


def remainder_labels(j):
    labels = labelled_positions(j)
    return [key for key in labels if key not in HEAD]


def glue(pred):
    """Indices 0..5 = a,b,e,c,u,v, then the old poset R."""
    p = [0, 1 << 0, 1 << 0, (1 << 1) | (1 << 2), 1 << 2, 1 << 4]
    p.extend((bits << 6) | (1 << 3) for bits in pred)
    return tuple(p)


def check_structure(j):
    p = closure(window_poset(PATTERN, 3, j))
    assert p == closure(chain_presentation(j)), ('position chains', j)
    labels = labelled_positions(j)
    head = [labels[key] for key in HEAD]
    rest_keys = remainder_labels(j)
    rest = [labels[key] for key in rest_keys]
    assert len(rest) == 6 * j - 4
    assert extension_count(restricted(p, head))[0] == 2
    assert all(p[y] & (1 << x) for x in head for y in rest), ('head', j)
    q = restricted(p, rest)
    if j >= 3:
        previous_keys = remainder_labels(j - 1)
        # Align R by the explicit subscript-decrement bijection, not a graph search.
        survivors = [(kind, t + 1) for kind, t in previous_keys]
        assert set(survivors) == set(rest_keys) - set(REMOVED)
        previous_labels = labelled_positions(j - 1)
        previous_q = restricted(window_poset(PATTERN, 3, j - 1),
                                [previous_labels[key] for key in previous_keys])
        actual_q = restricted(p, [labels[key] for key in REMOVED + tuple(survivors)])
        assert actual_q == closure(glue(previous_q)), ('exact gluing', j)
    return q


def fidelity():
    # Source Remark 2.7: prescribed starts 1,5; an additional start 3 is allowed.
    example = (1, 9, 2, 8, 3, 7, 4, 6, 5)
    small_pattern = (1, 5, 2, 4, 3)
    assert packs(example, small_pattern, 1, 2)
    extra = tuple(example[2:7])
    assert all((extra[a] < extra[b]) == (small_pattern[a] < small_pattern[b])
               for a in range(5) for b in range(5))
    print('FIDELITY Remark 2.7: 192837465 packs 15243 at starts 1,5 and additionally 3')
    # Source Proposition 5.25 (prefix 12, suffix 34): product binom(4i+2,2).
    for j in range(1, 6):
        value, _ = extension_count(window_poset((1, 2, 5, 6, 3, 4), 2, j))
        expected = prod(comb(4 * i + 2, 2) for i in range(j))
        assert value == expected
        print(f'FIDELITY Proposition 5.25: pattern=125634 overlap=2 j={j} count={value}')
    # Independent direct enumeration checks the poset-to-permutation count at n=9.
    brute = sum(packs(w, PATTERN, 3, 1) for w in permutations(range(1, 10)))
    assert brute == extension_count(window_poset(PATTERN, 3, 1))[0] == 1
    print('DIRECT permutations: target j=1, all 9! actual permutations, count=1')


def gluing_checks():
    cases = 0
    for s in range(5):
        edges = list(combinations(range(s), 2))
        for mask in range(1 << len(edges)):
            p = [0] * s
            for bit, (x, y) in enumerate(edges):
                if mask & (1 << bit):
                    p[y] |= 1 << x
            old = extension_count(p)[0]
            new = extension_count(glue(p))[0]
            assert new == (s + 3) ** 2 * old
            cases += 1
    print(f'GLUING: {cases} naturally labelled DAG presentations on 0..4 vertices checked')


def main():
    fidelity()
    actual_counts = []
    for j in range(1, 10):
        p = window_poset(PATTERN, 3, j)
        count, states = extension_count(p)
        expected = 1 if j == 1 else 30 * prod((6 * t - 7) ** 2 for t in range(3, j + 1))
        assert count == expected, (j, count, expected)
        if j >= 3:
            factor = 121 * comb(j - 3, 0) + 168 * comb(j - 3, 1) + 72 * comb(j - 3, 2)
            assert factor == (6 * j - 7) ** 2
            assert count == factor * actual_counts[-1]
        actual_counts.append(count)
        print(f'TARGET j={j} n={len(p)} count={count} ideal_states={states}')
    # Compare integer coefficients in the binomial-basis identity after multiplying by 2:
    # 2*(121+168*x+72*x*(x-1)/2) = 2*(6*x+11)^2, x=j-3.
    assert (242, 336 - 72, 72) == (2 * 121, 2 * 2 * 6 * 11, 2 * 36)
    print('EXACT polynomial identity: 121*C(j-3,0)+168*C(j-3,1)+72*C(j-3,2)=(6j-7)^2')
    for j in range(2, 13):
        q = check_structure(j)
        if j <= 5:
            assert 2 * extension_count(q)[0] == actual_counts[j - 1]
        print(f'STRUCTURE j={j}: chains, seven-element head, remainder/gluing OK; |Q|={len(q)}')
    gluing_checks()
    print('Finite exact evidence only; the universal claim has no Lean proof in this check.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

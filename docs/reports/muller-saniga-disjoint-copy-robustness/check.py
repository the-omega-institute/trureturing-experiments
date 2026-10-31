#!/usr/bin/env python3
"""Exact checks for arXiv:2609.19862v1, Sections 4.2 and 6.

All arithmetic is integer or Fraction arithmetic. No numerical limit or Lean
settlement is claimed. Python's integers also represent finite vertex sets.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import prod


# Entries are (exponent of i, unsigned label), in the order I,X,Y,Z.
TABLE = (
    ((0, 0), (0, 1), (0, 2), (0, 3)),
    ((0, 1), (0, 0), (1, 3), (3, 2)),
    ((0, 2), (3, 3), (0, 0), (1, 1)),
    ((0, 3), (1, 2), (3, 1), (0, 0)),
)
LETTERS = 'IXYZ'


def word(s):
    return tuple(LETTERS.index(c) for c in s)


def multiply(a, b):
    pairs = [TABLE[x][y] for x, y in zip(a, b)]
    return sum(p[0] for p in pairs) % 4, tuple(p[1] for p in pairs)


def anticommutes(a, b):
    p, c = multiply(a, b)
    q, d = multiply(b, a)
    assert c == d
    assert (p - q) % 4 in (0, 2)
    return (p - q) % 4 == 2


def graph(words):
    return [sum(1 << j for j, b in enumerate(words) if anticommutes(a, b))
            for a in words]


def support(adj):
    """Independent nonempty sets with even overlap with every neighborhood."""
    result = []
    for h in range(1, 1 << len(adj)):
        if any((h >> v) & 1 and adj[v] & h for v in range(len(adj))):
            continue
        if all((ngh & h).bit_count() % 2 == 0 for ngh in adj):
            result.append(h)
    return result


def signed_contexts(words):
    adj = graph(words)
    assert len(set(adj)) == len(adj) and all(adj), 'reduced graph'
    hs = support(adj)
    assert all(any((h >> v) & 1 for h in hs) for v in range(len(words)))
    out = {}
    for h in hs:
        phase, acc = 0, (0,) * len(words[0])
        for v, a in enumerate(words):
            if (h >> v) & 1:
                p, acc = multiply(acc, a)
                phase = (phase + p) % 4
        assert not any(acc) and phase in (0, 2)
        out[h] = 1 if phase == 0 else -1
    return out


def score(hs, assignment):
    return sum(s * (-1 if (h & assignment).bit_count() % 2 else 1)
               for h, s in hs.items())


def all_scores(hs, n):
    return [score(hs, a) for a in range(1 << n)]


def combine(hs, ns):
    """Explicitly construct every tuple of contexts/empty choices, except all empty."""
    offsets, off = [], 0
    for n in ns:
        offsets.append(off)
        off += n
    out = {}
    choices = [[(0, 1)] + list(h.items()) for h in hs]
    for t in product(*choices):
        mask = sum(h << j for (h, _), j in zip(t, offsets))
        if mask:
            assert mask not in out
            out[mask] = prod(s for _, s in t)
    return out


def walsh_scores(hs, n):
    """Direct objective on every valuation via exact Walsh-Hadamard transform."""
    a = [0] * (1 << n)
    for h, sign in hs.items():
        a[h] = sign
    width = 1
    while width < len(a):
        for j in range(0, len(a), 2 * width):
            for i in range(j, j + width):
                x, y = a[i], a[i + width]
                a[i], a[i + width] = x + y, x - y
        width *= 2
    return a


def check_small_signed_systems():
    # Every signed set family on two vertices: each nonempty set absent, +, or -.
    families = [{h: s for h, s in zip(range(1, 4), t) if s}
                for t in product((0, 1, -1), repeat=3)]
    count = 0
    for left, right in product(families, repeat=2):
        joined = combine([left, right], [2, 2])
        assert len(joined) == (len(left) + 1) * (len(right) + 1) - 1
        for a, b in product(range(4), repeat=2):
            assert 1 + score(joined, a | (b << 2)) == (
                1 + score(left, a)) * (1 + score(right, b))
            count += 1
    print('All signed systems on two vertices, every pair and valuation:', count)


def check_small_graphs():
    # All 8 labelled simple graphs on three vertices, including degenerate graphs.
    edges = [(0, 1), (0, 2), (1, 2)]
    graphs = []
    for mask in range(8):
        a = [0] * 3
        for i, (u, v) in enumerate(edges):
            if (mask >> i) & 1:
                a[u] |= 1 << v
                a[v] |= 1 << u
        graphs.append(a)
    for left, right in product(graphs, repeat=2):
        actual = set(support(left + [x << 3 for x in right]))
        expected = set(combine([{h: 1 for h in support(left)},
                                {h: 1 for h in support(right)}], [3, 3]))
        assert actual == expected
    print('Complete-support component decomposition: all 64 graph pairs on 3+3 vertices')


def main():
    check_small_signed_systems()
    check_small_graphs()
    square = list(map(word, ['XI', 'IX', 'XX', 'IZ', 'ZI', 'ZZ', 'XZ', 'ZX', 'YY']))
    pentagram = list(map(word, ['XII', 'IXI', 'IIX', 'ZII', 'IZI', 'IIZ',
                                'XXX', 'XZZ', 'ZXZ', 'ZZX']))
    doily = [p for p in product(range(4), repeat=2) if any(p)]
    cases = [('pentagram', pentagram, 5, 1, 4),
             ('square', square, 6, 1, 5),
             ('doily', doily, 15, 3, 10)]
    bases = {}
    for name, words, m, d, factor in cases:
        hs = signed_contexts(words)
        assert len(hs) == m
        scores = all_scores(hs, len(words))
        assert scores == walsh_scores(hs, len(words))
        counts = Counter(scores)
        assert sum(scores) == 0
        assert (m - max(scores)) // 2 == d and (m - max(scores)) % 2 == 0
        assert max(1 + x for x in scores) == factor
        assert max(abs(1 + x) for x in scores) == factor
        assert all(-m <= x <= m - 2 for x in scores)
        if name == 'pentagram':
            assert set(scores) == {-5, -1, 3}
        if name == 'square':
            assert set(scores) == {-4, 0, 4}
        if name == 'doily':
            assert min(scores) == -9 and max(scores) == 9
        bases[name] = hs, len(words), scores
        print(f'{name}: vertices={len(words)} contexts={m} degree={d} '
              f'epsilon={Fraction(2*d, m)} valuations={len(scores)}')
        print('  exact score multiplicities:', sorted(counts.items()))
        # The finite score-value Cartesian products retain negative local factors.
        for k in range(1, 9):
            products = {1}
            for _ in range(k):
                products = {p * (1 + x) for p in products for x in counts}
            global_b = max(products) - 1
            global_m = (m + 1) ** k - 1
            assert global_b == factor ** k - 1
            assert (global_m - global_b) % 2 == 0
            assert 0 <= Fraction(global_b, global_m) <= Fraction(
                (m - 1) ** k - 1, global_m)
        print('  score-value product maxima and contextual bound checked for k=1..8')

    # Explicit Pauli tensor constructions and direct objectives on every valuation.
    for name, words in [('pentagram', pentagram), ('square', square)]:
        hs, n, scores = bases[name]
        joined = combine([hs, hs], [n, n])
        padded = [w + (0,) * len(w) for w in words] + [
            (0,) * len(w) + w for w in words]
        for h, sign in joined.items():
            phase, acc = 0, (0,) * len(padded[0])
            for v, w in enumerate(padded):
                if (h >> v) & 1:
                    p, acc = multiply(acc, w)
                    phase = (phase + p) % 4
            assert not any(acc) and phase == (0 if sign == 1 else 2)
        direct = walsh_scores(joined, 2 * n)
        for assignment, actual in enumerate(direct):
            a, b = assignment & ((1 << n) - 1), assignment >> n
            assert 1 + actual == (1 + scores[a]) * (1 + scores[b])
        print(f'{name} two copies: tensor signs and all {len(direct)} '
              f'global valuations agree; degree={(len(joined)-max(direct))//2}')

    # Fidelity with §4.2/Table 2 and the associated deposit's exact integer fields.
    for name, k, reported_m, reported_d in [
            ('pentagram', 3, 215, 76), ('square', 3, 342, 109), ('doily', 2, 255, 78)]:
        hs, n, scores = bases[name]
        joined = combine([hs] * k, [n] * k)
        products = {prod(1 + x for x in t) - 1
                    for t in product(set(scores), repeat=k)}
        m, d = len(joined), (len(joined) - max(products)) // 2
        assert (m, d) == (reported_m, reported_d)
        print(f'Source fidelity: {k} {name}s contexts={m} degree={d} '
              f'epsilon={Fraction(2*d, m)}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact exhaustive checks for Shapiro arXiv:2610.14373v1, Conjecture 1(ii).

Meet-in-the-middle enumerates every half-word, without assuming a multiplicity
pattern or either proposed lemma. It joins by literal permutation product and
full support. No external packages, random choices, or floating-point arithmetic.
"""
from collections import Counter, defaultdict
from itertools import combinations, product
from math import comb


def evaluate(n, edges, word):
    """Images of the literal product tau[e1] ... tau[ek] (rightmost acts first)."""
    p = list(range(n))
    for e in word:
        a, b = edges[e]
        p[a], p[b] = p[b], p[a]
    return tuple(p)


def inverse(p):
    q = [0] * len(p)
    for i, x in enumerate(p):
        q[x] = i
    return tuple(q)


def decode(code, base, length):
    out = [0] * length
    for j in range(length - 1, -1, -1):
        code, out[j] = divmod(code, base)
    assert code == 0
    return tuple(out)


def halves(n, edges, length, with_mask=False):
    """Exhaust all |E|^length words, storing codes by product (and support)."""
    bucket = defaultdict(lambda: defaultdict(list)) if with_mask else defaultdict(list)
    p = list(range(n))
    base = len(edges)
    leaves = 0

    def visit(left, code, mask):
        nonlocal leaves
        if left == 0:
            key = tuple(p)
            if with_mask:
                bucket[key][mask].append(code)
            else:
                bucket[key].append(code)
            leaves += 1
            return
        for e, (a, b) in enumerate(edges):
            p[a], p[b] = p[b], p[a]
            visit(left - 1, code * base + e, mask | (1 << e))
            p[a], p[b] = p[b], p[a]

    visit(length, 0, 0)
    assert leaves == base ** length
    return bucket, leaves


def path_words(n):
    edges = tuple((i, i + 1) for i in range(n - 1))
    length = 2 * n - 3
    a = length // 2
    left, nl = halves(n, edges, a)
    right, nr = halves(n, edges, length - a)
    target = (n - 1, *range(1, n - 1), 0)
    words = set()
    for p, codes in left.items():
        inv = inverse(p)
        # p * q = target, so q = p^-1 * target.
        q = tuple(inv[x] for x in target)
        for x in codes:
            for y in right.get(q, ()):
                w = decode(x, n - 1, a) + decode(y, n - 1, length - a)
                assert evaluate(n, edges, w) == target
                assert w not in words
                words.add(w)
    return words, nl + nr


def interleaving_words(n):
    """Construct the two shuffle phases, independently of permutation search."""
    words = set()
    m = n - 2
    for a in range(m + 1):
        b = m - a
        for before_one in combinations(range(m), a):
            before_one = set(before_one)
            for after_one in combinations(range(m), b):
                after_one = set(after_one)
                p = list(range(n))
                word = []
                for phase, positions in ((0, before_one), (1, after_one)):
                    if phase:
                        x, y = p.index(0), p.index(n - 1)
                        assert y == x + 1
                        p[x], p[y] = p[y], p[x]
                        word.append(x)
                    for j in range(m):
                        if j in positions:
                            e = p.index(0)
                        else:
                            e = p.index(n - 1) - 1
                        assert 0 <= e < n - 1
                        p[e], p[e + 1] = p[e + 1], p[e]
                        word.append(e)
                w = tuple(word)
                assert w not in words
                words.add(w)
    return words


def check_path(n):
    words, searched = path_words(n)
    expected = comb(2 * n - 4, n - 2)
    assert len(words) == expected
    assert words == interleaving_words(n)
    meetings = Counter()
    for w in words:
        counts = Counter(w)
        assert len(counts) == n - 1
        assert sorted(counts.values()) == [1] + [2] * (n - 2)
        p = list(range(n))
        last_inv = 0
        meeting = None
        for e in w:
            x, y = p[e], p[e + 1]
            assert x == 0 or y == n - 1
            if (x, y) == (0, n - 1):
                assert meeting is None
                meeting = e
            p[e], p[e + 1] = y, x
            inv = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
            assert inv == last_inv + 1
            last_inv = inv
        assert meeting is not None
        meetings[meeting] += 1
        assert counts[meeting] == 1
    assert [meetings[a] for a in range(n - 1)] == [comb(n - 2, a) ** 2 for a in range(n - 1)]
    print(f'PATH n={n} halves_examined={searched} endpoint_words={len(words)} '
          f'meeting_distribution={[meetings[a] for a in range(n - 1)]} multiplicities=OK', flush=True)
    return words


def cycle_words(n):
    edges = tuple((i, (i + 1) % n) for i in range(n))
    h = n - 1
    bucket, searched = halves(n, edges, h, with_mask=True)
    full = (1 << n) - 1
    words = set()
    # Each second half is the reversal of a word with the same product as
    # the first half. All letters are involutions. This tests every join.
    for groups in bucket.values():
        for mask, left_codes in groups.items():
            for mask2, right_codes in groups.items():
                if mask | mask2 != full:
                    continue
                for x in left_codes:
                    u = decode(x, n, h)
                    for y in right_codes:
                        w = u + decode(y, n, h)[::-1]
                        assert w not in words
                        words.add(w)
    identity = tuple(range(n))
    assert all(evaluate(n, edges, w) == identity and len(set(w)) == n for w in words)
    # Direct Cartesian enumeration additionally checks the joining algorithm.
    if n <= 5:
        direct = {w for w in product(range(n), repeat=2 * h)
                  if len(set(w)) == n and evaluate(n, edges, w) == identity}
        assert direct == words
    return words, searched


def check_cycle(n, paths):
    words, searched = cycle_words(n)
    length = 2 * n - 2
    cat = comb(2 * (n - 2), n - 2) // (n - 1)
    assert cat * (n - 1) == comb(2 * (n - 2), n - 2)
    expected = n * (n - 1) ** 2 * cat
    assert len(words) == expected
    triples = set()
    singleton_edges = Counter()
    for w in words:
        counts = Counter(w)
        assert sorted(counts.values()) == [1, 1] + [2] * (n - 2)
        for e, count in counts.items():
            if count != 1:
                continue
            singleton_edges[e] += 1
            pos = w.index(e)
            rotated = w[pos + 1:] + w[:pos + 1]
            assert rotated[-1] == e
            u = rotated[:-1]
            # C_n minus e has vertex order e+1,e+2,...,e (mod n).
            # Its path edge j is the cycle edge (e+1+j) mod n.
            path_word = tuple((x - e - 1) % n for x in u)
            assert all(x < n - 1 for x in path_word)
            assert path_word in paths
            data = (e, pos, path_word)
            assert data not in triples
            triples.add(data)
            cut = length - 1 - pos
            assert rotated[cut:] + rotated[:cut] == w
    # Verify the inverse map on the whole product E x positions x path words.
    constructed = 0
    for e in range(n):
        for pos in range(length):
            for u in paths:
                assert (e, pos, u) in triples
                base = tuple((e + 1 + x) % n for x in u) + (e,)
                cut = length - 1 - pos
                w = base[cut:] + base[:cut]
                assert w in words and w[pos] == e and w.count(e) == 1
                constructed += 1
    assert len(triples) == constructed == 2 * len(words) == n * length * len(paths)
    assert set(singleton_edges.values()) == {length * len(paths)}
    fidelity = 'source_verified_range' if n <= 6 else 'extension'
    print(f'CYCLE n={n} half_words_examined={searched} full_identity_words={len(words)} '
          f'formula={expected} marked_pairs={len(triples)} '
          f'marked_per_edge={length * len(paths)} multiplicities=OK bijection=OK {fidelity}', flush=True)


def main():
    print('EXACT exhaustive meet-in-the-middle; integer permutations and support masks; no floating point', flush=True)
    for n in range(3, 9):
        paths = check_path(n)
        check_cycle(n, paths)
    print('FIDELITY: all simple cycle sizes 3..6 stated as verified in Conjecture 1(ii); extended through 8')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

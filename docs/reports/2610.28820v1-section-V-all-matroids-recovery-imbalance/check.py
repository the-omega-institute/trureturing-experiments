#!/usr/bin/env python3
"""Exact checks for arXiv:2610.28820v1; Python standard library only.

Database inputs are the n <= 8 portion of matroid-database 0.3 (GPL-3.0).
This program checks finite identities and an injection mechanism, not a limit.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import hashlib
import lzma

HERE = Path(__file__).resolve().parent


def bits(x):
    while x:
        b = x & -x
        yield b
        x ^= b


def rank_table(n, bases):
    return tuple(max((s & b).bit_count() for b in bases) for s in range(1 << n))


def exchange_ok(bases):
    for a in bases:
        for b in bases:
            for e in bits(a & ~b):
                if not any((a ^ e ^ f) in bases for f in bits(b & ~a)):
                    return False
    return True


def source_counts(n, ranks):
    return tuple(tuple(sum(1 for s in range(1 << n)
                           if not s & (1 << i) and s.bit_count() == k
                           and ranks[s | (1 << i)] == ranks[s])
                       for k in range(n)) for i in range(n))


def profile(n, counts):
    return tuple(sum((Fraction(c[k], comb(n - 1, k)) for k in range(1, n)),
                     Fraction(0)) for c in counts)


def dual(n, bases):
    e = (1 << n) - 1
    return frozenset(e ^ b for b in bases)


def representatives(n, r):
    if r == 0:
        yield frozenset({0})
        return
    if r == n:
        yield frozenset({(1 << n) - 1})
        return
    rp = min(r, n - r)
    path = HERE / 'data' / f'allr{rp}n{n:02}.txt.xz'
    masks = [s for s in range(1 << n) if s.bit_count() == rp]
    with lzma.open(path, 'rt', encoding='ascii') as stream:
        for line in stream:
            line = line.strip()
            assert len(line) == len(masks) and set(line) <= {'0', '*'}
            bases = frozenset(s for s, c in zip(masks, line) if c == '*')
            assert bases
            yield dual(n, bases) if r > n - r else bases


def permute_mask(mask, p):
    return sum(1 << p[b.bit_length() - 1] for b in bits(mask))


def automorphisms(n, bases):
    """Exact clone-class enumeration, with pair-incidence pruning.

    A clone class consists of elements whose transposition fixes the basis set.
    Every automorphism permutes these classes; each class permutation has either
    zero or product(|class|!) lifts. Test one canonical lift explicitly.
    """
    def transposition(i, j):
        p = list(range(n))
        p[i], p[j] = p[j], p[i]
        return all(permute_mask(b, p) in bases for b in bases)
    classes = []
    for i in range(n):
        for c in classes:
            if transposition(i, c[0]):
                c.append(i)
                break
        else:
            classes.append([i])
    degree = [sum(bool(b & (1 << i)) for b in bases) for i in range(n)]
    pairs = [[sum((b & (1 << i) != 0) and (b & (1 << j) != 0)
                  for b in bases) for j in range(n)] for i in range(n)]
    signatures = [(len(c), degree[c[0]], tuple(sorted(pairs[c[0]]))) for c in classes]
    groups = {}
    for k, sig in enumerate(signatures):
        groups.setdefault(sig, []).append(k)
    groups = list(groups.values())
    quotient_order = 0
    for choices in product(*(permutations(g) for g in groups)):
        q = list(range(len(classes)))
        for g, image in zip(groups, choices):
            for i, j in zip(g, image):
                q[i] = j
        if any(pairs[c[0]][d[0]] != pairs[classes[q[i]][0]][classes[q[j]][0]]
               for i, c in enumerate(classes) for j, d in enumerate(classes)):
            continue
        p = list(range(n))
        for i, c in enumerate(classes):
            for u, v in zip(c, classes[q[i]]):
                p[u] = v
        if all(permute_mask(b, p) in bases for b in bases):
            quotient_order += 1
    order = quotient_order
    for c in classes:
        order *= factorial(len(c))
    assert order > 0 and factorial(n) % order == 0
    return order


def exhaustive_labelled(n):
    total = balanced = 0
    for r in range(n + 1):
        masks = [s for s in range(1 << n) if s.bit_count() == r]
        for code in range(1, 1 << len(masks)):
            bases = frozenset(masks[i] for i in range(len(masks)) if code & (1 << i))
            if exchange_ok(bases):
                total += 1
                a = profile(n, source_counts(n, rank_table(n, bases)))
                balanced += len(set(a)) <= 1
    return balanced, total


def main():
    print('Arithmetic: integers and fractions.Fraction; no floating-point evidence')
    print('Database: matroid-database 0.3; all ranks; isomorphism representatives weighted by n!/|Aut(M)|')
    print('Input SHA256:')
    for p in sorted((HERE / 'data').glob('*.xz')):
        print(p.name, hashlib.sha256(p.read_bytes()).hexdigest())
    iso_expected = [1, 2, 4, 8, 17, 38, 98, 306, 1724]
    totals = []
    fidelity = relaxations = balanced_pairs = 0
    by_size = []
    for n in range(9):
        m = b = iso = 0
        decoded_images = {}
        for r in range(n + 1):
            for bases in representatives(n, r):
                iso += 1
                assert exchange_ok(bases)
                ranks = rank_table(n, bases)
                counts = source_counts(n, ranks)
                values = profile(n, counts)
                is_balanced = len(set(values)) <= 1
                aut = automorphisms(n, bases)
                weight = factorial(n) // aut
                m += weight
                if is_balanced:
                    b += weight
                # Independent brute-force group order on all cases n <= 5.
                if n <= 5:
                    brute = sum(all(permute_mask(s, p) in bases for s in bases)
                                for p in permutations(range(n)))
                    assert aut == brute
                if not 2 <= r <= n - 2:
                    continue
                h_family = [s for s in range(1 << n) if s.bit_count() == r
                            and ranks[s] == r - 1
                            and all(ranks[s ^ e] == r - 1 for e in bits(s))
                            and all(ranks[s | e] == r for e in bits(((1 << n) - 1) ^ s))]
                # Source lemma labelled lemma:alpha: sparse-paving a_i(s).
                sparse = (all(ranks[s] == r - 1 for s in range(1 << n) if s.bit_count() == r - 1)
                          and all(ranks[s] == r for s in range(1 << n) if s.bit_count() == r + 1))
                if sparse:
                    nonbases = [s for s in range(1 << n) if s.bit_count() == r and s not in bases]
                    assert set(nonbases) == set(h_family)
                    degrees = [sum(bool(h & (1 << i)) for h in nonbases) for i in range(n)]
                    for i in range(n):
                        for k in range(1, n):
                            expected = (0 if k <= r - 2 else degrees[i] if k == r - 1
                                        else comb(n - 1, r) - len(nonbases) + degrees[i] if k == r
                                        else comb(n - 1, k))
                            assert counts[i][k] == expected
                    assert is_balanced == (len(set(degrees)) == 1)
                    fidelity += 1
                for h in h_family:
                    newbases = bases | {h}
                    assert exchange_ok(newbases)
                    newranks = rank_table(n, newbases)
                    assert all(newranks[s] == ranks[s] + (s == h) for s in range(1 << n))
                    newcounts = source_counts(n, newranks)
                    newvalues = profile(n, newcounts)
                    for i in range(n):
                        delta = (-Fraction(1, comb(n - 1, r - 1)) if h & (1 << i)
                                 else Fraction(1, comb(n - 1, r)))
                        assert newvalues[i] == values[i] + delta
                        changed = [(k, newcounts[i][k] - counts[i][k])
                                   for k in range(n) if newcounts[i][k] != counts[i][k]]
                        assert changed == ([(r - 1, -1)] if h & (1 << i) else [(r, 1)])
                    relaxations += 1
                    if is_balanced:
                        assert len(set(newvalues)) == 2
                        decoded_h = sum(1 << i for i in range(n) if newvalues[i] == min(newvalues))
                        assert decoded_h == h and newbases - {decoded_h} == bases
                        key = tuple(sorted(newbases))
                        assert key not in decoded_images or decoded_images[key] == (bases, h)
                        decoded_images[key] = (bases, h)
                        balanced_pairs += 1
        assert iso == iso_expected[n], (n, iso)
        totals.append((b, m))
        if n <= 5:
            assert exhaustive_labelled(n) == (b, m)
        by_size.append((n, iso, b, m, str(Fraction(b, m))))
        print('COUNT n=%d iso=%d balanced=%d all=%d ratio=%s' % by_size[-1], flush=True)
    print('FIDELITY sparse-paving formula and regularity equivalence:', fidelity, 'isomorphism representatives')
    print('RELAXATION rank-table, closure-count and rational-profile identities:', relaxations, 'pairs')
    print('INJECTION minimum-profile decoding and reconstruction:', balanced_pairs, 'balanced pairs')
    print('INDEPENDENT exhaustive labelled enumeration and automorphism cross-check: n<=5')
    print('Scope: database completeness is an external input; finite checks do not prove the asymptotic claim')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

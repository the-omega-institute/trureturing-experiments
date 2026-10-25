#!/usr/bin/env python3
"""Exact finite checks for arXiv:2610.29846v1, Section 7 item 11.

Uses Python integers only. Original packing counts come from all pairwise
window comparisons, independently of the weighted recurrence. Finite checks
do not establish the universal assertion or its recurrence bridge.
"""
from functools import lru_cache
from itertools import permutations
from math import comb

PATTERN = (1, 4, 3, 2, 5)
SOURCE_COUNTS = (
    1, 6, 71, 1266, 30206, 902796, 32420011, 1359292626,
    65164480466, 3515569641156, 210779736073446,
    13903319821066836, 1000559812125494076, 78012524487061315416,
    6550837823204594551731, 589404446176366002280146,
)


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def order_isomorphic(word, pattern):
    return all((word[a] < word[b]) == (pattern[a] < pattern[b])
               for a in range(len(pattern)) for b in range(len(pattern)))


def packs(word, pattern, overlap, occurrences):
    k = len(pattern)
    return (len(word) == k + (occurrences - 1) * (k - overlap)
            and all(order_isomorphic(word[t * (k - overlap):
                                          t * (k - overlap) + k], pattern)
                    for t in range(occurrences)))


def source_predecessors(n):
    """All strict comparisons from Definition 2.6, without a grid rewrite."""
    size = 5 + (n - 1) * 3
    pred = [0] * size
    for t in range(n):
        for a in range(5):
            for b in range(5):
                if PATTERN[a] < PATTERN[b]:
                    pred[3 * t + b] |= 1 << (3 * t + a)
    return pred


def linear_extensions(pred):
    full = (1 << len(pred)) - 1

    @lru_cache(None)
    def visit(ideal):
        if ideal == full:
            return 1
        total = 0
        remaining = full ^ ideal
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            vertex = bit.bit_length() - 1
            if pred[vertex] & ~ideal == 0:
                total += visit(ideal | bit)
        return total

    return visit(0), visit.cache_info().currsize


def grid_predecessors(n, v):
    """Two length-n chains, with v internal nodes in rung 1, then one each."""
    pred = [0] * (3 * n - 1 + v)
    for t in range(n - 1):
        pred[t + 1] |= 1 << t
        pred[n + t + 1] |= 1 << (n + t)
    nxt = 2 * n
    for t in range(n):
        chain = [t] + list(range(nxt, nxt + (v if t == 0 else 1))) + [n + t]
        nxt += v if t == 0 else 1
        for a, b in zip(chain, chain[1:]):
            pred[b] |= 1 << a
    check(nxt == len(pred), 'grid size')
    return pred


def weighted_rows(max_n, final_k):
    """B_1(k)=1; B_n(k)=sum_{i=1}^k i B_{n-1}(i+2), with exact integers."""
    rows = [None, [1] * (final_k + 2 * (max_n - 1) + 1)]
    for n in range(2, max_n + 1):
        cap = final_k + 2 * (max_n - n)
        row = [0] * (cap + 1)
        running = 0
        for i in range(1, cap + 1):
            running += i * rows[n - 1][i + 2]
            row[i] = running
        rows.append(row)
    return rows


def ballot(n, r):
    top = 2 * n + r - 2
    return comb(top, n - 1) - (comb(top, n - 2) if n >= 2 else 0)


def power_of_two(m):
    return m > 0 and m & (m - 1) == 0


def main():
    # Fidelity: the paper's example explicitly allows an additional occurrence.
    example = (1, 9, 2, 8, 3, 7, 4, 6, 5)
    check(packs(example, (1, 5, 2, 4, 3), 1, 2), 'Remark 2.7 prescribed copies')
    check(order_isomorphic(example[2:7], (1, 5, 2, 4, 3)), 'Remark 2.7 extra copy')
    print('FIDELITY Remark 2.7: starts 1, 5 prescribed; start 3 also present')

    rows = weighted_rows(256, 65)
    for n, expected in enumerate(SOURCE_COUNTS, 1):
        check(rows[n][3] == expected, f'source table term {n}')
    print('FIDELITY Section 7(11): all 16 displayed exact counts reproduced')
    print('SOURCE_COUNTS', ','.join(map(str, SOURCE_COUNTS)))

    for n in (1, 2):
        size = 3 * n + 2
        brute = sum(packs(w, PATTERN, 2, n) for w in permutations(range(1, size + 1)))
        check(brute == rows[n][3], f'actual permutation enumeration n={n}')
        print(f'PERMUTATIONS n={n} count={brute}')

    for n in range(1, 11):
        original, ideals = linear_extensions(source_predecessors(n))
        grid, _ = linear_extensions(grid_predecessors(n, 1))
        check(original == grid == rows[n][3], f'source/grid/recurrence n={n}')
        print(f'ORIGINAL_WINDOWS n={n} count={original} ideals={ideals} grid=match recurrence=match')

    for n in range(1, 6):
        for v in range(6):
            f, _ = linear_extensions(grid_predecessors(n, v))
            check(f == rows[n][v + 2], f'general grid recurrence n={n},v={v}')
    print('GRID_BRIDGE n=1..5 v=0..5: all 30 grid counts match B_n(v+2)')

    vectors = ((1, 1, 1, 1, 1), (0, 1, 3, 1, 0), (0, 1, 1, 1, 0))
    for n in range(1, 257):
        vector = vectors[min(n - 1, 2)]
        for k in range(66):
            check(rows[n][k] % 5 == vector[k % 5], f'mod-5 n={n},k={k}')
        for r in range(21):
            check(rows[n][2 * r + 1] % 2 == ballot(n, r) % 2,
                  f'ballot parity n={n},r={r}')
        cat = comb(2 * n, n) // (n + 1)
        check(ballot(n, 1) == cat, f'ballot/Catalan n={n}')
        check(cat % 2 == int(power_of_two(n + 1)), f'Catalan parity n={n}')
        want = 1 if power_of_two(n + 1) else 6
        check(rows[n][3] % 10 == want, f'mod-10 n={n}')
    print('MOD5 n=1..256 k=0..65: u -> v -> w -> w verified')
    print('BALLOT_PARITY n=1..256 r=0..20: all 5376 entries agree')
    print('CATALAN ballot identity and binary parity n=1..256: agree')
    print('MOD10 n=1..256: all residues agree; exceptions=1,3,7,15,31,63,127,255')
    print('SCOPE: original window counts checked through n=10; recurrence through n=256; no infinite proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

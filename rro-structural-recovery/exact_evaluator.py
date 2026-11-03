"""Direct source-capacity evaluator for arbitrary finite indexed masks.

Input JSON: {"cache": [{"0": 0, "7": 1}, {}, ...], "theta": "1/3"}.
Coordinates are zero-based within each original sample (left 0..6, right 7..14).
No theta returns exact rational coefficients Z_0,Z_1 of the affine expression.
Run: python3 exact_evaluator.py INPUT.json

This is an unrestricted finite mathematical evaluator with exponential loops,
not a practical performance guarantee for arbitrary input sizes. Its type-map
search prunes impossible assignments; the theory's arithmetic bound describes
its separately specified simple lexicographic evaluation schedule.
"""
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
from exact_checks import (M, literal, twin_geometry, capacity_embeddings,
                          merge_masks, label_dp, falling)


def partition_words(n, capacity):
    if n == 0:
        yield ()
        return
    word = [0]*n
    while True:
        yield tuple(word)
        for j in range(n-1, 0, -1):
            limit = min(capacity-1, max(word[:j])+1)
            if word[j] < limit:
                word[j] += 1
                word[j+1:] = [0]*(n-j-1)
                break
        else:
            return


def block_lists(word):
    b = [[] for _ in range(max(word, default=-1)+1)]
    for j, x in enumerate(word):
        b[x].append(j)
    return b


def normalize_cache(raw):
    cache = []
    for row in raw:
        entries = tuple(sorted((int(h), bit) for h, bit in row.items()))
        if not all(0 <= h < 15 and bit in (0, 1) for h, bit in entries):
            raise ValueError('Each cache coordinate must be 0..14 with a binary value.')
        cache.append(entries)
    return tuple(cache)


def normalizer(cache):
    n = len(cache)
    geometry = [twin_geometry(i, literal(i)) for i in (0, 1)]
    totals = [Fraction(0), Fraction(0)]
    for p in partition_words(n, 72):
        r = max(p, default=-1)+1
        l = label_dp(tuple(merge_masks(B, cache, 0, 7) for B in block_lists(p)), 7)
        if l == 0:
            continue
        for q in partition_words(n, 243):
            s = max(q, default=-1)+1
            right = label_dp(tuple(merge_masks(B, cache, 7, 8) for B in block_lists(q)), 8)
            if right == 0:
                continue
            edges = tuple(sorted(set(zip(p, q))))
            denominator = M**n*falling(128, r)*falling(256, s)
            for i in (0, 1):
                j = capacity_embeddings(r, s, edges, geometry[i])
                totals[i] += Fraction(j*l*right, denominator)
    return tuple(totals)


def main():
    obj = json.loads(Path(sys.argv[1]).read_text())
    cache = normalize_cache(obj['cache'])
    z0, z1 = normalizer(cache)
    result = {'Z0': str(z0), 'Z1': str(z1), 'symbolic': '(1-theta)*Z0+theta*Z1'}
    if 'theta' in obj:
        theta = Fraction(obj['theta'])
        if not 0 <= theta <= 1:
            raise ValueError('theta must belong to [0,1].')
        result['Z'] = str((1-theta)*z0+theta*z1)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()

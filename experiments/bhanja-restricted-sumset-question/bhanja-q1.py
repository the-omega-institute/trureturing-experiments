"""Bhanja (Integers 25 (2025) A12, arXiv:2307.03977v2) Section 4, Question 1:
is min |h^(<=r)A| = min |hA| over k-subsets A of every abelian group G?

h^(<=r)A = sums a_1+...+a_h of elements of A using at most r distinct elements.
Usage: python3 bhanja-q1.py
Checks:
  1. the family G = (Z/2)^m, h odd >= 3, 2^(m-1) < k < 2^m: every k-subset has h^(<=2)A = A and hA = G
     (exhaustive for m = 2, 3 and h = 3, 5);
  2. a survey of both minima for r = 2, h = 3 over small abelian groups (products of cyclic groups of order <= 12).
Exits nonzero if any expected reading fails.
"""
import itertools, sys

def group(mods):
    els = list(itertools.product(*[range(m) for m in mods]))
    add = lambda x, y: tuple((a + b) % m for a, b, m in zip(x, y, mods))
    zero = tuple(0 for _ in mods)
    return els, add, zero

def sumset(A, h, r, add, zero):
    out = set()
    for vals in itertools.combinations_with_replacement(A, h):
        if len(set(vals)) <= r:
            s = zero
            for v in vals:
                s = add(s, v)
            out.add(s)
    return out

ok = True
for m in (2, 3):
    els, add, zero = group((2,) * m)
    for h in (3, 5):
        for k in range(2 ** (m - 1) + 1, 2 ** m):
            for A in itertools.combinations(els, k):
                if sumset(A, h, 2, add, zero) != set(A) or sumset(A, h, h, add, zero) != set(els):
                    ok = False
            print(f"family (Z/2)^{m} h={h} k={k}: checked all {sum(1 for _ in itertools.combinations(els, k))} subsets")

def minima(mods, k, h, r):
    els, add, zero = group(mods)
    mu = min(len(sumset(A, h, h, add, zero)) for A in itertools.combinations(els, k))
    mur = min(len(sumset(A, h, r, add, zero)) for A in itertools.combinations(els, k))
    return mu, mur

survey = [(2, 2), (4,), (2, 4), (3,), (5,), (6,), (7,), (8,), (2, 2, 2), (3, 3), (9,), (10,), (12,), (2, 6)]
for mods in survey:
    n = 1
    for x in mods: n *= x
    diffs = []
    for k in range(1, n):
        mu, mur = minima(mods, k, 3, 2)
        if mu != mur:
            diffs.append((k, mur, mu))
    print(f"G=Z{mods} h=3 r=2: k with min|h^(<=2)A| < min|hA|: {diffs}")
sys.exit(0 if ok else 1)

#!/usr/bin/env python3
"""Checks for Mahalingam-Pandoh, arXiv:2108.02729v1, Conjecture 5.6 (maximum number of distinct
non-empty palindromic scattered subwords of a length-n word).  Usage: python3 check.py
Exits 0 and prints `bad= 0` when every check passes."""
import itertools, random, sys
from functools import lru_cache

bad = 0
def ok(cond, msg):
    global bad
    if not cond:
        bad += 1
        print("FAIL", msg)

def sp_brute(w):
    """Distinct non-empty palindromic subsequences, by enumerating all index subsets."""
    s = set()
    n = len(w)
    for mask in range(1, 1 << n):
        u = tuple(w[i] for i in range(n) if mask >> i & 1)
        if u == u[::-1]:
            s.add(u)
    return len(s)

def sp_dp(w):
    """Distinct palindromic subsequences by the standard interval recursion (counts the empty one, then subtract)."""
    w = tuple(w)
    @lru_cache(maxsize=None)
    def f(i, j):  # number of distinct palindromic subsequences of w[i..j], including empty
        if i > j:
            return 1
        total = 1
        for a in set(w[i:j + 1]):
            lo = w.index(a, i, j + 1)
            hi = j - w[i:j + 1][::-1].index(a)
            total += 1                      # "a"
            if hi > lo:
                total += f(lo + 1, hi - 1)  # "a p a" with p possibly empty
        return total
    return f(0, len(w) - 1) - 1

def M(n):
    if n == 0:
        return 0
    m = n // 2
    return 2 ** (m + 1) - 2 if n % 2 == 0 else 3 * 2 ** m - 2

def canonical_words(n):
    """Restricted growth strings: one representative per relabelling class."""
    def rec(prefix, mx):
        if len(prefix) == n:
            yield tuple(prefix); return
        for c in range(mx + 2):
            yield from rec(prefix + [c], max(mx, c))
    if n == 0:
        yield (); return
    yield from rec([0], 0)

# 1. the two counters agree
for n in range(1, 9):
    for w in canonical_words(n):
        ok(sp_brute(w) == sp_dp(w), f"brute/dp w={w}")

# 2. exhaustive maxima per (n, q) and the source Table 3
table3 = {1: [1], 2: [2, 2], 3: [3, 4, 3], 4: [4, 6, 6, 4], 5: [5, 9, 10, 8, 5], 6: [6, 12, 14, 14, 10, 6],
          7: [7, 17, 21, 22, 18, 12, 7], 8: [8, 22, 28, 30, 30, 22, 14, 8], 9: [9, 30, 41, 45, 46, 38, 26, 16, 9],
          10: [10, 38, 54, 60, 62, 62, 46, 30, 18, 10]}
for n in range(1, 13):
    best = [0] * (n + 1)
    for w in canonical_words(n):
        q = len(set(w))
        v = sp_dp(w)
        best[q] = max(best[q], v)
    if n in table3:
        ok(best[1:] == table3[n], f"table3 n={n}: {best[1:]}")
    top = max(best)
    ok(top == M(n), f"max n={n}: {top} vs {M(n)}")
    qs = [q for q in range(1, n + 1) if best[q] == top]
    want = [(n + 1) // 2] if n % 2 else [n // 2, n // 2 + 1]
    ok(qs == want, f"argmax q n={n}: {qs} vs {want}")
    # the source's proved formula for q >= n/2
    for q in range(1, n + 1):
        if 2 * q >= n:
            ok(best[q] == 2 ** (n - q) * (2 * q - n + 2) - 2, f"q>=n/2 formula n={n} q={q}")
    print(f"n={n} max={top} argmax_q={qs}")

# 3. Lemma 1 (outer-letter decomposition) on random words
rng = random.Random(20261011)
def lemma1_rhs(w):
    w = list(w)
    tot = 0
    for a in set(w):
        idx = [i for i, x in enumerate(w) if x == a]
        if len(idx) == 1:
            tot += 1
        else:
            tot += 2 + sp_dp(w[idx[0] + 1: idx[-1]])
    return tot
for _ in range(3000):
    n = rng.randint(1, 14); k = rng.randint(1, 5)
    w = [rng.randrange(k) for _ in range(n)]
    ok(sp_dp(w) == lemma1_rhs(w), f"lemma1 w={w}")

# 4. M recurrences, monotonicity, discrete convexity, M(t) >= t with equality iff t <= 2
for t in range(2, 60):
    ok(M(t) == 2 * M(t - 2) + 2, f"M rec t={t}")
for t in range(0, 60):
    ok(M(t + 1) - M(t) == 2 ** (t // 2), f"M increment t={t}")
    ok(M(t) >= t and ((M(t) == t) == (t <= 2)), f"M>=t t={t}")
    if t >= 1:
        ok(M(t + 1) - M(t) >= M(t) - M(t - 1), f"M convex t={t}")

# 5. Lemma 2 (convex endpoint matching) for f = M, exhaustively on small n
def matchings(pts):
    if not pts:
        yield []; return
    a = pts[0]
    for i in range(1, len(pts)):
        rest = pts[1:i] + pts[i + 1:]
        for m in matchings(rest):
            yield [(a, pts[i])] + m
for n in range(2, 11):
    for r in range(1, n // 2 + 1):
        bound = sum(M(n - 2 * i) for i in range(1, r + 1))
        for pts in itertools.combinations(range(1, n + 1), 2 * r):
            for m in matchings(list(pts)):
                ok(sum(M(h - l - 1) for l, h in m) <= bound, f"lemma2 n={n} m={m}")

# 6. the construction a_1..a_r c_1..c_k a_r..a_1 attains M(2r+k) for k in {0,1,2}
for r in range(0, 7):
    for k in (0, 1, 2):
        if 2 * r + k == 0:
            continue
        w = list(range(r)) + list(range(r, r + k)) + list(range(r))[::-1]
        ok(sp_dp(w) == M(2 * r + k), f"construction r={r} k={k}")

print("bad=", bad)
sys.exit(1 if bad else 0)

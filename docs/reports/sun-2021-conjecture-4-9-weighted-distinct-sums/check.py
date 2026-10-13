#!/usr/bin/env python3
"""Exact check of Z.-W. Sun, arXiv:1811.10503v6, Conjecture 4.9 (J. Algebraic Combin. 54 (2021)) in small cases.
S(A) = { sum_{k=1}^n k a_k : a_1..a_n distinct in A }; conjectured |S(A)| >= min(p(F), (|A|-n) n(n+1)/2 + n(n^2-1)/6 + 1)
for |A| >= n + [n = 3] and p(F) > n + 1. Reports: the counterexample n = 2, A = {0,1,2} over Q and over F_p for 5 <= p < 30;
all failing pairs (n, |A|) among arithmetic progressions for n <= 6, n <= |A| <= n + 4; and the minimum of |S(A)| over all
subsets of {0..9} for n = 2, 3. Exit 0 iff the only failing pair found is (2, 3)."""
import itertools, sys
def S(A, n, p=None):
    out = set()
    for a in itertools.permutations(A, n):
        v = sum((k + 1) * a[k] for k in range(n))
        out.add(v % p if p else v)
    return out
def bound(m, n):
    return (m - n) * n * (n + 1) // 2 + n * (n * n - 1) // 6 + 1
fails = []
print(f"Q: n=2 A={{0,1,2}} S={sorted(S((0,1,2),2))} |S|={len(S((0,1,2),2))} bound={bound(3,2)}")
for p in [5, 7, 11, 13, 17, 19, 23, 29]:
    s = S((0, 1, 2), 2, p)
    print(f"F_{p}: |S|={len(s)} bound=min({p},{bound(3,2)})={min(p,bound(3,2))}")
for n in range(1, 7):
    for m in range(n + (1 if n == 3 else 0), n + 5):
        s = len(S(tuple(range(m)), n)); b = bound(m, n)
        if s < b: fails.append((n, m, s, b))
for n in (2, 3):
    for m in range(n + (1 if n == 3 else 0), n + 3):
        best = min(len(S(A, n)) for A in itertools.combinations(range(10), m))
        if best < bound(m, n): fails.append((n, m, best, bound(m, n), 'subsets0..9'))
print("failing:", fails)
ok = all(f[0] == 2 and f[1] == 3 for f in fails) and any(f[0] == 2 and f[1] == 3 for f in fails)
print("bad=" + ("0" if ok else "1")); sys.exit(0 if ok else 1)

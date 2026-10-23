#!/usr/bin/env python3
"""Checks for Chamberland and Dilcher, arXiv:2510.26291v1, Conjecture 2.1.

d_n(l) = floor(sqrt(2ln)) - floor(sqrt((2l-1)n)) (2.3). Fidelity: Table 1 (n = 33) and the proved identity (2.4),
sum_{l <= (n-1)/2} d_n(l) = (n-1)/2. Conjecture: for each l with d_n(l) >= 2, the explicit block
{m - k + l : A_l + 1 <= k <= B_l - 1} (m = (n-1)/2, A_l = floor(sqrt((2l-1)n)), B_l = floor(sqrt(2ln))) is a run of
d_n(l) - 1 consecutive zeros inside 1..m; the blocks are pairwise disjoint and cover every zero. Exits 1 on failure.
"""
from math import isqrt
import sys

def d(n, l):
    return isqrt(2 * l * n) - isqrt((2 * l - 1) * n)

ok = True
row = [d(33, l) for l in range(1, 17)]
t1 = row == [3, 2, 2, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1]
print("fidelity Table 1 (n = 33):", row, t1)
ok &= t1
bad24 = [n for n in range(1, 4002, 2) if sum(d(n, l) for l in range(1, (n - 1) // 2 + 1)) != (n - 1) // 2]
print(f"fidelity (2.4) for odd n <= 4001: {len(bad24)} failures")
ok &= not bad24
N = 20001; fails = 0; blocks = 0
for n in range(3, N + 1, 2):
    m = (n - 1) // 2
    dv = [None] + [d(n, l) for l in range(1, m + 1)]
    used = set()
    for lam in range(1, m + 1):
        if dv[lam] >= 2:
            A, B = isqrt((2 * lam - 1) * n), isqrt(2 * lam * n)
            blk = [m - k + lam for k in range(A + 1, B)]
            good = len(blk) == dv[lam] - 1 and all(1 <= l <= m and dv[l] == 0 for l in blk) and not used & set(blk)
            good &= all(blk[i + 1] == blk[i] - 1 for i in range(len(blk) - 1))
            fails += not good; blocks += 1; used |= set(blk)
    fails += used != {l for l in range(1, m + 1) if dv[l] == 0}
print(f"conjecture blocks, odd n in [3, {N}]: {blocks} surplus indices, {fails} failures")
ok &= fails == 0
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)

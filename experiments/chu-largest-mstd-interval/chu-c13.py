"""Exact checks for Chu (Integers 20 (2020) A87) Conjecture 1.3: for n >= 16 the largest
MSTD (sum-dominant) subset of {0,...,n-1} has n - 7 elements.

1. Upper bound, exhaustive in the tested range: no subset with at least n - 6 elements
   (complement of at most 6 points) is MSTD, for 16 <= n <= NMAX.
2. Construction: A_n = I_n minus {3,5,6,n-7,n-6,n-5,n-3} is MSTD with |A_n+A_n| = 2n-2,
   |A_n-A_n| = 2n-3, for 16 <= n <= NMAX.
3. The six-point lemma |P|+|Q| <= 2|C| for all 1024 pairs X, Y of subsets of {0..5} containing 0.
Sets are bit masks; A+A and A-A are computed exactly by shifts.
"""
from itertools import combinations
import sys

NMAX = 40

def sums(mask, n):
    s = 0
    m = mask
    while m:
        b = m & -m
        i = b.bit_length() - 1
        s |= mask << i
        m ^= b
    return s

def diffs(mask, n):
    # differences shifted by n-1 to stay nonnegative
    rev = 0
    for i in range(n):
        if mask >> i & 1:
            rev |= 1 << (n - 1 - i)
    return sums_pair(mask, rev)

def sums_pair(a, b):
    s = 0
    m = b
    while m:
        bit = m & -m
        s |= a << (bit.bit_length() - 1)
        m ^= bit
    return s

def mstd(mask, n):
    return bin(sums(mask, n)).count("1") > bin(diffs(mask, n)).count("1")

ok = True
for n in range(16, NMAX + 1):
    full = (1 << n) - 1
    for k in range(0, 7):
        for H in combinations(range(n), k):
            mask = full
            for h in H:
                mask &= ~(1 << h)
            if mstd(mask, n):
                ok = False
                print("UPPER BOUND FAIL", n, H)
    A = full
    for h in (3, 5, 6, n - 7, n - 6, n - 5, n - 3):
        A &= ~(1 << h)
    s, d = bin(sums(A, n)).count("1"), bin(diffs(A, n)).count("1")
    if bin(A).count("1") != n - 7 or (s, d) != (2 * n - 2, 2 * n - 3):
        ok = False
        print("CONSTRUCTION FAIL", n, s, d)
print(f"upper bound and construction checked for 16 <= n <= {NMAX}: {ok}")

lem = True
I = 0b111111
for X in range(64):
    for Y in range(64):
        if not (X & 1 and Y & 1):
            continue
        P = bin(sums_pair(X, X) & I).count("1")
        Q = bin(sums_pair(Y, Y) & I).count("1")
        C = bin(sums_pair(X, Y) & I).count("1")
        if P + Q > 2 * C:
            lem = False
print("six-point lemma, all 1024 pairs:", lem)
ok &= lem
print("ALL_OK" if ok else "FAILED")
sys.exit(0 if ok else 1)

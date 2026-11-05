#!/usr/bin/env python3
from fractions import Fraction
from functools import lru_cache
from math import factorial

# The PARI recurrence on A380944, with its two boundary values.
@lru_cache(maxsize=None)
def b(v: int, k: int) -> int:
    if k == 0:
        return 1
    if v == 0:
        return 0
    if v & 1:
        return b(v // 2, k) + b(v // 2, k - 1)
    return ((v // 2).bit_count() - k + 1) * (
        b(v // 2, k) + b(v // 2, k - 1)
    )

def a(v: int) -> int:
    return b(v, v.bit_count())

@lru_cache(maxsize=None)
def stirling_second(n: int, k: int) -> int:
    if n == 0 and k == 0:
        return 1
    if n == 0 or k == 0:
        return 0
    return k * stirling_second(n - 1, k) + stirling_second(n - 1, k - 1)

def c(n: int, k: int) -> int:
    return sum(stirling_second(k + i, k) for i in range(max(0, n - k + 1)))

def conjectured_value(m: int, n: int) -> int:
    return (n + 1) ** m - n * factorial(n) * c(m, n + 1)

published = [
    1, 1, 2, 1, 3, 2, 3, 1, 4, 3, 5, 2, 9, 3, 4, 1,
    5, 4, 7, 3, 13, 5, 7, 2, 23, 9, 12, 3, 16, 4, 5, 1,
    6, 5, 9, 4, 17, 7, 10, 3, 31, 13, 18, 5, 25, 7, 9, 2,
    53, 23, 32, 9, 44, 12, 15, 3, 64, 16, 20, 4, 25, 5, 6, 1,
    7, 6, 11, 5, 21, 9, 13, 4, 39, 17, 24, 7, 34, 10,
]
assert len(published) == 78
assert [a(i) for i in range(len(published))] == published
print("fidelity: 78 published A380944 terms reproduced")

pairs = 0
for n in range(13):
    for m in range(21):
        v = 2**m * (2**n - 1)
        assert a(v) == conjectured_value(m, n), (m, n, v, a(v), conjectured_value(m, n))
        pairs += 1
print(f"claim: {pairs} pairs with 0 <= n <= 12 and 0 <= m <= 20")

# Exact inclusion-exclusion bridge used by the proposed proof.  Fractions keep
# the check explicitly in an exact field while all values remain integral.
lemmas = 0
for n in range(13):
    for j in range(21):
        lhs = sum(Fraction((-1) ** r * __import__('math').comb(n, r) * (n + 1 - r) ** j) for r in range(n + 1))
        rhs = Fraction(factorial(n) * stirling_second(j + 1, n + 1))
        assert lhs == rhs, (n, j, lhs, rhs)
        lemmas += 1
print(f"lemma-4: {lemmas} exact inclusion-exclusion instances")
print("ALL_CHECKS_OK")

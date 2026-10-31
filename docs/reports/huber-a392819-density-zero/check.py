#!/usr/bin/env python3
"""Exact finite checks for Huber's OEIS A392819 density-zero conjecture.

Uses only Python's standard library. Bounds are checked with integer arithmetic
or Fraction; finite checks do not establish a limit or a Lean theorem.
"""
from bisect import bisect_left
from fractions import Fraction
from math import gcd, isqrt, lcm, prod

LIMIT = 100_000
SOURCE_TERMS = [
    1, 4, 12, 48, 54, 60, 64, 90, 108, 120, 132, 180, 192, 224,
    240, 270, 336, 360, 432, 448, 528, 540, 660, 672, 720, 729,
    840, 918, 960, 1000, 1080, 1104, 1140, 1188, 1320, 1344,
    1400, 1458, 1530, 1680, 1728, 1740, 1836, 1980, 1984, 2112,
    2160, 2256, 2280, 2520, 2610,
]


def factor_trial(n):
    factors = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            factors.append((p, e))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def smallest_prime_factors(limit):
    spf = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if spf[p] == p:
            for n in range(p * p, limit + 1, p):
                if spf[n] == n:
                    spf[n] = p
    return spf


spf = smallest_prime_factors(LIMIT)
primes = [p for p in range(2, LIMIT + 1) if spf[p] == p]


def factor_sieve(n):
    factors = []
    while n > 1:
        p = spf[n]
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        factors.append((p, e))
    return factors


factors = [[]] + [factor_sieve(n) for n in range(1, LIMIT + 1)]
for n in range(1, 3001):
    assert factors[n] == factor_trial(n), ("independent factorization", n)
print("Independent trial-division/SPF factorization agreement: 1 <= n <= 3000")

is_term = [False] + [
    all(n % (p + e) == 0 for p, e in factors[n])
    for n in range(1, LIMIT + 1)
]
observed = [n for n in range(1, SOURCE_TERMS[-1] + 1) if is_term[n]]
assert observed == SOURCE_TERMS, ("source prefix", observed)
assert is_term[54] and is_term[840]
print(f"Fidelity: all {len(SOURCE_TERMS)} source terms, and no omissions, through 2610")
print("Fidelity: source examples 54 = 2*3^3 and 840 = 2^3*3*5*7 reproduced")

for n in range(1, LIMIT + 1):
    individual = is_term[n]
    assert individual == (n % lcm(*(p + e for p, e in factors[n])) == 0)
    if individual:
        for p, e in factors[n]:
            assert (n % (p * p) == 0) if e >= 2 else (n % (p * (p + 1)) == 0)
            if e == 1:
                assert gcd(p, p + 1) == 1
assert is_term[1] and not is_term[0]
assert is_term[12] and 12 % prod(p + e for p, e in factors[12]) != 0
assert is_term[54] and 54 % prod(e + 1 for _, e in factors[54]) != 0
print(f"Individual/lcm equivalence and prime-square-or-successor-product cover: all 1 <= n <= {LIMIT}")
print("Excluded reductions: 12 disproves product divisibility; 54 disproves containment in refactorable numbers")

prefix = [0] * (LIMIT + 1)
for n in range(1, LIMIT + 1):
    prefix[n] = prefix[n - 1] + is_term[n]

for y in (2, 3, 5, 10, 20, 50):
    constant = 2 ** bisect_left(primes, y)
    smooth_prefix = [0] * (LIMIT + 1)
    for n in range(1, LIMIT + 1):
        smooth = all(p < y for p, _ in factors[n])
        smooth_prefix[n] = smooth_prefix[n - 1] + smooth
        assert smooth_prefix[n] <= constant * isqrt(n), ("smooth bound", n, y)
        # Equivalent to c(n)/n <= constant*floor(sqrt(n))/n + 2/(y-1).
        assert (prefix[n] - constant * isqrt(n)) * (y - 1) <= 2 * n, ("count bound", n, y)
    for x in (1, 10, 100, 1000, 10_000, LIMIT):
        relevant = [p for p in primes if y <= p <= isqrt(x)]
        rough_bound = sum(x // (p * p) + x // (p * (p + 1)) for p in relevant)
        rough_count = sum(is_term[n] and any(p >= y for p, _ in factors[n]) for n in range(1, x + 1))
        assert rough_count <= rough_bound
        reciprocal_sum = sum((Fraction(1, p * p) + Fraction(1, p * (p + 1)) for p in relevant), Fraction())
        assert reciprocal_sum <= Fraction(2, y - 1)
        assert prefix[x] <= smooth_prefix[x] + rough_bound
    print(f"Cutoff y={y}: smooth and total count bounds hold for every 1 <= X <= {LIMIT}; divisor-cover counts checked at six X")

for x in (100, 1000, 10_000, LIMIT):
    print(f"Observed counting ratio at X={x}: {prefix[x]}/{x} = {Fraction(prefix[x], x)}")
print("Finite checks support the stated argument; they do not prove natural density zero")
print("ALL_CHECKS_OK")

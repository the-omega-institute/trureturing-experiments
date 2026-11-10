#!/usr/bin/env python3
"""Exact finite checks for Peter Bala's three conjectures in OEIS A099237.

The universal conclusion requires the mathematical argument; these computations
check source fidelity, all prime powers at most 200, and its intermediate facts.
Uses the Python standard library only. License: Apache-2.0.
"""
from math import comb, gcd, isqrt


def a(n: int) -> int:
    assert n >= 0
    return sum(comb(n * (n - j), j) for j in range(n + 1))


def prime(p: int) -> bool:
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def valuation(p: int, n: int) -> int:
    assert prime(p) and n > 0
    t = 0
    while n % p == 0:
        n //= p
        t += 1
    return t


# Published revision #21, OEIS A099237, DATA at indices 0 through 22.
source_terms = [
    1, 1, 3, 10, 45, 251, 1624, 11908, 97545, 880660, 8664546,
    92096731, 1050304775, 12778138842, 165033693175, 2253204163256,
    32401745953105, 489207829112931, 7733130368443057,
    127664099576228184, 2196149923000824756, 39288005888532450302,
    729576586543010269779,
]
assert [a(n) for n in range(len(source_terms))] == source_terms
print(f"FIDELITY: all {len(source_terms)} published terms (n=0..22) agree")
# A different indexing, also stated in the entry, checks endpoint conventions.
for n in range(201):
    assert a(n) == sum(comb(n * j, n - j) for j in range(n + 1))
assert a(0) == a(1) == 1
print("REINDEXING: literal formula equals sum_j C(n*j,n-j), n=0..200")

cases = []
term_count = 0
primes = [p for p in range(2, 201) if prime(p)]
for p in primes:
    # k=0: p**0=1, including every source clause's endpoint.
    assert (a(1) - 1) % p == 0
    k, n = 1, p
    while n <= 200:
        modulus = p ** (k + 1)
        assert n == p ** k
        # Iterated Frobenius: (1+X)^n = 1+X^n over F_p.
        assert all(comb(n, r) % p == 0 for r in range(1, n))
        assert comb(0, n) == 0
        total = 1
        for j in range(1, n):
            t = valuation(p, j)
            assert 0 <= t < k
            u = j // (p ** t)
            v = p ** (k - t) - u
            assert u > 0 and v > 0 and u % p != 0
            assert j == p ** t * u and n - j == p ** t * v
            assert (v + u) % p == 0
            m = n * (n - j)
            b = comb(m - 1, j - 1)
            c = comb(m, j)
            assert m > 0 and m == p ** (k + t) * v
            assert (b - (-1) ** (j - 1)) % p == 0
            assert j * c == m * b  # absorption, before cancellation
            assert u * c == n * v * b  # exact p**t cancellation
            assert gcd(u, modulus) == 1
            assert (u * c - n * u * (-1) ** j) % modulus == 0
            assert (c - n * (-1) ** j) % modulus == 0
            total += c
            term_count += 1
        assert total == a(n)
        alternating_sum = sum((-1) ** j for j in range(1, n))
        assert alternating_sum == (-1 if p == 2 else 0)
        assert (total - (1 + n * alternating_sum)) % modulus == 0
        if p == 2:
            assert (total - 1) % n == 0  # source dyadic clause
            assert (total - (1 - n)) % modulus == 0  # consumed stronger helper
        else:
            assert (total - 1) % modulus == 0  # both odd-prime clauses
        cases.append((p, k, n))
        print(f"PRIME_POWER: p={p} k={k} n={n} a(n)%p^(k+1)={total % modulus} OK")
        k += 1
        n *= p
print(f"ENDPOINT: k=0 checked for {len(primes)} primes; a(1)=1")
print(f"TOTAL: {len(cases)} prime-power cases, {term_count} interior terms, n<=200")
print("INTERMEDIATES: Frobenius coefficients, row residue, absorption, valuation and unit cancellation OK")
assert a(7) - 1 == 49 * 243 and (a(7) - 1) % (7 ** 3) != 0
print("SHARPNESS: a(7)-1=49*243 is not divisible by 7^3")
print("ALL_CHECKS_OK")

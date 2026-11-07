#!/usr/bin/env python3
"""Exact finite checks for Li's static exit sets and the square-deficit route.

This program checks finite identities and implications. It does not establish
an asymptotic density or refute a big-O assertion by finite computation.
"""
from collections import Counter
from fractions import Fraction
from math import isqrt, prod


def factor_trial(n):
    assert n >= 1
    factors = Counter()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] += 1
            n //= d
        d += 1
    if n > 1:
        factors[n] += 1
    return factors


def tau_trial(n):
    return prod(e + 1 for e in factor_trial(n).values())


def active_and_exit(k):
    if k == 1:
        return [], {0}
    active = [(j, tau_trial(k * k - j)) for j in range(1, 2 * k - 1)]
    active = [(j, t) for j, t in active if t >= j]
    return [j for j, _ in active], {t - j for j, t in active}


def transfer(k, r):
    x = k * k + r
    while x < (k + 1) ** 2:
        x += tau_trial(x)
    return x - (k + 1) ** 2


# Li, section 5, exact table: active deficits, distinct exits, transferred exits.
TABLE = {
    2: ([1, 2], {0, 1}, {0}),
    3: ([1, 2, 3], {0, 1, 3}, {2}),
    4: ([1, 2, 4], {2, 3}, {0, 7}),
    5: ([1, 2, 3, 4, 5], {0, 1, 7}, {2}),
    6: ([1, 2, 3, 4, 6], {1, 2, 3}, {1, 2}),
    7: ([1, 2, 3, 4, 5, 7], {0, 1, 2, 9}, {0, 2, 5}),
    8: ([1, 2, 4, 8], {0, 2, 5, 8}, {0, 3, 5}),
    9: ([1, 2, 3, 4, 5, 6, 9], {0, 1, 3, 5, 9}, {2, 8}),
    10: ([1, 2, 4, 10], {2, 4, 5, 8}, {1, 2, 15}),
    11: ([1, 2, 3, 4, 5, 7, 9], {1, 2, 15}, {0, 1, 2}),
    12: ([1, 2, 3, 4, 6, 8, 12], {0, 1, 2, 3, 8}, {3, 4}),
}
assert active_and_exit(1)[1] == {0}
for k, (expected_a, expected_e, expected_t) in TABLE.items():
    active, exits = active_and_exit(k)
    assert (active, exits) == (expected_a, expected_e), k
    assert {transfer(k, r) for r in exits} == expected_t, k
    print(f"FIDELITY k={k} active={active} E={sorted(exits)} A(E)={sorted(expected_t)}")
for k in range(2, 257):
    _, exits = active_and_exit(k)
    assert len(exits) >= 2
    if k % 8 in (3, 5):
        assert len(exits) >= 3
print("FIDELITY Proposition 7.5 checked for every 2 <= k <= 256")

# Independent divisor implementation, via a smallest-prime-factor sieve.
LIMIT = 100008
spf = list(range(LIMIT + 1))
for p in range(2, isqrt(LIMIT) + 1):
    if spf[p] == p:
        for n in range(p * p, LIMIT + 1, p):
            if spf[n] == n:
                spf[n] = p


def factor_sieve(n):
    assert 1 <= n <= LIMIT
    out = Counter()
    while n > 1:
        p = spf[n]
        out[p] += 1
        n //= p
    return out


for n in range(1, 4097):
    assert factor_sieve(n) == factor_trial(n)
    assert prod(e + 1 for e in factor_sieve(n).values()) == tau_trial(n)
print("INDEPENDENT divisor implementations agree for every 1 <= n <= 4096")

# Full finite residue enumeration: independence is across primes, not shifts.
P = (5, 7, 11)
M = prod(p * p for p in P)
b = {p: Fraction(p - 1, p * p) for p in P}
mu = sum(b.values(), Fraction())
variance = sum((bp * (1 - bp) for bp in b.values()), Fraction())
for h in (1, 2):
    counts = Counter()
    first = second = 0
    pairs = Counter()
    for k in range(M):
        indicators = {p: int((k - h) % p == 0 and (k - h) % (p * p) != 0) for p in P}
        x = sum(indicators.values())
        first += x
        second += x * x
        counts.update(p for p in P if indicators[p])
        for p in P:
            for q in P:
                if p < q and indicators[p] and indicators[q]:
                    pairs[p, q] += 1
    assert all(Fraction(counts[p], M) == b[p] for p in P)
    assert all(Fraction(pairs[p, q], M) == b[p] * b[q] for p in P for q in P if p < q)
    assert Fraction(first, M) == mu
    assert Fraction(second, M) - mu * mu == variance <= mu
print(f"RESIDUES M={M} h=1,2 mean={mu} variance={variance}; exact CRT pair counts agree")

# Check the square-deficit mechanism on every k <= 100000 for each J <= 4.
# n=(k-h)(k+h) is factored using the independently checked sieve, so factoring
# n itself or sieving an interval up to k^2 is unnecessary.
for J in range(1, 5):
    a = 1
    while 2 ** a <= J * J:
        a += 1
    k0 = max(J + 1, (J * J + 3) // 2)
    prime_successes = common_power_successes = 0
    first_prime = first_common = None
    for k in range(k0, 100001):
        counts, offsets = [], []
        for h in range(1, J + 1):
            left = factor_sieve(k - h)
            right = factor_sieve(k + h)
            exponents = left + right
            assert prod(p ** e for p, e in exponents.items()) == k * k - h * h
            t = prod(e + 1 for e in exponents.values())
            x = sum(e == 1 and p > 2 * J for p, e in left.items())
            for p, e in left.items():
                if e == 1 and p > 2 * J:
                    assert p not in right and exponents[p] == 1
            assert t % (2 ** x) == 0
            counts.append(x)
            offsets.append((t, t - h * h))
        all_primes = min(counts) >= a
        common_power = all(t % (2 ** a) == 0 for t, _ in offsets)
        if all_primes:
            prime_successes += 1
            first_prime = first_prime or k
            assert common_power
        if common_power:
            common_power_successes += 1
            first_common = first_common or k
            rs = []
            for h, (t, r) in enumerate(offsets, 1):
                assert 1 <= h * h <= 2 * k - 2
                assert t >= 2 ** a > J * J >= h * h
                assert r >= 0 and r % (2 ** a) == (-h * h) % (2 ** a)
                rs.append(r)
            assert len(set(rs)) == J
    if J <= 2:
        assert prime_successes > 0
    assert common_power_successes > 0
    print(f"SQUARES J={J} a={a} k0={k0} cutoff=100000 "
          f"prime_premise={prime_successes} first={first_prime} "
          f"common_power_premise={common_power_successes} first={first_common}")
print("SCOPE: finite identities and implications only; no numerical claim of asymptotic divergence")
print("ALL_CHECKS_OK")

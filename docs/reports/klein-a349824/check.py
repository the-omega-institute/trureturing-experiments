#!/usr/bin/env python3
"""Exact finite checks for Mike Klein's OEIS A349824 conjecture.

Source: https://oeis.org/A349824, revision 44. These finite checks are
experimental evidence; the universal descent argument requires a separate proof.
Uses only Python's integer arithmetic and standard library.
"""
from collections import Counter
from functools import lru_cache
from math import isqrt

BOUND = 200_000
SIEVE_LIMIT = 2 * BOUND + 20
SOURCE_VALUES = (
    0, 0, 2, 3, 8, 5, 10, 7, 18, 12, 14, 11, 21, 13, 18, 16, 32, 17,
    24, 19, 27, 20, 26, 23, 36, 20, 30, 27, 33, 29, 30, 31, 50, 28,
    38, 24, 40, 37, 42, 32, 44, 41, 36, 43, 45, 33, 50, 47, 55, 28,
    36, 40, 51, 53, 44, 32, 52, 44, 62, 59, 48, 61, 66, 39, 72, 36,
    48, 67, 63, 52, 42, 71, 60, 73,
)


def trial_factors(n):
    """Prime factor multiset by trial division; empty at 0 and 1."""
    assert n >= 0
    factors = []
    d = 2
    while d * d <= n and n > 1:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def trial_F(n):
    factors = trial_factors(n)
    return len(factors) * sum(factors)


spf = list(range(SIEVE_LIMIT + 1))
for p in range(2, isqrt(SIEVE_LIMIT) + 1):
    if spf[p] == p:
        for m in range(p * p, SIEVE_LIMIT + 1, p):
            if spf[m] == m:
                spf[m] = p
omega = [0] * (SIEVE_LIMIT + 1)
sopfr = [0] * (SIEVE_LIMIT + 1)
for n in range(2, SIEVE_LIMIT + 1):
    p = spf[n]
    omega[n] = omega[n // p] + 1
    sopfr[n] = sopfr[n // p] + p


@lru_cache(maxsize=None)
def F(n):
    return omega[n] * sopfr[n] if n <= SIEVE_LIMIT else trial_F(n)


def prime(n):
    if n <= SIEVE_LIMIT:
        return n >= 2 and spf[n] == n
    return n >= 2 and len(trial_factors(n)) == 1


def terminal(n):
    return n in (0, 27, 30, 28, 33) or prime(n)


def trajectory(n):
    seen = set()
    path = [n]
    while not terminal(path[-1]):
        assert path[-1] not in seen, ("unlisted cycle", path)
        seen.add(path[-1])
        path.append(F(path[-1]))
    return path


def main():
    # Fidelity checks use source data and source examples, not generated targets.
    assert tuple(F(n) for n in range(74)) == SOURCE_VALUES
    for n in range(10_001):
        assert F(n) == trial_F(n), ("independent factorization mismatch", n)
    example4 = (4, 8, 18, 24, 36, 40, 44, 45, 33, 28, 33, 28, 33, 28)
    example6 = (6, 10, 14, 18, 24, 36, 40, 44, 45, 33, 28, 33, 28, 33, 28)
    for example in (example4, example6):
        assert all(F(a) == b for a, b in zip(example, example[1:]))
    assert F(0) == F(1) == 0
    assert F(27) == 27 and F(30) == 30 and F(28) == 33 and F(33) == 28
    print("FIDELITY_OK: all 74 source DATA terms and both example trajectories")
    print("CROSSCHECK_OK: sieve and independent trial division agree for 0..10000")

    # Complete helper table for the general strong-induction argument.
    table = [F(n) for n in range(73)]
    assert max(table) <= 72
    base_paths = [trajectory(n) for n in range(73)]
    base_max = max(len(path) - 1 for path in base_paths)
    print("BASE_TABLE_F_0_TO_72:", ",".join(map(str, table)))
    print("BASE_TERMINATION_OK: n=0..72; maximum steps to target =", base_max)

    descent_count = 0
    for n in range(73, BOUND + 1):
        twice_prime = n % 2 == 0 and prime(n // 2)
        if not prime(n) and not twice_prime:
            assert F(n) < n, ("one-step descent failed", n, F(n))
            descent_count += 1
    print("ONE_STEP_DESCENT_OK: composite n=73..200000 excluding twice-primes; cases =", descent_count)

    odd_count = 0
    for q in range(39, BOUND + 1, 2):
        if not prime(q):
            assert F(2 * q) <= q + 15, ("odd composite bound failed", q, F(2 * q))
            odd_count += 1
    print("ODD_COMPOSITE_BOUND_OK: odd composite q=39..200000; cases =", odd_count)

    first_drop = Counter()
    for p in range(37, BOUND + 1):
        if prime(p):
            n = 2 * p
            x = n
            for j in range(1, 4):
                x = F(x)
                if x < n:
                    first_drop[j] += 1
                    break
            else:
                raise AssertionError(("three-step descent failed", p, x))
            if prime(p + 2):
                assert p % 3 == 2 and (p + 4) % 3 == 0 and not prime(p + 4)
    print("TWICE_PRIME_DESCENT_OK: prime p=37..200000; first-drop counts =", dict(sorted(first_drop.items())))

    endpoint_counts = Counter()
    max_steps = 0
    max_start = 0
    for n in range(BOUND + 1):
        path = trajectory(n)
        endpoint_counts["prime" if prime(path[-1]) else str(path[-1])] += 1
        steps = len(path) - 1
        if steps > max_steps:
            max_steps, max_start = steps, n
    assert sum(endpoint_counts.values()) == BOUND + 1
    print("ORBIT_CHECK_OK: every n=0..200000 reaches the source target set")
    print("ENDPOINT_COUNTS:", dict(sorted(endpoint_counts.items())))
    print("MAXIMUM_STEPS:", max_steps, "FIRST_START:", max_start)
    print("SCOPE: finite exact checks; no claim of a proof for all natural numbers")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

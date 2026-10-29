#!/usr/bin/env python3
"""Exact fidelity and finite checks for the OEIS A156834 preregistration.

Sources: https://oeis.org/A156834/internal (revision 20), and
el Abdalaoui, Dahmoune, Ziadi, arXiv:1301.3751v1, section 5.2.
The OEIS prefix and the paper's first eleven triangle rows are source data.
All computed evidence uses Python integer arithmetic; no external packages.
These finite checks do not prove an assertion for all natural numbers.
"""
from math import comb, gcd, isqrt

LIMIT = 10_000
TRIANGLE_LIMIT = 256
OEIS_PREFIX = (
    1, 2, 3, 5, 5, 12, 7, 17, 19, 30, 11, 63, 13, 56, 99, 89,
    17, 154, 19, 269, 237, 132, 23, 509, 301, 182, 379, 783, 29,
    1230, 31, 881, 813, 306, 2125, 2431, 37, 380, 1299, 4157,
    41, 4822, 43, 3695, 6175, 552, 47, 8529, 5587, 6266, 2787,
)
PAPER_ROWS = (
    (1,), (1, 1), (1, 0, 1), (1, 2, 0, 1), (1, 0, 0, 0, 1),
    (1, 3, 3, 0, 0, 1), (1, 0, 0, 0, 0, 0, 1),
    (1, 4, 0, 4, 0, 0, 0, 1), (1, 0, 6, 0, 0, 0, 0, 0, 1),
    (1, 5, 0, 0, 5, 0, 0, 0, 0, 1),
    (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def is_prime_trial(n):
    return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))


def totient_trial(n):
    remaining, result = n, n
    for p in range(2, isqrt(n) + 1):
        if remaining % p == 0:
            result = result // p * (p - 1)
            while remaining % p == 0:
                remaining //= p
    if remaining > 1:
        result = result // remaining * (remaining - 1)
    return result


def main():
    divisors = [[] for _ in range(LIMIT + 1)]
    phi = list(range(LIMIT + 1))
    for d in range(1, LIMIT + 1):
        for n in range(d, LIMIT + 1, d):
            divisors[n].append(d)
    for p in range(2, LIMIT + 1):
        if phi[p] == p:
            for n in range(p, LIMIT + 1, p):
                phi[n] = phi[n] // p * (p - 1)
    for n in range(1, LIMIT + 1):
        require(phi[n] == totient_trial(n), f"totient trial mismatch at {n}")
    for n in range(1, 201):
        count = sum(gcd(u, n) == 1 for u in range(1, n + 1))
        require(phi[n] == count, f"totient definition mismatch at {n}")

    # Construct Pascal's triangle by addition, then move column k to rows
    # k, 2k, 3k, ...; this directly implements the source's interleaving.
    triangle = [[0] * (n + 1) for n in range(TRIANGLE_LIMIT + 1)]
    pascal = [1]
    for i in range(1, TRIANGLE_LIMIT + 1):
        for k in range(1, i + 1):
            n = (i - k + 1) * k
            if n <= TRIANGLE_LIMIT:
                triangle[n][k] = pascal[k - 1]
        pascal = [1] + [pascal[j - 1] + pascal[j] for j in range(1, len(pascal))] + [1]
    for n in range(1, TRIANGLE_LIMIT + 1):
        for k in range(1, n + 1):
            formula = comb(n // k + k - 2, k - 1) if n % k == 0 else 0
            require(triangle[n][k] == formula, f"Proposition 6 mismatch at {(n, k)}")
    for n, expected in enumerate(PAPER_ROWS, 1):
        require(tuple(triangle[n][1:]) == expected, f"published triangle row mismatch at {n}")

    values = [0] * (LIMIT + 1)
    primes = composites = 0
    for n in range(1, LIMIT + 1):
        require(sum(phi[d] for d in divisors[n]) == n, f"totient sum mismatch at {n}")
        excess = 0
        strict_divisors = []
        for d in divisors[n]:
            b = comb(d + n // d - 2, d - 1)
            require(b >= 1, f"binomial positivity mismatch at {(n, d)}")
            require((b == 1) == (d == 1 or d == n), f"binomial equality mismatch at {(n, d)}")
            if 1 < d < n:
                require(b >= 2 and phi[d] > 0, f"proper divisor strictness mismatch at {(n, d)}")
                strict_divisors.append(d)
            values[n] += phi[d] * b
            excess += phi[d] * (b - 1)
        require(values[n] - n == excess, f"excess identity mismatch at {n}")
        if n <= TRIANGLE_LIMIT:
            matrix_value = sum(phi[k] * triangle[n][k] for k in range(1, n + 1))
            require(values[n] == matrix_value, f"matrix-vector definition mismatch at {n}")
        if n > 1:
            prime = is_prime_trial(n)
            require((values[n] == n) == prime, f"conjecture mismatch at {n}")
            if prime:
                primes += 1
                require(not strict_divisors and excess == 0, f"prime excess mismatch at {n}")
            else:
                composites += 1
                require(strict_divisors and excess > 0, f"composite strict inequality mismatch at {n}")
    require(tuple(values[1:len(OEIS_PREFIX) + 1]) == OEIS_PREFIX, "OEIS prefix mismatch")
    require(values[1] == 1 and not is_prime_trial(1), "excluded boundary mismatch")
    print(f"FIDELITY: all {len(OEIS_PREFIX)} OEIS A156834 published terms reproduced")
    print("PREFIX:", ",".join(map(str, values[1:len(OEIS_PREFIX) + 1])))
    print("FIDELITY: paper section 5.2 rows 1..11 reproduced")
    print(f"FIDELITY: Pascal interleaving equals Proposition 6 for all 1 <= k <= n <= {TRIANGLE_LIMIT}")
    print(f"FIDELITY: matrix-times-totient equals published divisor sum for n=1..{TRIANGLE_LIMIT}")
    print(f"EXACT: totient sieve equals trial factorization for n=1..{LIMIT}; coprime count for n=1..200")
    print(f"EXACT: totient divisor sum, binomial equality criterion, and excess identity hold for n=1..{LIMIT}")
    print(f"SETTLEMENT CHECK: n=2..{LIMIT}; primes={primes}; composites={composites}; mismatches=0")
    print("BOUNDARY: a(1)=1 and 1 is not prime; n>1 is essential")
    print("SCOPE: finite exact checks; the unbounded claim is reserved for a later Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

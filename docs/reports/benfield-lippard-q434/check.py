#!/usr/bin/env python3
"""Exact finite checks for Benfield--Lippard, arXiv:2407.20048v2, Question 4.34.

These checks corroborate a proposed universal classification; they do not prove it.
Only Python 3's standard library and integer arithmetic are used.
"""

from collections import Counter
from functools import lru_cache
from math import gcd, lcm


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def v2(n):
    require(n > 0, ("v2 domain", n))
    return (n & -n).bit_length() - 1


@lru_cache(maxsize=None)
def statistics(k, m):
    """Return (least pair period, least positive zero, zero count per period).

    A zero at index 0 is counted, and the zero at the return index is excluded.
    The invertible transition permutes the m*m pairs, so the search is bounded.
    The quotient period/rank is checked, rather than used to count zeros.
    """
    require(m > 0, ("modulus", m))
    x, y = 0, 1 % m
    first_zero = None
    zero_count = 0
    for n in range(1, m * m + 1):
        zero_count += x == 0
        x, y = y, (k * y + x) % m
        if x == 0 and first_zero is None:
            first_zero = n
        if (x, y) == (0, 1 % m):
            require(first_zero is not None, ("missing zero", k, m))
            require(n == first_zero * zero_count,
                    ("period product", k, m, n, first_zero, zero_count))
            return n, first_zero, zero_count
    raise AssertionError(("no return in finite state space", k, m))


def matrix_multiply(a, b, m):
    return ((a[0] * b[0] + a[1] * b[2]) % m,
            (a[0] * b[1] + a[1] * b[3]) % m,
            (a[2] * b[0] + a[3] * b[2]) % m,
            (a[2] * b[1] + a[3] * b[3]) % m)


def matrix_power(k, n, m):
    a = (k % m, 1 % m, 1 % m, 0)
    result = (1 % m, 0, 0, 1 % m)
    while n:
        if n % 2:
            result = matrix_multiply(result, a, m)
        a = matrix_multiply(a, a, m)
        n //= 2
    return result


def prime_divisors(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            yield p
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        yield n


def crosscheck_matrix_order(k, m, period):
    identity = (1 % m, 0, 0, 1 % m)
    require(matrix_power(k, period, m) == identity, ("matrix period", k, m))
    for p in prime_divisors(period):
        require(matrix_power(k, period // p, m) != identity,
                ("matrix minimality", k, m, period, p))


def check_case(k, m, counts):
    a = v2(m)
    s = m >> a
    r = v2(k)
    t = max(1, a - r + 1)
    period, rank, order = statistics(k, m)
    odd_period, odd_rank, odd_order = statistics(k, s)
    c = v2(odd_period)
    expected = 2 if c > t else 1
    require(order == expected, ("threshold", k, m, a, s, r, t, c, order))
    require(order in (1, 2), ("order range", k, m, order))
    require(rank % 2 == 0, ("even rank", k, m, rank))
    require(period == lcm(2 ** t, odd_period), ("period CRT", k, m))
    require(rank == lcm(2 ** t, odd_rank), ("rank CRT", k, m))
    require(gcd(2 ** a, s) == 1, ("coprime factors", m))
    if order == 1:
        require((period % 4 != 0) == (a <= r), ("second clause", k, m))
    if period % 4 == 0 and rank % 2 == 0:
        require((order == 2) == (c > t), ("first clause", k, m))
    if s > 1:
        b = v2(odd_rank)
        require((odd_order == 1 and b == c == 1)
                or (odd_order == 2 and b >= 1 and c == b + 1)
                or (odd_order == 4 and b == 0 and c == 2),
                ("odd component alternatives", k, s, odd_period, odd_rank))
    crosscheck_matrix_order(k, m, period)
    counts[order] += 1


def main():
    # Source fidelity: the displayed period table following Theorem 4.28.
    table = {2: (2, 4, 8, 16), 4: (2, 2, 4, 8),
             6: (2, 4, 8, 16), 8: (2, 2, 2, 4)}
    for k, periods in table.items():
        for a, expected in enumerate(periods, 1):
            require(statistics(k, 2 ** a) == (expected, expected, 1),
                    ("source table", k, a, expected))
    require(statistics(2, 8) == (8, 8, 1), "source counterexample omega_2(8)")
    require(statistics(2, 10)[2] == 2, "source omega_2(10)")
    require(statistics(2, 40)[2] == 1, "source omega_2(lcm(8,10))")
    print("FIDELITY: 16 published dyadic table entries and 3 mixed-order examples OK")

    dyadic_cases = 0
    for k in range(2, 130, 2):
        for a in range(1, 11):
            expected = 2 ** max(1, a - v2(k) + 1)
            require(statistics(k, 2 ** a) == (expected, expected, 1),
                    ("dyadic formula", k, a))
            dyadic_cases += 1
    print(f"DYADIC: {dyadic_cases} cases, even K=2..128, a=1..10 OK")

    dense = Counter()
    for k in range(2, 98, 2):
        for m in range(4, 1025, 4):
            check_case(k, m, dense)
    print(f"DENSE: {sum(dense.values())} pairs, even K=2..96, 4|m, 4<=m<=1024 OK")
    print(f"DENSE order counts: {dict(sorted(dense.items()))}")

    structured = Counter()
    for k in range(2, 66, 2):
        for a in range(2, 9):
            for s in range(1, 100, 2):
                check_case(k, 2 ** a * s, structured)
    print(f"STRUCTURED: {sum(structured.values())} triples, even K=2..64, a=2..8, odd s=1..99 OK")
    print(f"STRUCTURED order counts: {dict(sorted(structured.items()))}")

    # Independent recurrence over unreduced integers checks the valuation lemma.
    valuation_cases = 0
    for k in range(2, 130, 2):
        x, y = 0, 1
        for n in range(1, 129):
            x, y = y, k * y + x
            if n % 2:
                require(x % 2 == 1, ("odd-index parity", k, n))
            else:
                require(v2(x) == v2(k) + v2(n) - 1, ("U valuation", k, n))
            valuation_cases += 1
    print(f"INTEGER VALUATIONS: {valuation_cases} cases, even K=2..128, n=1..128 OK")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

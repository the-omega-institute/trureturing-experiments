#!/usr/bin/env python3
"""Exact finite checks for Kurkov's composite-term conjecture, OEIS A248667.

Source: https://oeis.org/A248667, revision 53 (2026-03-12).
All arithmetic is integral or fractions.Fraction; there is no floating point.
These are finite checks of a proposed general argument, not a proof.
"""
from fractions import Fraction
from math import gcd, isqrt, prod

POLYNOMIAL_BOUND = 100
MEMBERSHIP_BOUND = 3000
SOURCE_COEFFICIENTS = (
    (1,),
    (2, 2),
    (5, 12, 9),
    (16, 68, 112, 64),
    (65, 420, 1125, 1375, 625),
    (326, 2910, 11124, 21600, 20736, 7776),
)
SOURCE_TERMS = (
    1, 3, 7, 9, 11, 17, 19, 21, 23, 27, 29, 31, 33, 41, 43, 47, 49,
    51, 53, 57, 59, 61, 63, 67, 69, 71, 73, 77, 79, 81, 83, 87, 89,
    93, 97, 99, 101, 103, 107, 109, 113, 119, 121, 123, 127, 129,
    131, 133, 137, 139, 141, 147, 149, 151, 153, 157, 159, 161, 163, 167,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def multiply_linear(coefficients, constant, linear):
    """Multiply an ascending coefficient vector by constant + linear*X."""
    result = [0] * (len(coefficients) + 1)
    for degree, coefficient in enumerate(coefficients):
        result[degree] += constant * coefficient
        result[degree + 1] += linear * coefficient
    return result


def numerator_coefficients(n):
    """Direct sum of denominator tails, without using the content criterion."""
    require(n >= 1, "source domain n >= 1")
    result = [0] * n
    tail = [1]
    for j in range(n - 1, -1, -1):
        for degree, coefficient in enumerate(tail):
            result[degree] += coefficient
        if j:
            tail = multiply_linear(tail, j, n)
    return result


def evaluate(coefficients, x):
    value = 0
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def denominator_times_source_sum(n, x):
    """Independently evaluate the literal rational expression in the source."""
    denominator = 1
    total = Fraction(1)
    for i in range(1, n):
        denominator *= n * x + i
        total += Fraction(1, denominator)
    return denominator * total


def prime_divisors(n):
    """Trial division, returning each prime divisor once."""
    result = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return result


def main():
    factorial_sum = [1]
    for m in range(1, MEMBERSHIP_BOUND):
        factorial_sum.append(m * factorial_sum[-1] + 1)
    member = [False] + [
        gcd(n, factorial_sum[n - 1]) == 1
        for n in range(1, MEMBERSHIP_BOUND + 1)
    ]

    for n, expected in enumerate(SOURCE_COEFFICIENTS, 1):
        require(tuple(numerator_coefficients(n)) == expected,
                f"source polynomial mismatch n={n}")
    require(tuple(n for n in range(1, SOURCE_TERMS[-1] + 1) if member[n])
            == SOURCE_TERMS, "source initial terms mismatch")
    print("FIDELITY: source polynomials n=1..6 and all 60 displayed sequence terms match")

    polynomial_members = set()
    for n in range(1, POLYNOMIAL_BOUND + 1):
        coefficients = numerator_coefficients(n)
        require(coefficients[0] == factorial_sum[n - 1], f"constant n={n}")
        require(coefficients[-1] == n ** (n - 1), f"leading coefficient n={n}")
        require(all(c % n == 0 for c in coefficients[1:]),
                f"positive-degree divisibility n={n}")
        content = gcd(*coefficients)
        require((content == 1) == member[n], f"content criterion n={n}")
        if content == 1:
            polynomial_members.add(n)
        for x in (0, 1, 2):
            require(evaluate(coefficients, x) == denominator_times_source_sum(n, x),
                    f"literal numerator n={n}, x={x}")
    print(f"POLYNOMIALS: n=1..{POLYNOMIAL_BOUND}; exact content criterion, constant and leading coefficients")
    print(f"RATIONAL FIDELITY: n=1..{POLYNOMIAL_BOUND}, x in {{0,1,2}}; unreduced numerator identity")

    reset_pairs = 0
    divisor_pairs = 0
    composite_members = 0
    witnesses = []
    for n in range(1, MEMBERSHIP_BOUND + 1):
        primes = prime_divisors(n)
        for p in primes:
            require((factorial_sum[n - 1] - factorial_sum[p - 1]) % p == 0,
                    f"reset congruence n={n}, p={p}")
            reset_pairs += 1
        prime_support = all(factorial_sum[p - 1] % p != 0 for p in primes)
        require(member[n] == prime_support, f"prime-support criterion n={n}")
        if member[n]:
            for d in range(1, isqrt(n) + 1):
                if n % d == 0:
                    require(member[d] and member[n // d],
                            f"divisor closure n={n}, d={d}")
                    divisor_pairs += 1
            if n > 1 and not (len(primes) == 1 and primes[0] == n):
                p = primes[0]
                q = n // p
                require(p > 1 and q > 1 and member[p] and member[q] and p * q == n,
                        f"composite factorization n={n}")
                composite_members += 1
                if len(witnesses) < 10:
                    witnesses.append(f"{n}={p}*{q}")
                if n <= POLYNOMIAL_BOUND:
                    require(p in polynomial_members and q in polynomial_members,
                            f"independent polynomial factorization n={n}")
    print(f"RESET: {reset_pairs} prime-divisor pairs through n={MEMBERSHIP_BOUND}")
    print(f"PRIME SUPPORT: n=1..{MEMBERSHIP_BOUND}; {sum(member)} members")
    print(f"DIVISOR CLOSURE: {divisor_pairs} divisor pairs for those members")
    print(f"FACTORIZATION: all {composite_members} composite members through n={MEMBERSHIP_BOUND}")
    print("EXAMPLES: " + ", ".join(witnesses))
    require(member[9] and 9 == 3 * 3, "equal factors must be allowed")
    require(factorial_sum[2] == 5, "3 is an allowed prime")
    powers = []
    power = 1
    while power <= MEMBERSHIP_BOUND:
        require(member[power], f"power-of-three membership n={power}")
        powers.append(power)
        power *= 3
    print("ENUMERATION EXISTENCE CHECK: powers of 3 in range are members: "
          + ",".join(map(str, powers)))
    print("SCOPE: polynomial membership checked directly through 100; 101..3000 use Taylor's criterion")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite checks of the conjecture in OEIS A269454; no proof is claimed."""
from itertools import combinations
from math import gcd

PRIME_BOUND = 10_000
BRUTE_BOUND = 200
SOURCE_TERMS = (
    5, 11, 59, 83, 107, 179, 227, 347, 467, 563, 587, 1019, 1187,
    1283, 1307, 1523, 1619, 1907, 2027, 2099, 2459, 2579, 2819,
    2963, 3203, 3467, 3779, 3803, 3947, 4139, 4259, 4283, 4547,
    4787, 5099, 5387, 5483, 5507, 5939, 6659, 6779, 6827, 6899,
    7187, 7523,
)


def prime(n: int) -> bool:
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def prime_factors(n: int) -> tuple[int, ...]:
    assert n >= 1
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return tuple(factors)


def totient(n: int) -> int:
    assert n >= 1
    value = n
    for p in prime_factors(n):
        value -= value // p
    return value


def order_two(n: int) -> int:
    """Reduce Euler's exponent; certify the exact order by all prime factors."""
    assert n > 1 and gcd(2, n) == 1
    d = totient(n)
    assert pow(2, d, n) == 1
    for f in prime_factors(d):
        while d % f == 0 and pow(2, d // f, n) == 1:
            d //= f
    assert pow(2, d, n) == 1
    assert all(pow(2, d // f, n) != 1 for f in prime_factors(d))
    return d


def brute_order_two(n: int) -> int:
    """Visit successive residues without using totient or CRT."""
    assert n > 1 and gcd(2, n) == 1
    residue = 1
    for k in range(1, n + 1):
        residue = 2 * residue % n
        if residue == 1:
            return k
    raise AssertionError(("no positive order", n))


def safe_prime(p: int) -> bool:
    return prime(p) and prime((p - 1) // 2)


def main() -> None:
    # Source fidelity: every displayed %S/%T/%U term, with no omitted term.
    source_recomputed = tuple(
        p for p in range(2, SOURCE_TERMS[-1] + 1)
        if safe_prime(p) and p % 8 != 7
    )
    assert source_recomputed == SOURCE_TERMS
    assert safe_prime(5) and not safe_prime(3)
    print(f"FIDELITY: all {len(SOURCE_TERMS)} displayed A269454 terms through {SOURCE_TERMS[-1]} reproduced exactly")

    # Source fidelity: the stated universal half-totient upper bound.
    odd_primes = [p for p in range(3, 102, 2) if prime(p)]
    upper_pairs = 0
    for p, q in combinations(odd_primes, 2):
        n = p * q
        phi = totient(n)
        assert phi == sum(gcd(a, n) == 1 for a in range(1, n + 1))
        assert phi % 2 == 0
        assert pow(2, phi // 2, n) == 1
        assert order_two(n) <= phi // 2
        upper_pairs += 1
    print(f"FIDELITY: universal half-totient exponent and bound checked on {upper_pairs} distinct odd-prime pairs with p,q <= 101")

    allowed = [p for p in range(2, PRIME_BOUND + 1)
               if safe_prime(p) and p % 8 != 7]
    assert all(order_two(p) == p - 1 for p in allowed)
    checked = with_five = brute_checked = 0
    for p, q in combinations(allowed, 2):
        n = p * q
        phi = totient(n)
        actual = order_two(n)
        assert phi == (p - 1) * (q - 1)
        assert phi % 2 == 0
        assert actual == phi // 2, (p, q, actual, phi // 2)
        checked += 1
        with_five += p == 5 or q == 5
        if q <= BRUTE_BOUND:
            assert actual == brute_order_two(n)
            brute_checked += 1
    print(f"SETTLEMENT CHECK: {len(allowed)} admitted safe primes <= {PRIME_BOUND}; all {checked} unordered distinct pairs satisfy the conjectured identity")
    print(f"BOUNDARY: {with_five} pairs containing 5; {brute_checked} pairs with p,q <= {BRUTE_BOUND} independently checked by successive residues")
    assert order_two(5 * 11) == 20
    assert order_two(7 * 23) == 33 < totient(7 * 23) // 2
    print("EXAMPLES: ord_55(2)=20=phi(55)/2; excluded pair (7,23) has ord_161(2)=33<phi(161)/2=66")
    print("SCOPE: finite exact checks only; the universal assertion remains unproved by this program")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

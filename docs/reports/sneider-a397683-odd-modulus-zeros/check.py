#!/usr/bin/env python3
"""Exact checks for the odd-modulus zero conjecture in OEIS A397683.

Only Python 3's integer arithmetic and standard library are used. These finite
checks are evidence for a preregistration, not a proof of the conjecture.
"""
from functools import lru_cache
from math import gcd, lcm

SOURCE = "https://oeis.org/A397683/internal (revision 38)"
SOURCE_TERMS = (
    1, 0, 1, 0, 0, 4, 1, 1, 5, 6, 1, 8, 3, 0, 9, 6, 4, 3, 5, 1,
    7, 10, 1, 20, 0, 5, 15, 7, 5, 28, 11, 5, 9, 12, 9, 14, 8, 3, 11, 20,
    0, 40, 15, 9, 24, 17, 6, 11, 21, 4, 23, 30, 3, 52, 21, 4, 27, 19, 7, 15,
    36, 7, 15, 20, 7, 38, 33, 3, 43, 32, 21, 15, 41, 5, 51, 26, 5, 13, 19, 15,
)


@lru_cache(maxsize=None)
def factor(n):
    assert n >= 1
    out = []
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            out.append((d, e))
        d = 3 if d == 2 else d + 2
    if n > 1:
        out.append((n, 1))
    return tuple(out)


@lru_cache(maxsize=None)
def phi(n):
    z = n
    for p, _ in factor(n):
        z = z // p * (p - 1)
    return z


def order(u, m):
    """Reduce Euler's exponent by prime factors; require an actual unit."""
    assert m > 1 and gcd(u, m) == 1
    d = phi(m)
    for q, _ in factor(d):
        while d % q == 0 and pow(u, d // q, m) == 1:
            d //= q
    assert pow(u, d, m) == 1
    assert all(pow(u, d // q, m) != 1 for q, _ in factor(d))
    return d


def order_direct(u, m):
    """Independent repeated multiplication, used only for the small cross-check."""
    assert m > 1 and gcd(u, m) == 1
    v = 1
    for d in range(1, phi(m) + 1):
        v = v * u % m
        if v == 1:
            return d
    raise AssertionError((u, m))


@lru_cache(maxsize=None)
def f(m):
    assert m >= 1
    if m == 1:
        return 1
    orders = {u: order(u, m) for u in range(1, m) if gcd(u, m) == 1}
    return sum(orders[u] == orders[u + 1] for u in range(1, m - 1)
               if u in orders and u + 1 in orders)


def exceptional(m):
    fs = factor(m)
    return len(fs) == 1 and fs[0][0] in (3, 7)


@lru_cache(maxsize=None)
def primitive_pair(p):
    assert p >= 5 and p != 7 and factor(p) == ((p, 1),)
    for u in range(1, p - 1):
        if order(u, p) == order(u + 1, p) == p - 1:
            return u
    raise AssertionError(("no consecutive primitive pair", p))


@lru_cache(maxsize=None)
def lifted_pair(p):
    u = primitive_pair(p)
    bad_left = [t for t in range(p) if pow(u + t * p, p - 1, p * p) == 1]
    bad_right = [t for t in range(p) if pow(u + 1 + t * p, p - 1, p * p) == 1]
    assert len(bad_left) == len(bad_right) == 1
    for t in range(p):
        if t not in bad_left and t not in bad_right:
            return u + t * p
    raise AssertionError(("no joint lift", p))


def crt(residues):
    """Assemble one integer from pairwise coprime positive moduli."""
    x, modulus = 0, 1
    for u, q in residues:
        assert q > 1 and gcd(modulus, q) == 1
        t = ((u - x) * pow(modulus, -1, q)) % q
        x += modulus * t
        modulus *= q
        assert 0 <= x < modulus
    assert all(x % q == u % q for u, q in residues)
    return x, modulus


def construct(m):
    assert m > 1 and m % 2 == 1 and not exceptional(m)
    residues = []
    maximal_order = 1
    for p, e in factor(m):
        q = p ** e
        u = (4 if p == 3 else 3) if p in (3, 7) else (
            primitive_pair(p) if e == 1 else lifted_pair(p))
        residues.append((u, q))
        maximal_order = lcm(maximal_order, (p - 1) * p ** (e - 1))
    x, modulus = crt(residues)
    assert modulus == m and 1 <= x < m - 1
    assert gcd(x, m) == gcd(x + 1, m) == 1
    assert order(x, m) == order(x + 1, m) == maximal_order
    return x, maximal_order


def main():
    print("Source:", SOURCE)
    checked_orders = 0
    for m in range(3, 162, 2):
        for u in range(1, m):
            if gcd(u, m) == 1:
                assert order(u, m) == order_direct(u, m)
                checked_orders += 1
    print("Independent order algorithms agree:", checked_orders, "unit residues; odd 3 <= m <= 161")
    actual_terms = tuple(f(2 * n - 1) for n in range(1, len(SOURCE_TERMS) + 1))
    assert actual_terms == SOURCE_TERMS
    valid11 = [u for u in range(1, 10) if order(u, 11) == order(u + 1, 11)]
    assert valid11 == [3, 4, 6, 7] and f(11) == 4 and f(1) == 1
    print("Fidelity: all 80 published terms agree; f(11)=4 with x=", valid11, "; f(1)=1", sep="")

    zeros = []
    for m in range(1, 2002, 2):
        assert (f(m) == 0) == exceptional(m), (m, f(m))
        if f(m) == 0:
            zeros.append(m)
    print("Full counts: all 1001 positive odd m <= 2001; zero moduli:", zeros)

    prime_power_count = 0
    for p in range(3, 5001, 2):
        if factor(p) != ((p, 1),):
            continue
        q = p
        while q <= 5000:
            assert (f(q) == 0) == (p in (3, 7)), (p, q)
            prime_power_count += 1
            q *= p
    print("Source prime-power classification:", prime_power_count, "odd prime powers <= 5000")

    for e in range(1, 9):
        assert (order(4, 3 ** e), order(5, 3 ** e)) == (3 ** (e - 1), 2 * 3 ** (e - 1))
        assert (order(3, 7 ** e), order(4, 7 ** e)) == (6 * 7 ** (e - 1), 3 * 7 ** (e - 1))
    print("Exceptional local order formulas: both pairs at every exponent 1..8")

    lift_count = 0
    for p in range(5, 300, 2):
        if p == 7 or factor(p) != ((p, 1),):
            continue
        w = lifted_pair(p)
        for e in range(2, 6):
            assert order(w, p ** e) == order(w + 1, p ** e) == (p - 1) * p ** (e - 1)
            lift_count += 1
    print("Joint primitive lifting:", lift_count, "prime/exponent instances; p < 300, exponents 2..5")

    crt_count = 0
    pure_exceptions = []
    for m in range(3, 20000, 2):
        if exceptional(m):
            assert f(m) == 0
            pure_exceptions.append(m)
        else:
            construct(m)
            crt_count += 1
    print("All odd 1 <= m < 20000: f(1)=1;", crt_count,
          "CRT maximal-order witnesses; zero counts for", len(pure_exceptions), "pure powers")
    print("Pure powers below 20000:", pure_exceptions)
    assert construct(21) == (10, 6)
    print("Example modulo 21: x=10; both orders=6")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

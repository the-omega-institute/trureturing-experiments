#!/usr/bin/env python3
"""Exact finite checks for Ripà--Di Pietro, Conjecture 3.1, (3.4).

No numerical check certifies the eventual statement at all heights.
The modular tower evaluator does not use the proposed speed formula.
"""
from functools import cache
from math import gcd


@cache
def factor(n):
    assert n >= 1
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1 if p == 2 else 2
    if n > 1:
        out.append((n, 1))
    return tuple(out)


def valuation(n, p):
    assert n > 0 and p > 1
    e = 0
    while n % p == 0:
        n //= p
        e += 1
    return e


@cache
def tower_cap(a, h, cap):
    """min(T_h(a), cap), using only exact bounded integer arithmetic."""
    assert a >= 2 and h >= 0 and cap >= 1
    if h == 0:
        return min(1, cap)
    # Smallest exponent whose a-th power reaches cap.
    e, x = 0, 1
    while x < cap:
        x *= a
        e += 1
    exponent = tower_cap(a, h - 1, max(1, e))
    return cap if exponent >= e else pow(a, exponent)


@cache
def tower_mod(a, h, m):
    """T_h(a) mod m, via prime powers, Euler for units, and exact cutoffs."""
    assert a >= 2 and h >= 0 and m >= 1
    if m == 1:
        return 0
    if h == 0:
        return 1 % m
    residues = []
    for p, k in factor(m):
        pk = p ** k
        if a % p:
            phi = (p - 1) * p ** (k - 1)
            exponent = tower_mod(a, h - 1, phi)
            value = pow(a, exponent, pk)
        else:
            threshold = (k + valuation(a, p) - 1) // valuation(a, p)
            exponent = tower_cap(a, h - 1, threshold)
            value = 0 if exponent == threshold else pow(a, exponent, pk)
        residues.append((pk, value))
    # CRT, with inverses only between coprime factors.
    answer = sum(v * (m // pk) * pow(m // pk, -1, pk) for pk, v in residues) % m
    assert all(answer % pk == v for pk, v in residues)
    return answer


def local_slope(a, p):
    assert a > 1 and a % p
    return valuation(a * a - 1, 2) - 1 if p == 2 else valuation(a ** (p - 1) - 1, p)


def predicted_speed(r, a):
    units = [p for p, _ in factor(r) if a % p]
    assert units
    return min(local_slope(a, p) for p in units)


def capped_delta_order(a, p, h, precision):
    m = p ** precision
    residue = (tower_mod(a, h + 1, m) - tower_mod(a, h, m)) % m
    return precision if residue == 0 else valuation(residue, p)


def main():
    # Independent direct-integer checks of the modular evaluator, including nonunits.
    direct_count = 0
    for a in range(2, 9):
        t = 1
        for h in range(4):
            for m in range(1, 121):
                assert tower_mod(a, h, m) == t % m, (a, h, m)
                direct_count += 1
            if h < 3:
                t = a ** t
    for a in range(2, 15):
        for h in range(3):
            t = 1
            for _ in range(h):
                t = a ** t
            for cap in range(1, 80):
                assert tower_cap(a, h, cap) == min(t, cap)
    print(f"Evaluator cross-check: {direct_count} direct integer residues; capped towers also agree")

    # Source, journal p. 478: the seven residues shown for a = 27057, n = 25.
    expected = [
        27057,
        4777971988283521495227057,
        545271142050361495227057,
        6344863520050361495227057,
        6449963520050361495227057,
        1449963520050361495227057,
        1449963520050361495227057,
    ]
    got = [tower_mod(27057, h, 10 ** 25) for h in range(1, 8)]
    assert got == expected, (got, expected)
    print("Fidelity: journal's 27057 residues at heights 1..7 modulo 10^25 reproduced")
    for h, value in enumerate(got, 1):
        print(f"  h={h}: {value:025d}")
    tower_mod.cache_clear()

    # Full rectangular finite domain: all composite squarefree r <= 120,
    # all 2 <= a <= 150 with r not dividing a, all heights 8 <= h <= 12.
    # Orders >= precision are lower bounds only. A minimum strictly below
    # precision is an exact order, even when some prime branches saturate.
    radices = [r for r in range(2, 121)
               if len(factor(r)) >= 2 and all(e == 1 for _, e in factor(r))]
    primes = sorted({p for r in radices for p, _ in factor(r)})
    precision = 128
    cases = noncoprime = slope_checks = prime_formula_checks = 0
    strict_examples = []
    for a in range(2, 151):
        orders = {p: {h: capped_delta_order(a, p, h, precision) for h in range(7, 13)}
                  for p in primes}
        # Source's established prime-radix formulas, checked independently.
        for p in primes:
            if a % p:
                c = local_slope(a, p)
                for h in range(8, 13):
                    assert orders[p][h] < precision
                    assert orders[p][h] - orders[p][h - 1] == c, (a, p, h)
                    prime_formula_checks += 1
        for r in radices:
            if a % r == 0:
                continue
            nr = {h: min(orders[p][h] for p, _ in factor(r)) for h in range(7, 13)}
            assert max(nr.values()) < precision, (r, a, 'increase precision')
            v = predicted_speed(r, a)
            for h in range(8, 13):
                assert nr[h] - nr[h - 1] == v, (r, a, h, nr, v)
                slope_checks += 1
            factor_min = min(predicted_speed(r, p) for p, _ in factor(a))
            assert v >= factor_min, (r, a, v, factor_min)
            if v > factor_min and len(strict_examples) < 6:
                strict_examples.append((r, a, v, factor_min))
            cases += 1
            noncoprime += gcd(r, a) > 1
        tower_mod.cache_clear()
    print(f"Domain: {len(radices)} composite squarefree radices 2..120; bases 2..150 with r not dividing a")
    print(f"Expected identity: {slope_checks} exact observed slope checks at heights 8..12, precision {precision}")
    print(f"Source prime-radix formulas: {prime_formula_checks} exact observed slope checks")
    print(f"Inequality: {cases} cases, including {noncoprime} noncoprime cases")
    print(f"Strict examples (r, a, V(a), minimum factor speed): {strict_examples}")
    print("Finite agreement supports the proposed general argument; it is not a proof of eventual stabilization")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()

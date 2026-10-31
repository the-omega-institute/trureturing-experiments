#!/usr/bin/env python3
"""Exact checks of Davis Smith's second conjecture in OEIS A324017.

Source: https://oeis.org/A324017/internal, revision 91.
The modular tower evaluator uses prime powers, Euler reduction on units,
valuation cutoffs on nonunits, and CRT. It never uses the conjectured residue.
Only Python integers are used. Finite checks are evidence, not a proof.
"""
from functools import lru_cache


@lru_cache(maxsize=None)
def prime_powers(modulus):
    assert modulus >= 1
    remaining = modulus
    factors = []
    p = 2
    while p * p <= remaining:
        if remaining % p == 0:
            power = 1
            degree = 0
            while remaining % p == 0:
                remaining //= p
                power *= p
                degree += 1
            factors.append((p, degree, power))
        p += 1
    if remaining > 1:
        factors.append((remaining, 1, remaining))
    return tuple(factors)


@lru_cache(maxsize=None)
def tower_capped(base, height, cap):
    """Return min(base^^height, cap), without constructing a large tower."""
    assert base >= 1 and height >= 0 and cap >= 0
    if cap <= 1:
        return cap
    if height == 0 or base == 1:
        return 1
    threshold, value = 0, 1
    while value < cap:
        value *= base
        threshold += 1
    exponent = tower_capped(base, height - 1, threshold)
    if exponent == threshold:
        return cap
    return min(pow(base, exponent), cap)


@lru_cache(maxsize=None)
def tower_mod(base, height, modulus):
    """Return base^^height modulo modulus, by an independent generic algorithm."""
    assert base >= 1 and height >= 0 and modulus >= 1
    if modulus == 1:
        return 0
    if height == 0:
        return 1 % modulus
    result = 0
    for p, degree, power in prime_powers(modulus):
        if base % p:
            phi = (power // p) * (p - 1)
            exponent = tower_mod(base, height - 1, phi)
            residue = pow(base, exponent, power)
        else:
            valuation, remaining = 0, base
            while remaining % p == 0:
                valuation += 1
                remaining //= p
            threshold = (degree + valuation - 1) // valuation
            exponent = tower_capped(base, height - 1, threshold)
            residue = 0 if exponent == threshold else pow(base, exponent, power)
        complement = modulus // power
        result += residue * complement * pow(complement, -1, power)
    return result % modulus


def tower_exact(base, height):
    """Literal source recursion, called only on deliberately bounded cases."""
    value = 1
    for _ in range(height):
        value = pow(base, value)
    return value


def main():
    # Cross-check the generic evaluator against fully constructed integers.
    generic_cases = 0
    for base in range(1, 8):
        for height in range(4):
            value = tower_exact(base, height)
            for modulus in range(1, 129):
                assert tower_mod(base, height, modulus) == value % modulus
                generic_cases += 1
            for cap in range(65):
                assert tower_capped(base, height, cap) == min(value, cap)
    print(f"EVALUATOR: {generic_cases} exact tower/modulus comparisons; capped checks OK")

    # Fidelity: first three displayed rows, with no claim-specific formula.
    source_rows = (
        (1, 3, 5, 7, 9, 11, 13, 15),
        (1, 11, 29, 55, 89, 131, 181, 239),
        (1, 59, 29, 119, 289, 563, 965, 1519),
    )
    for height, expected in enumerate(source_rows, 1):
        actual = tuple(tower_mod(2*n - 1, height, (2*n)**height) for n in range(1, 9))
        assert actual == expected, (height, actual, expected)
        print(f"FIDELITY A({height},1..8): {actual}")

    # If X = q^2*k + r with 0 <= r < q^2, then (X//q)%q = r//q.
    # This extracts the digit from a generic modular computation.
    for n, height, expected in ((7, 3, 12), (8, 4, 14)):
        q = 2*n
        residue = tower_mod(q - 1, height, q*q)
        actual = residue // q
        assert actual == expected
        print(f"FIDELITY quotient example: n={n} m={height} residue_mod_q2={residue} digit={actual}")

    # Literal quotient checks, avoiding the modular evaluator altogether.
    exact_cases = 0
    for n in range(2, 201):
        q = 2*n
        value = tower_exact(q - 1, 2)
        assert (value // q) % q == q - 2
        exact_cases += 1
    for n in range(2, 5):
        q = 2*n
        value = tower_exact(q - 1, 3)
        assert (value // q) % q == q - 2
        exact_cases += 1
    print(f"CLAIM literal integers: {exact_cases} cases (n=2..200,m=2; n=2..4,m=3)")

    modular_cases = 0
    for n in range(2, 201):
        q = 2*n
        for height in range(2, 13):
            residue = tower_mod(q - 1, height, q*q)
            assert residue // q == q - 2, (n, height, residue)
            assert residue == q*q - q - 1, (n, height, residue)
            modular_cases += 1
    print(f"CLAIM generic modular evaluator: {modular_cases} cases (n=2..200,m=2..12)")

    lemma_cases = 0
    for q in range(4, 401, 2):
        for exponent in range(1, 402, 2):
            assert pow(q - 1, exponent, q*q) == (q*exponent - 1) % (q*q)
            lemma_cases += 1
    print(f"OUTLINE odd-exponent binomial residue: {lemma_cases} cases")

    boundary_failures = []
    for n in range(2, 21):
        q = 2*n
        actual = (tower_exact(q - 1, 1) // q) % q
        assert actual == 0 and actual != q - 2
        boundary_failures.append(n)
    print(f"BOUNDARY m=1 fails for n=2..20: {len(boundary_failures)} cases")
    print("LIMIT: finite exact checks only; the unbounded claim has no Lean proof here")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()

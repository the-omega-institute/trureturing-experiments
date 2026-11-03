#!/usr/bin/env python3
"""Exact finite checks for Bayart--Matheron Question 5.7 (arXiv:2605.23736v1).

This is a preregistration check, not a proof of an infinite-dimensional claim.
Only Python integers and fractions.Fraction enter assertions.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod


def decode(k, radices):
    digits = []
    for radix in radices:
        digits.append(k % radix)
        k //= radix
    return tuple(digits)


def encode(x, radices):
    result, place = 0, 1
    for digit, radix in zip(x, radices):
        result += digit * place
        place *= radix
    return result


def add_with_carry(x, y, radices):
    carry, z = 0, []
    for a, b, radix in zip(x, y, radices):
        total = a + b + carry
        z.append(total % radix)
        carry = total // radix
    return tuple(z)


def odometer(x, radices):
    digits = list(x)
    for i, radix in enumerate(radices):
        if digits[i] + 1 < radix:
            digits[i] += 1
            break
        digits[i] = 0
    return tuple(digits)


def source_fidelity():
    """Source Fact 1.2 and the basic-cylinder preimage formula preceding (1.2)."""
    cases = 0
    for length in range(1, 4):
        for radices in product((2, 3, 4), repeat=length):
            modulus = prod(radices)
            for j in range(modulus):
                x = decode(j, radices)
                assert encode(x, radices) == j
                state = x
                for k in range(2 * modulus):
                    assert state == add_with_carry(x, decode(k, radices), radices)
                    assert encode(state, radices) == (j + k) % modulus
                    state = odometer(state, radices)
                    cases += 1
                predecessor = decode((j - 1) % modulus, radices)
                if j == 0:
                    assert predecessor == tuple(r - 1 for r in radices)
                else:
                    ell = next(i for i, a in enumerate(x) if a > 0)
                    expected = tuple(radices[i] - 1 if i < ell else
                                     x[i] - 1 if i == ell else x[i]
                                     for i in range(length))
                    assert predecessor == expected
    print(f"FIDELITY Fact 1.2 mixed-radix iterate cases: {cases}; cylinder preimages: OK")


def probabilities(length):
    # Source indexing i >= 1: a_i = 1/(4(i+1)).
    return [(Q(1, 2) + Q(1, 4 * (i + 1)),
             Q(1, 2) - Q(1, 4 * (i + 1))) for i in range(1, length + 1)]


def prefix_checks():
    """Exhaust all atoms and all cyclic powers, not just a chosen pair."""
    ratios_checked = 0
    for length in range(1, 9):
        marginals = probabilities(length)
        modulus = 2 ** length
        q = [prod(marginals[i][bit] for i, bit in
                  enumerate(decode(j, (2,) * length))) for j in range(modulus)]
        assert sum(q) == 1 and min(q) > 0
        one_step = max(q[(j - 1) % modulus] / q[j] for j in range(modulus))
        assert one_step == Q(5, 3)
        all_powers = max(q[(j - k) % modulus] / q[j]
                         for j in range(modulus) for k in range(modulus))
        ratios_checked += modulus ** 2
        assert all_powers == q[0] / q[-1] == Q(2 * length + 3, 3)
        # The test cylinder E_N consists of all ones, k_N = 2^N - 1.
        assert (modulus - 1 - (modulus - 1)) % modulus == 0
        uniform = Q(1, modulus)
        moment_by_atoms = sum(uniform ** 2 / mass for mass in q)
        moment_by_product = prod(Q(1, 1) / (1 - 4 * Q(1, 4 * (i + 1)) ** 2)
                                 for i in range(1, length + 1))
        assert moment_by_atoms == moment_by_product <= Q(4, 3)
        print(f"N={length}: one-step={one_step}; all-power maximum={all_powers}; "
              f"inverse-density L2 moment={moment_by_atoms}")
    print(f"COUNTEREXAMPLE all-atom/all-power ratios checked: {ratios_checked}")


def analytic_ingredients():
    """Check exact telescoping and product estimates through N=128.

    Their universal proofs and the passage to L2 are given in the draft;
    finite verification is not represented as an infinite proof.
    """
    ratio_product = Q(1)
    lower_product = Q(1)
    square_sum = Q(0)
    for i in range(1, 129):
        a = Q(1, 4 * (i + 1))
        assert 0 < a <= Q(1, 8)
        assert Q(3, 8) <= Q(1, 2) - a < Q(1, 2) + a <= Q(5, 8)
        factor = (1 + 2 * a) / (1 - 2 * a)
        assert factor == Q(2 * i + 3, 2 * i + 1)
        ratio_product *= factor
        assert ratio_product == Q(2 * i + 3, 3)
        preceding_inverse = Q(3, 2 * i + 1)
        derivative = preceding_inverse * factor
        assert derivative == Q(3 * (2 * i + 3), (2 * i + 1) ** 2) <= Q(5, 3)
        u = 4 * a ** 2
        assert u <= Q(1, 4 * i * (i + 1))
        square_sum += u
        lower_product *= 1 - u
        assert square_sum <= Q(1, 4) * (1 - Q(1, i + 1))
        assert lower_product >= 1 - square_sum > Q(3, 4)
        assert Q(1) / lower_product < Q(4, 3)
    # Polynomial nonnegativity behind derivative <= 5/3 for every i >= 1:
    # 5(2i+1)^2 - 9(2i+3) = (i-1)(20i+22).
    for i in range(1, 129):
        assert 5 * (2 * i + 1) ** 2 - 9 * (2 * i + 3) == (i - 1) * (20 * i + 22)
    print("Exact telescoping and inverse-density estimates: N=1..128 OK")
    print("Scope: finite checks; analytic limits and functional extension require the subsequent proof")


def main():
    source_fidelity()
    prefix_checks()
    analytic_ingredients()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact sanity checks for Eldar's conjecture in OEIS A179875 (revision 13).

Direct gcd-filtered integer sums are compared with independently computed
Möbius signs and prime-factor formulas. This is a finite check, not a proof.
Source fidelity: A179875's displayed prefix and B(10) = B(11) = 7;
also Petridi arXiv:1612.07632v2, Section 1, first and second moments.
"""

from fractions import Fraction
from math import gcd, prod

LIMIT = 2000
SOURCE_PREFIX = [
    1, 6, 10, 22, 46, 58, 65, 69, 77, 82, 106, 129, 166, 178, 185,
    194, 210, 221, 226, 237, 254, 262, 265, 309, 321, 330, 346, 358,
    365, 382, 398, 417, 437, 454, 462, 466, 469, 473, 478, 482, 493,
    497, 502,
]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def factor(n):
    """Trial division; retain every prime exponent, including the n=1 case."""
    result = []
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            result.append((d, e))
        d += 1
    if n > 1:
        result.append((n, 1))
    return result


def mobius_sieve(limit):
    """Prime sieve, independently of trial factorization and residue sums."""
    mu = [1] * (limit + 1)
    mu[0] = 0
    composite = [False] * (limit + 1)
    for p in range(2, limit + 1):
        if composite[p]:
            continue
        for k in range(p, limit + 1, p):
            mu[k] = -mu[k]
            if k > p:
                composite[k] = True
        for k in range(p * p, limit + 1, p * p):
            mu[k] = 0
    return mu


def main():
    mu = mobius_sieve(LIMIT)
    means = {}
    errors = {}
    for n in range(1, LIMIT + 1):
        residues = [k for k in range(1, n + 1) if gcd(k, n) == 1]
        s1 = sum(residues)
        s2 = sum(k * k for k in residues)
        require(s1 > 0, f"denominator positivity n={n}")
        means[n] = Fraction(s2, s1)
        factors = factor(n)
        radical = prod(p for p, _ in factors)
        sign = (-1) ** len(factors)
        error = Fraction(sign * radical, n)
        errors[n] = error
        phi = prod((p - 1) * p ** (e - 1) for p, e in factors)
        require(phi == len(residues), f"totient fidelity n={n}")
        require(-1 <= error <= 1, f"sharp error bound n={n}")
        require((error == 1) == (mu[n] == 1), f"positive endpoint n={n}")
        require((error == -1) == (mu[n] == -1), f"negative endpoint n={n}")
        require(3 * means[n] == 2 * n + error, f"normalized moment n={n}")
        if n > 1:
            require(s1 == Fraction(n * phi, 2), f"Petridi first moment n={n}")
            require(s2 == Fraction(phi * (2 * n * n + sign * radical), 6),
                    f"Petridi second moment n={n}")
        else:
            require(residues == [1] and means[n] == 1, "OEIS n=1 convention")

    equal_indices = []
    for n in range(1, LIMIT):
        equal = means[n] == means[n + 1]
        predicted = mu[n] == 1 and mu[n + 1] == -1
        require(equal == predicted, f"Eldar biconditional n={n}")
        require(3 * (means[n + 1] - means[n]) ==
                2 + errors[n + 1] - errors[n], f"adjacent identity n={n}")
        require(means[n] <= means[n + 1], f"monotonicity n={n}")
        if equal:
            equal_indices.append(n)

    require([n for n in equal_indices if n <= SOURCE_PREFIX[-1]] == SOURCE_PREFIX,
            "OEIS A179875 displayed prefix")
    require(means[10] == means[11] == 7, "OEIS example B(10)=B(11)=7")
    require(means[1] == means[2] == 1, "first equal pair")
    require(means[6] == means[7] == Fraction(13, 3), "nonintegral equal pair")
    print(f"Fidelity: reproduced all {len(SOURCE_PREFIX)} displayed A179875 terms through 502")
    print("Fidelity: B(10)=B(11)=7; B(1)=B(2)=1; B(6)=B(7)=13/3")
    print(f"Fidelity: Petridi Section 1 moments verified for 2 <= n <= {LIMIT}")
    print(f"Exact normalized moments and sharp endpoint tests: {LIMIT} positive indices")
    print(f"Adjacent biconditionals and monotonicity: {LIMIT - 1} pairs; "
          f"{len(equal_indices)} equal pairs")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()

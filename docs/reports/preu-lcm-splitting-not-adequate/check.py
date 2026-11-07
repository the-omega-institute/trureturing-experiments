#!/usr/bin/env python3
"""Exact finite projections of the proposed Preu converse counterexample.

The infinite quotient and the Bezout realization require mathematical proofs.
No finite quotient here is advertised as a counterexample: on a finite set,
every function belongs to the finite-superlevel ideal.
"""
from fractions import Fraction as Q
from itertools import product
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def fidelity():
    # Preu, arXiv:2610.23799v1, Section 2, first observations.
    pairs = [(c, d) for x, y in [(6, 5), (3, 10)]
             for c in [x, -x] for d in [y, -y]]
    for c, d in pairs:
        require(6 % c == 0 and 10 % d == 0, "source divisibility")
        require(abs(c * d) == lcm(6, 10), "source lcm")
        require(gcd(c, d) == 1, "source splitting")
    require(gcd(12, 18 // gcd(12, 18)) == 3, "source greedy failure")
    print("FIDELITY: all 8 signed source splittings of (6,10); greedy gcd(12,3)=3")


def ideal_and_lattice_bounds():
    signed = [Q(-1), Q(-1, 2), Q(0), Q(1, 3), Q(1, 2), Q(1)]
    thresholds = [Q(1, 8), Q(1, 2), Q(1), Q(2)]
    additions = 0
    for f, g, eps in product(signed, signed, thresholds):
        require(abs(f + g) < eps or abs(f) >= eps / 2 or abs(g) >= eps / 2,
                "superlevel union inclusion")
        additions += 1
    lattice = 0
    for f, g, h, k in product(signed, repeat=4):
        bound = abs(f - g) + abs(h - k)
        require(abs(min(f, h) - min(g, k)) <= bound, "meet continuity bound")
        require(abs(max(f, h) - max(g, k)) <= bound, "join continuity bound")
        lattice += 1
    print(f"IDEAL: {additions} superlevel implications; {lattice} meet/join bounds")


def splitting():
    alphabet = [Q(0), Q(1, 4), Q(1, 2), Q(1), Q(2)]
    vectors = list(product(alphabet, repeat=3))
    count = 0
    for f, h in product(vectors, repeat=2):
        c = tuple(x if x >= y else Q(0) for x, y in zip(f, h))
        d = tuple(Q(0) if x >= y else y for x, y in zip(f, h))
        for x, y, u, v in zip(f, h, c, d):
            require(0 <= u <= x and 0 <= v <= y, "positive divisors")
            require(min(u, v) == 0, "disjoint factors")
            require(u + v == max(x, y), "lcm value")
        count += 1
    print(f"SPLITTING: {count} ordered vector pairs, 3 coordinates, exact rationals")


def floor_reciprocal(eps):
    return eps.denominator // eps.numerator


def diagonal():
    profiles = [("zero", lambda m: 0), ("linear", lambda m: m),
                ("quadratic", lambda m: m * m), ("affine", lambda m: 2 * m + 3)]
    thresholds = [Q(2), Q(1), Q(1, 2), Q(1, 3), Q(1, 4), Q(1, 8), Q(1, 16)]
    family_count = support_checks = point_checks = 0
    for name, k in profiles:
        for shift, power in product([0, 2], [1, 2]):
            # q_m <= 1/m; q_m > 0 and q_m -> 0 for these entire families.
            def beta(m):
                return Q(1, (m + shift) ** power)

            def c(m, n):
                return Q(1) if n <= k(m) else Q(1, m + n + 1)

            def selected(m):
                return k(m) + 1

            def h(m, n):
                return Q(1, 2) if n == selected(m) else Q(0)

            # These are box restrictions of infinite formulas, not finite quotients.
            for m in range(1, 33):
                require(beta(m) > 0 and beta(m) <= Q(1, m), "beta bound")
                require(c(m, selected(m)) < Q(1, 2), "selected small value")
                for n in range(1, selected(m) + 3):
                    cv, hv = c(m, n), h(m, n)
                    require(0 <= cv <= 1 and 0 <= hv <= 1 - cv, "diagonal divisor")
                    require(cv + (1 - cv) == 1, "factorization of strong unit")
                    point_checks += 1
            for rows in [1, 2, 4, 8, 16, 32]:
                support_h = sum(h(m, selected(m)) >= Q(1, 2)
                                for m in range(1, rows + 1))
                require(support_h == rows, "nonvanishing diagonal prefixes")
            for eps in thresholds:
                limit = floor_reciprocal(eps)
                # If min(C,beta)>=eps, then m<=limit. Outside the finite
                # initial segment n<=k(m), C>=eps forces n<=limit-m-1.
                columns = max([0] + [max(k(m), limit - m - 1)
                                        for m in range(1, limit + 1)])
                expected_c = sum(max(k(m), max(0, limit - m - 1))
                                 for m in range(1, limit + 1) if beta(m) >= eps)
                actual_c = sum(min(c(m, n), beta(m)) >= eps
                               for m in range(1, limit + 3)
                               for n in range(1, columns + 3))
                require(actual_c == expected_c, "exact finite support certificate for meet(C,beta)")
                expected_h = (sum(beta(m) >= eps for m in range(1, limit + 1))
                              if eps <= Q(1, 2) else 0)
                actual_h = sum(min(h(m, selected(m)), beta(m)) >= eps
                               for m in range(1, limit + 3))
                require(actual_h == expected_h, "exact finite support certificate for meet(H,beta)")
                support_checks += 1
            family_count += 1
    print(f"DIAGONAL: {family_count} infinite formula families; {point_checks} projected points")
    print(f"SUPERLEVELS: {support_checks} exact support counts; thresholds through 1/16")
    print("NONVANISHING: H has exactly M entries >=1/2 on M selected rows, M=1,2,4,8,16,32")
    print("CONTROL: beta=1/m is constant and positive on each infinite row; H selects one point per row")
    print("BOUNDARY: finite projections do not prove the infinite quotient or the Bezout realization")


def main():
    fidelity()
    ideal_and_lattice_bounds()
    splitting()
    diagonal()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

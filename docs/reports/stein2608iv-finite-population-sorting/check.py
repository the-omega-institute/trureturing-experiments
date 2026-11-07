#!/usr/bin/env python3
"""Exact finite-population SAP checks for Stein, arXiv:2608.21955v1.

Enumerate labelled permutations, including all their multiplicities when values
repeat. Integrate the independent uniform thresholds and fair singleton coin
exactly, rather than simulate them. E_A and E_B are ratios of expected assigned
weight to expected assigned count, as in the paper's fractional-marble definition.
Only Python's standard library is needed. No floating-point arithmetic is used.
"""
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations
from math import factorial


def enumerate_population(w):
    n = len(w)
    assert n > 0 and all(Q(0) < t < Q(1) for t in w)
    k, d = divmod(n, 2)
    m = sum(w, Q(0)) / n
    variance = sum(((t - m) ** 2 for t in w), Q(0)) / n
    ca = cb = wa = wb = Q(0)
    pair_x = [Q(0)] * k
    pair_u = [Q(0)] * k
    pair_xu = [Q(0)] * k
    last_weight = Q(0)
    permutation_count = 0
    for perm in permutations(range(n)):
        permutation_count += 1
        c_a = c_b = Q(d, 2)
        w_a = w_b = Q(0)
        for j in range(k):
            x, u = w[perm[2 * j]], w[perm[2 * j + 1]]
            pair_x[j] += x
            pair_u[j] += u
            pair_xu[j] += x * u
            c_a += x
            c_b += 1 - x
            w_a += x * u
            w_b += (1 - x) * u
        if d:
            singleton = w[perm[-1]]
            last_weight += singleton
            w_a += singleton / 2
            w_b += singleton / 2
        ca += c_a
        cb += c_b
        wa += w_a
        wb += w_b
    assert permutation_count == factorial(n)
    ca, cb, wa, wb = (t / permutation_count for t in (ca, cb, wa, wb))
    assert ca > 0 and cb > 0
    assert ca == k * m + Q(d, 2)
    assert cb == k * (1 - m) + Q(d, 2)
    ea, eb = wa / ca, wb / cb
    if n == 1:
        assert ea == m == eb == w[0]
        assert wa == wb == w[0] / 2
    else:
        cross = m * m - variance / (n - 1)
        for j in range(k):
            assert pair_x[j] / permutation_count == m
            assert pair_u[j] / permutation_count == m
            assert pair_xu[j] / permutation_count == cross
        assert wa == k * cross + Q(d, 2) * m
        assert wb == k * (m - cross) + Q(d, 2) * m
        assert ea == m - k * variance / ((n - 1) * (k * m + Q(d, 2)))
        assert eb == m + k * variance / ((n - 1) * (k * (1 - m) + Q(d, 2)))
    if d:
        assert last_weight / permutation_count == m
    assert ea <= m <= eb
    assert (ea == eb) == (variance == 0)
    if variance > 0:
        assert ea < m < eb
    return ea, m, eb, permutation_count


def fidelity():
    # Source Corollary 1.1: exact formulas, not its noisy simulation table.
    s, ell = Q(1, 4), Q(3, 4)
    ea, m, eb, _ = enumerate_population((s, ell))
    assert ea == 2 * s * ell / (s + ell) == Q(3, 8)
    assert eb == (s + ell - 2 * s * ell) / (2 - s - ell) == Q(5, 8)
    print(f"FIDELITY Corollary 1.1: S={s} L={ell} EA={ea} m={m} EB={eb}")
    # A stated numerical source datum, Section IV's final table row (page 20).
    # Its displayed empirical EA/EB come from 1,000 random trials and are not
    # exact theoretical expectations; only its exact population mean is used.
    w = (Q(8, 100),) * 2 + (Q(97, 100),) * 2
    ea, m, eb, _ = enumerate_population(w)
    assert m == Q(525, 1000) == Q(21, 40)
    print(f"FIDELITY Section IV final table: S=2/25 #S=2 L=97/100 #L=2 m={m}")
    print(f"EXACT expectation for that population: EA={ea} EB={eb}")


def main():
    fidelity()
    population_count = permutation_count = odd_count = constant_count = 0
    grid = (Q(1, 4), Q(1, 2), Q(3, 4))
    for n in range(1, 8):
        cases = perms = 0
        for w in combinations_with_replacement(grid, n):
            *_, checked = enumerate_population(w)
            cases += 1
            perms += checked
            odd_count += n % 2
            constant_count += len(set(w)) == 1
        population_count += cases
        permutation_count += perms
        print(f"GRID N={n} populations={cases} labelled_permutations={perms}")
    # Separate two-species multiplicity family and its variance identity.
    family_cases = family_perms = 0
    s, ell = Q(1, 5), Q(4, 5)
    for n in range(1, 8):
        for a in range(n + 1):
            b = n - a
            w = (s,) * a + (ell,) * b
            mean = sum(w, Q(0)) / n
            var = sum(((t - mean) ** 2 for t in w), Q(0)) / n
            assert var == Q(a * b, n * n) * (ell - s) ** 2
            *_, checked = enumerate_population(w)
            family_cases += 1
            family_perms += checked
    print(f"GRID_TOTAL populations={population_count} labelled_permutations={permutation_count} "
          f"odd_populations={odd_count} constant_populations={constant_count}")
    print(f"TWO_SPECIES_TOTAL populations={family_cases} labelled_permutations={family_perms}")
    print("Verified covariance count, expected counts/weights, both closed forms, singleton rule, "
          "weak inequalities, strictness for nonconstant populations and constant equality.")
    print("Finite checks support the proposed general argument; no Lean proof is claimed.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

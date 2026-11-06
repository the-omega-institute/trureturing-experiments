#!/usr/bin/env python3
"""Exact finite checks for Sedlak et al., §V.C.

All values are fractions.  For diagonal testers the covering SDP reduces to the
finite LP whose exact optimum is sum_x max(rho_x, sigma_x); this script checks
that identity over a small exhaustive rational family and checks the paper's
orthogonal qubit example.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import product


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def split_effects(masses: tuple[int, ...], outcomes: int):
    """All outcome arrays with integer numerator masses per input coordinate."""
    per_coordinate = [list(compositions(m, outcomes)) for m in masses]
    for columns in product(*per_coordinate):
        yield tuple(tuple(columns[x][j] for x in range(len(masses)))
                    for j in range(outcomes))


def covering_optimum(rho: tuple[Fraction, ...], sigma: tuple[Fraction, ...]) -> Fraction:
    # For a diagonal tester, every output diagonal entry has row total rho[x]
    # and column total sigma[x].  A nonnegative joint table needs, and can
    # attain, total max(rho[x], sigma[x]); the common H uses that total.
    return sum(max(r, s) for r, s in zip(rho, sigma))


def robustness_from_cover(t: Fraction) -> Fraction:
    return 1 - 1 / t


def products_zero(A, B) -> bool:
    return all(a_i * b_i == 0 for a in A for b in B for a_i, b_i in zip(a, b))


def check_source_fidelity() -> None:
    # T_V and T_H from the source: output/input basis order (V,H) x (V,H).
    A = (
        (Fraction(1), Fraction(0), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1), Fraction(0)),
    )
    B = (
        (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
        (Fraction(0), Fraction(1), Fraction(0), Fraction(0)),
    )
    rho = (Fraction(1), Fraction(0))
    sigma = (Fraction(0), Fraction(1))
    assert all(sum(row[out * 2 + x] for row in A) == rho[x]
               for out in range(2) for x in range(2))
    assert all(sum(row[out * 2 + x] for row in B) == sigma[x]
               for out in range(2) for x in range(2))
    assert products_zero(A, B)
    t = covering_optimum(rho, sigma)
    assert t == 2 and robustness_from_cover(t) == Fraction(1, 2)
    print("FIDELITY_OK source qubit orthogonal testers: t*=2, R_t=1/2")


def check_small_diagonal_family() -> None:
    checked = 0
    for d0, d1, p, q, denominator in product(
        (1, 2), (1, 2), (1, 2), (1, 2), (2, 4)
    ):
        states = list(compositions(denominator, d0))
        for r_num in states:
            for s_num in states:
                rho = tuple(Fraction(n, denominator) for n in r_num)
                sigma = tuple(Fraction(n, denominator) for n in s_num)
                for A_num in split_effects(r_num, p):
                    for B_num in split_effects(s_num, q):
                        A = tuple(tuple(Fraction(n, denominator) for n in row)
                                  for row in A_num)
                        B = tuple(tuple(Fraction(n, denominator) for n in row)
                                  for row in B_num)
                        t = covering_optimum(rho, sigma)
                        R = robustness_from_cover(t)
                        orth_states = all(r * s == 0 for r, s in zip(rho, sigma))
                        orth_effects = products_zero(A, B)
                        assert (R == Fraction(1, 2)) == orth_states
                        assert orth_states == orth_effects
                        checked += 1
    print(f"SMALL_DIAGONAL_OK cases={checked} exact rational covering identity")


if __name__ == "__main__":
    check_source_fidelity()
    check_small_diagonal_family()
    print("ALL_CHECKS_OK")

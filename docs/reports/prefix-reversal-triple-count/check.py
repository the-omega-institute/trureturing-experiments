#!/usr/bin/env python3
"""Exact finite checks for Blanco et al., arXiv:2511.16959v3, Conjecture 1.

Python 3; SymPy 1.14.0. No numerical approximation or randomized group test.
All pairs are checked for every 3 <= n <= 24, without obstruction pruning.
The finite error bound is fixed at C=2, N=6 before running this program;
it is a finite consistency test, not a verification of an asymptotic claim.
"""

from collections import deque
from fractions import Fraction
from math import factorial

import sympy
from sympy.combinatorics import Permutation, PermutationGroup


MAX_N = 24
COEFFICIENTS = (Fraction(1, 13), Fraction(1, 9), Fraction(1, 10), Fraction(1, 7))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def reversal(n, length):
    """Zero-based conjugate of the source's one-based r_length."""
    require(2 <= length <= n, "invalid source reversal length")
    return tuple(length - 1 - i if i < length else i for i in range(n))


def compose(left, right):
    """Right-to-left composition, as in Section 2 of the source."""
    return tuple(left[right[i]] for i in range(len(left)))


def bfs_order(generators):
    """Independent whole-subgroup enumeration; used only for n <= 7."""
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    pending = deque([identity])
    while pending:
        current = pending.popleft()
        for generator in generators:
            candidate = compose(generator, current)
            if candidate not in seen:
                seen.add(candidate)
                pending.append(candidate)
    return len(seen)


def main():
    print("Python 3, SymPy " + sympy.__version__)
    print("Exact deterministic Schreier-Sims orders; tuple-BFS cross-check for n <= 7.")
    print("Scope: every pair 2 <= k < m < n, for every 3 <= n <= 24.")
    orders = {}
    generating = {}
    total_pairs = 0
    bfs_pairs = 0
    involutions = 0
    for n in range(3, MAX_N + 1):
        reversals = {length: reversal(n, length) for length in range(2, n + 1)}
        for r in reversals.values():
            require(compose(r, r) == tuple(range(n)), "reversal is not involutive")
            involutions += 1
        generating[n] = set()
        for m in range(3, n):
            for k in range(2, m):
                generators = [reversals[length] for length in (n, m, k)]
                group = PermutationGroup([Permutation(list(g), size=n) for g in generators])
                order = int(group.order())
                require(factorial(n) % order == 0, "subgroup order does not divide n!")
                orders[n, m, k] = order
                total_pairs += 1
                if n <= 7:
                    require(order == bfs_order(generators), "BFS/Schreier-Sims mismatch")
                    bfs_pairs += 1
                if order == factorial(n):
                    generating[n].add((m, k))

    # Source Section 4.2 gives these eight triples, not just their count.
    expected_seven = {(4, 2), (5, 2), (6, 2), (4, 3), (6, 3), (5, 4), (6, 4), (6, 5)}
    require(generating[7] == expected_seven, "source n=7 generating-triple list mismatch")
    print("FIDELITY: source Section 4.2 eight n=7 triples reproduced exactly:")
    print(sorted((7, m, k) for m, k in generating[7]))

    # Source Theorem 4 is a second, uniform-in-k finite fidelity check.
    boundary_cases = 0
    for n in range(5, MAX_N + 1):
        for k in range(2, n - 2):
            require(((n - 2, k) in generating[n]) == ((n - k) % 2 == 1),
                    "source Theorem 4 mismatch")
            boundary_cases += 1
    print("FIDELITY: Theorem 4, m=n-2 and opposite parity:", boundary_cases, "cases")

    triangle_cases = 0
    shift_four_cases = 0
    shift_two_cases = 0
    for (n, m, k), order in orders.items():
        if m + k < (n if n % 2 == 0 else n - 1):
            require(order < factorial(n), "triangular non-generation mismatch")
            triangle_cases += 1
        if order == factorial(n) and n + 4 <= MAX_N:
            require((m + 4, k + 4) in generating[n + 4], "shift-by-four mismatch")
            shift_four_cases += 1
        if order == factorial(n) and n % 4 in (0, 1) and n + 2 <= MAX_N:
            require((m + 2, k + 2) in generating[n + 2], "shift-by-two mismatch")
            shift_two_cases += 1
    print("STRUCTURE: parity-correct triangular exclusions:", triangle_cases)
    print("STRUCTURE: forward shifts by four:", shift_four_cases,
          "; parity-restricted shifts by two:", shift_two_cases)

    # Test exactly the proposed leading constants on the entire finite range.
    # Passing cannot establish existence of C_q,N_q valid for unbounded n.
    print("n,residue,f(n),c_q,error,error/n")
    errors = {q: [] for q in range(4)}
    for n in range(3, MAX_N + 1):
        count = len(generating[n])
        q = n % 4
        error = abs(Fraction(count) - COEFFICIENTS[q] * n * n)
        print(f"{n},{q},{count},{COEFFICIENTS[q]},{error},{error / n}")
        if n >= 6:
            require(error <= 2 * n, "preregistered finite C=2,N=6 consistency check failed")
            errors[q].append(error / n)
    for q in range(4):
        print(f"FINITE_ERROR q={q}: max |f(n)-c_q*n^2|/n = {max(errors[q])}")
    print("CHECKED_PAIRS", total_pairs, "BFS_CROSS_CHECKS", bfs_pairs,
          "INVOLUTION_CHECKS", involutions)
    print("LIMIT: finite checks establish neither the conjecture nor its negation.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

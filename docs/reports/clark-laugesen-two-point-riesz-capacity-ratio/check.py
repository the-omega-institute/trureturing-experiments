#!/usr/bin/env python3
"""Exact finite checks for Clark–Laugesen's two-point capacity-ratio conjecture.

Only integers and fractions.Fraction are used. No numerical optimization or
floating-point comparison is used. This is finite evidence, not a proof of the
compact-set conjecture.
"""

from fractions import Fraction as Q
from itertools import combinations, product


def solve_unique(matrix, rhs):
    """Gauss–Jordan solve over Q; return None for a singular system."""
    n = len(rhs)
    aug = [[Q(x) for x in row] + [Q(b)] for row, b in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if aug[i][col]), None)
        if pivot is None:
            return None
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [v / scale for v in aug[col]]
        for i in range(n):
            if i != col and aug[i][col]:
                scale = aug[i][col]
                aug[i] = [v - scale * w for v, w in zip(aug[i], aug[col])]
    return [aug[i][-1] for i in range(n)]


def quadratic(matrix, weights):
    return sum((weights[i] * weights[j] * matrix[i][j]
                for i in range(len(matrix)) for j in range(len(matrix))), Q(0))


def simplex_max(matrix):
    """Exact maximum of w^T A w on the closed probability simplex.

    Enumerate all nonempty supports, solving A_S w = c*1, sum(w)=1.
    A maximum of smallest support has an invertible bordered system: a null
    direction would have sum zero and A_S v=0, so one can move to the boundary
    with the same objective. Thus singular systems can be skipped. All retained
    candidates are actual feasible measures, including boundary zero masses.
    """
    n = len(matrix)
    best = Q(-1)
    best_weights = None
    for size in range(1, n + 1):
        for support in combinations(range(n), size):
            bordered = [[matrix[i][j] for j in support] + [-1] for i in support]
            bordered.append([1] * size + [0])
            solution = solve_unique(bordered, [0] * size + [1])
            if solution is None or any(w < 0 for w in solution[:-1]):
                continue
            weights = [Q(0)] * n
            for i, w in zip(support, solution[:-1]):
                weights[i] = w
            assert sum(weights) == 1
            value = quadratic(matrix, weights)
            assert value == solution[-1]
            if value > best:
                best, best_weights = value, weights
    assert best_weights is not None
    return best, best_weights


def matrix_from_edges(n, edges):
    matrix = [[0] * n for _ in range(n)]
    for (i, j), value in zip(combinations(range(n), 2), edges):
        matrix[i][j] = matrix[j][i] = value
    return matrix


def main():
    # Source fidelity: Clark–Laugesen introduction, regular k-point formula:
    # Cap_p = ((k-1)/k)^(-1/p) * diameter. Equivalently U_t=(k-1)/k*D^t.
    # D=base^4 and t=a/4 make D^t=base^a exactly integral.
    fidelity = 0
    for n, base, a in product(range(2, 7), (1, 2), range(1, 9)):
        matrix = [[base ** a if i != j else 0 for j in range(n)] for i in range(n)]
        maximum, _ = simplex_max(matrix)
        assert maximum == Q(n - 1, n) * base ** a
        assert quadratic(matrix, [Q(1, n)] * n) == maximum
        fidelity += 1
    print(f"FIDELITY regular-simplex energy formula: {fidelity} exact cases")

    # Independent check of the source's explicit two-point calculation,
    # including lambda=0,1 (zero masses are allowed).
    mass_checks = 0
    for denominator in range(1, 17):
        for numerator in range(denominator + 1):
            mass = Q(numerator, denominator)
            energy = 2 * mass * (1 - mass)
            assert Q(1, 2) - energy == 2 * (mass - Q(1, 2)) ** 2
            assert energy <= Q(1, 2)
            mass_checks += 1
    print(f"FIDELITY two-point mass formula: {mass_checks} rational masses")

    # Complete enumeration of symmetric zero-diagonal kernels D_ij in {0,1,16}
    # on 2,3,4 labels. These include nonmetric kernels. For t=a/4, their powered
    # entries are exactly base_ij^a with base_ij in {0,1,2}.
    exponents = (1, 2, 3, 4, 5, 6, 7, 8, 12)
    kernels = energy_checks = ratio_checks = in_open_range = 0
    boundary_optima = 0
    for n in range(2, 5):
        for edges in product((0, 1, 2), repeat=n * (n - 1) // 2):
            largest = max(edges)
            if largest == 0:
                continue
            kernels += 1
            energies = {}
            for a in exponents:
                powered = matrix_from_edges(n, [base ** a for base in edges])
                maximum, weights = simplex_max(powered)
                assert maximum > 0
                assert Q(largest ** a, 2) <= maximum <= largest ** a
                assert quadratic(powered, weights) == maximum
                boundary_optima += int(any(w == 0 for w in weights))
                energies[a] = maximum
                energy_checks += 1
            for a, b in combinations(exponents, 2):
                ur, us = energies[a], energies[b]
                # Check optimized-kernel comparison separately from the ratio.
                assert us <= largest ** (b - a) * ur
                # Raise Cs/Cr <= 2^(1/r-1/s) to positive power a*b/4:
                # Us^a <= 2^(b-a) Ur^b. No roots or approximation are needed.
                assert us ** a <= 2 ** (b - a) * ur ** b
                if n == 2:
                    assert us ** a == 2 ** (b - a) * ur ** b
                ratio_checks += 1
                in_open_range += int(b < 8)
    assert kernels == sum(3 ** (n * (n - 1) // 2) - 1 for n in range(2, 5))
    print(f"KERNELS exhaustive nonzero symmetric kernels: {kernels}")
    print(f"OPTIMIZATION exact simplex maxima: {energy_checks}")
    print(f"OPTIMIZATION returned measures with zero masses: {boundary_optima}")
    print(f"RATIO exact exponent-pair comparisons: {ratio_checks}")
    print(f"RATIO comparisons with 0 < r < s < 2: {in_open_range}")
    print("SCOPE finite kernels only; compact-set and arbitrary-exponent proof remains unformalized")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

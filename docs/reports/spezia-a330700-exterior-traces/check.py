#!/usr/bin/env python3
"""Exact finite checks for Spezia's OEIS A330700, second conjecture line.

Run with Python 3 and SymPy. All mathematical computations use integers or
fractions.Fraction. This program supplies finite evidence, not a uniform proof.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb
import sys

import sympy
from sympy import Matrix, Symbol, zeros


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def entry(i, j):
    """Literal A330287 formula, with positive one-based indices."""
    return i * j if i < 3 or j < 3 else 2 * (i + j) - 4


def matrix(n):
    return [[entry(i, j) for j in range(1, n + 1)]
            for i in range(1, n + 1)]


def determinant(rows):
    """Independent exact Gaussian determinant, including the empty matrix."""
    n = len(rows)
    require(all(len(row) == n for row in rows), "square determinant input")
    a = [[Fraction(x) for x in row] for row in rows]
    answer = Fraction(1)
    for col in range(n):
        pivot = next((r for r in range(col, n) if a[r][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[pivot], a[col] = a[col], a[pivot]
            answer = -answer
        pivot_value = a[col][col]
        answer *= pivot_value
        for r in range(col + 1, n):
            ratio = a[r][col] / pivot_value
            for c in range(col + 1, n):
                a[r][c] -= ratio * a[col][c]
            a[r][col] = Fraction(0)
    return answer


def minor(a, rows, columns):
    return determinant([[a[i][j] for j in columns] for i in rows])


def permanent(a):
    result = 0
    for p in permutations(range(len(a))):
        term = 1
        for i, j in enumerate(p):
            term *= a[i][j]
        result += term
    return result


def main():
    print("Python", sys.version.split()[0], "SymPy", sympy.__version__)
    print("Source: A330700 revision 23; A330287 revision 31; A319840 revision 105")
    require(matrix(1) == [[1]], "A330287 M(1)")
    require(matrix(2) == [[1, 2], [2, 4]], "A330287 M(2)")
    require(matrix(3) == [[1, 2, 3], [2, 4, 6], [3, 6, 8]], "A330287 M(3)")
    source_permanents = [1, 1, 8, 208, 11488, 1093056, 158972160, 32734095360]
    computed_permanents = [permanent(matrix(n)) for n in range(8)]
    require(computed_permanents == source_permanents, "A330287 permanents n=0..7")
    print("FIDELITY M(1), M(2), M(3): matched verbatim examples")
    print("FIDELITY A330287 permanents n=0..7:", computed_permanents)

    perimeter_cases = 0
    for i in range(1, 41):
        for j in range(1, 41):
            boundary = sum(r in (0, i - 1) or c in (0, j - 1)
                           for r in range(i) for c in range(j))
            require(entry(i, j) == boundary, f"A319840 boundary {i},{j}")
            perimeter_cases += 1
    print("FIDELITY geometric perimeter counts, 1<=i,j<=40:", perimeter_cases)

    x = Symbol("x")
    coefficient_checks = 0
    factorization_entries = 0
    zero_space_checks = 0
    sharp_rank_checks = 0
    second_traces = {}
    principal_checks = 0
    for n in range(41):
        a = matrix(n)
        v = list(range(1, n + 1))
        w = [max(i - 2, 0) for i in v]
        for i in range(n):
            for j in range(n):
                require(a[i][j] == v[i] * v[j] - w[i] * w[j],
                        f"factorization n={n},i={i+1},j={j+1}")
                factorization_entries += 1
        dense = Matrix(a) if n else zeros(0, 0)
        polynomial = dense.charpoly(x).as_poly()
        for k in range(3, n + 1):
            require(polynomial.nth(n - k) == 0, f"charpoly n={n},k={k}")
            coefficient_checks += 1
        for k in range(max(3, n + 1), 43):
            require(comb(n, k) == 0, f"zero exterior space n={n},k={k}")
            zero_space_checks += 1
        if n >= 2:
            second_traces[n] = int(polynomial.nth(n - 2))
        if n >= 3:
            require(minor(a, (0, 2), (0, 2)) == -1, f"rank witness n={n}")
            sharp_rank_checks += 1
        if n <= 10:
            for k in range(3, n + 1):
                total = Fraction(0)
                for s in combinations(range(n), k):
                    d = minor(a, s, s)
                    require(d == 0, f"principal minor n={n},k={k},S={s}")
                    total += d
                    principal_checks += 1
                require(total == (-1) ** k * polynomial.nth(n - k),
                        f"principal-trace/charpoly crosscheck n={n},k={k}")
    print("SETTLEMENT factorization entries, 0<=n<=40:", factorization_entries)
    print("SETTLEMENT higher characteristic coefficients, 3<=k<=n<=40:", coefficient_checks)
    print("SETTLEMENT zero-space cases, 0<=n<=40, 3<=k<=42, k>n:", zero_space_checks)
    print("SETTLEMENT all principal minors of order >=3, 0<=n<=10:", principal_checks)
    print("NONDEGENERACY principal minor {1,3}=-1, 3<=n<=40:", sharp_rank_checks)

    compound_entries = 0
    for n in range(3, 8):
        a = matrix(n)
        for k in range(3, n + 1):
            subsets = list(combinations(range(n), k))
            for rows in subsets:
                for columns in subsets:
                    require(minor(a, rows, columns) == 0,
                            f"compound entry n={n},k={k},rows={rows},cols={columns}")
                    compound_entries += 1
    print("SETTLEMENT every higher compound-matrix entry, 3<=n<=7:", compound_entries)
    require(second_traces[3] == -5, "separate first conjecture diagnostic")
    a3 = (3 - 1) * 3 * (2 * 3 ** 2 + 4 * 3 - 1) // 6
    require(a3 == 29 and abs(second_traces[3]) != a3,
            "first conjecture's source indexing mismatch")
    print("SCOPE separate first conjecture: traceExt_2(M(3))=-5, a(3)=29; not nominated")
    print("Finite exact checks only; no Lean proof or unbounded verification.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

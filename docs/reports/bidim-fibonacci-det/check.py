#!/usr/bin/env python3
"""Exact finite checks for Costa--Catarino's bidimensional Fibonacci matrix."""

LIMIT = 100


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def add(z, w):
    return z[0] + w[0], z[1] + w[1]


def neg(z):
    return -z[0], -z[1]


def sub(z, w):
    return add(z, neg(w))


def mul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def det(matrix):
    return sub(mul(matrix[0][0], matrix[1][1]),
               mul(matrix[0][1], matrix[1][0]))


def main():
    # Construct the source sequence from Definition 2.1, independently of (2.6).
    side = LIMIT + 3
    grid = [[(0, 0) for _ in range(side)] for _ in range(side)]
    grid[0][0], grid[1][0] = (0, 0), (1, 0)
    grid[0][1], grid[1][1] = (0, 1), (1, 1)
    for m in range(2):
        for n in range(2, side):
            grid[m][n] = add(grid[m][n - 1], grid[m][n - 2])
    for m in range(2, side):
        for n in range(side):
            grid[m][n] = add(grid[m - 1][n], grid[m - 2][n])

    fib = [0, 1]
    for _ in range(side + 2):
        fib.append(fib[-1] + fib[-2])
    for m in range(side):
        for n in range(side):
            require(grid[m][n] == (fib[m] * fib[n + 1], fib[m + 1] * fib[n]),
                    f"Lemma 2.5 failed at {(m, n)}")
            if n >= 2:
                require(grid[m][n] == add(grid[m][n - 1], grid[m][n - 2]),
                        f"second-coordinate recurrence failed at {(m, n)}")
    print(f"recurrences_and_lemma_2_5: {side * side} Gaussian-integer entries OK")

    def matrix(m, n):
        return ((grid[m + 2][n + 2], neg(grid[m + 1][n + 1])),
                (grid[m + 1][n + 1], neg(grid[m][n])))

    examples = {(0, 0): (0, 2), (1, 0): (-7, 1),
                (0, 1): (7, 1), (1, 1): (0, -4)}
    for indices, expected in examples.items():
        actual = det(matrix(*indices))
        require(actual == expected, f"Example 4.1 failed at {indices}")
        print(f"fidelity_Example_4_1: det{indices} = {actual} (real, imaginary) OK")
    require(matrix(1, 0) == (((4, 3), (-1, -2)), ((1, 2), (-1, 0))),
            "matrix (3.3) failed")

    t = [fib[r + 1] * fib[r + 2] for r in range(LIMIT + 1)]
    for r in range(LIMIT):
        p = fib[r] * fib[r + 3]
        require(t[r] - p == (-1) ** r, f"source Lemma 1.3 failed at {r}")
        require(t[r] > 0, f"T positivity failed at {r}")
        require(t[r + 1] - t[r] == fib[r + 2] ** 2 > 0,
                f"T strict increase failed at {r}")
    print(f"scalar_identities_and_strict_increase: 0 <= r < {LIMIT} OK")

    minimum = None
    for m in range(LIMIT):
        for n in range(LIMIT):
            value = det(matrix(m, n))
            expanded = sub(mul(grid[m + 1][n + 1], grid[m + 1][n + 1]),
                           mul(grid[m + 2][n + 2], grid[m][n]))
            require(value == expanded, f"determinant expansion failed at {(m, n)}")
            predicted_imaginary = (-1) ** m * t[n] + (-1) ** n * t[m]
            require(value[1] == predicted_imaginary,
                    f"mixed-index imaginary identity failed at {(m, n)}")
            require(value[1] != 0 and value != (0, 0),
                    f"nonvanishing failed at {(m, n)}")
            transposed_indices = det(matrix(n, m))
            require(transposed_indices == (-value[0], value[1]),
                    f"source Proposition 4.5 failed at {(m, n)}")
            magnitude = abs(value[1])
            minimum = magnitude if minimum is None else min(minimum, magnitude)
    print(f"mixed_identity_and_nonvanishing: {LIMIT ** 2} pairs, 0 <= m,n < {LIMIT} OK")
    print(f"boundary_pairs_included: {2 * LIMIT - 1}; minimum_abs_imaginary: {minimum}")
    for n in range(LIMIT):
        require(det(matrix(n, n)) == (0, 2 * (-1) ** n * t[n]),
                f"source Theorem 4.3 formula failed at {n}")
    print(f"fidelity_Theorem_4_3: {LIMIT} diagonal determinants OK")

    # Record an actual printed discrepancy, never use that formula as a premise.
    m = 1
    printed_axis_formula = (-fib[m] * (fib[m + 2] + fib[m + 4]), 2 * (-1) ** m)
    actual = det(matrix(m, 0))
    require(printed_axis_formula == (-7, -2) and actual == (-7, 1),
            "source formula discrepancy was not reproduced")
    print(f"SOURCE_DISCREPANCY_Proposition_4_6: printed={printed_axis_formula}, "
          f"Definition_3_1_and_Example_4_1={actual}")
    print("scope: finite exact checks only; no proof of the universal conjecture")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

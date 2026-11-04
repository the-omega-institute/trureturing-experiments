#!/usr/bin/env python3
"""Exact checks for the four-sided Fibonacci hitting-time conjecture."""

from fractions import Fraction

from sympy import Rational, sqrt


def fibonacci(n: int) -> list[int]:
    values = [0, 1]
    while len(values) <= n:
        values.append(values[-1] + values[-2])
    return values


def renewal_masses(limit: int) -> list[Fraction]:
    masses = [Fraction(1)]
    for n in range(1, limit + 1):
        masses.append(
            sum((masses[n - step] if n >= step else Fraction(0)) for step in range(1, 5))
            / 4
        )
    return masses


def first_passage_without_target(threshold: int, targets: set[int]) -> Fraction:
    """Probability of crossing threshold before a positive Fibonacci hit."""
    prefixes: dict[tuple[int, int], int] = {(0, 0): 1}
    result = Fraction(0)
    for length in range(threshold + 1):
        current = [
            (sum_so_far, prefix_length, count)
            for (sum_so_far, prefix_length), count in prefixes.items()
            if prefix_length == length
        ]
        for sum_so_far, prefix_length, count in current:
            for roll in range(1, 5):
                next_sum = sum_so_far + roll
                if next_sum in targets:
                    continue
                if next_sum >= threshold:
                    result += Fraction(count, 4 ** (prefix_length + 1))
                else:
                    key = (next_sum, prefix_length + 1)
                    prefixes[key] = prefixes.get(key, 0) + count
    return result


def main() -> None:
    # Fidelity: the paper’s four-sided heuristic is (3/5)φ < 1, checked in
    # exact quadratic arithmetic, with φ = (1 + √5)/2.
    ratio = Rational(3, 5) * (1 + sqrt(5)) / 2
    assert ratio < 1
    print(f"source_ratio_exact={ratio}")

    # Exact renewal recurrence for the unrestricted {1,2,3,4}-walk.
    masses = renewal_masses(256)
    expected_initial = {
        6: Fraction(1589, 4096),
        7: Fraction(6665, 16384),
        8: Fraction(26925, 65536),
        9: Fraction(102625, 262144),
    }
    for index, value in expected_initial.items():
        assert masses[index] == value
    assert all(value >= Fraction(5, 13) for value in masses[6:])
    print("renewal_minorant_exact=verified_through_n=256")

    # Exact Fibonacci potential and geometric survival checks from the proof.
    fib = fibonacci(15)
    q = Fraction(8, 13)

    def potential(k: int) -> int:
        return 169 * fib[k - 1] + 104 * fib[k - 2]

    for k in range(8, 14):
        assert potential(k) == fib[k - 1] + q * potential(k + 1)
    for k in range(8, 15):
        partial = sum(
            (q ** (j - 8)) * fib[j - 1] for j in range(8, k)
        )
        assert partial <= 3029
    print("fibonacci_potential_exact=verified_k=8..14")

    # Independent finite-prefix enumeration checks the conditional survival
    # contraction on the same unrestricted roll history.
    targets = set(fib[1:])
    passage = {k: first_passage_without_target(fib[k], targets) for k in range(8, 13)}
    for k in range(8, 12):
        assert passage[k + 1] <= q * passage[k]
    print("conditional_survival_exact=verified_k=8..11")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

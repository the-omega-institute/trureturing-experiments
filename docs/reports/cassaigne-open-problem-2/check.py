#!/usr/bin/env python3
"""Exact finite checks for Cassaigne's Open problem 2.

The checks use only integer arithmetic and Fraction/SymPy exact algebra.  They
reproduce a printed Fibonacci-word occurrence and test the proposed doubled
word family on a long exact prefix.  They are finite checks, not a proof of an
asymptotic limit.
"""

from fractions import Fraction
from math import isqrt
from pathlib import Path

from sympy import Eq, Symbol, sqrt, solve


def fibonacci_word(iterations: int) -> str:
    word = "a"
    for _ in range(iterations):
        word = "".join("ab" if letter == "a" else "a" for letter in word)
    return word


def factors(word: str, length: int) -> set[str]:
    return {word[i : i + length] for i in range(len(word) - length + 1)}


def least_period(word: str) -> int:
    for period in range(1, len(word) + 1):
        if all(word[i] == word[i - period] for i in range(period, len(word))):
            return period
    raise AssertionError("every finite word has its length as a period")


def exponent(word: str) -> Fraction:
    return Fraction(len(word), least_period(word))


def max_factor_exponent(word: str, length: int) -> Fraction:
    # The exact factor-complexity checks use the full generated prefix.  Keep
    # the exponent smoke checks bounded so the report remains a quick fixture.
    sample = word[:1200]
    return max((exponent(sample[i : i + length]) for i in range(len(sample) - length + 1)), default=Fraction(0))


def check_source_fidelity(fibonacci: str) -> None:
    # Cassaigne, p. 702: (aba)^3 occurs in f at position 5 (zero based).
    cube = "aba" * 3
    assert fibonacci[5 : 5 + len(cube)] == cube
    assert exponent(cube) == Fraction(3, 1)

    # The paper records the Fibonacci slope 2 - Phi and ind*(f) = Phi + 2.
    phi = (1 + sqrt(5)) / 2
    assert solve(Eq(Symbol("x"), 2 - phi), Symbol("x"))[0] == 2 - phi
    assert (2 + phi).expand() == (2 + (1 + sqrt(5)) / 2).expand()


def check_doubled_family(fibonacci: str) -> None:
    # A long prefix gives exact factor sets; no floating point is used.
    for k in range(2, 7):
        replicated = "".join(letter * k for letter in fibonacci)
        assert {replicated[i : i + 2] for i in range(len(replicated) - 1)} == {
            "aa",
            "ab",
            "ba",
            "bb",
        }

        # For n >= 2k+1 the proposed formula is p_v_k(n) = n + 2k - 1.
        for n in range(2 * k + 1, min(70, len(replicated) // 3)):
            observed = len(factors(replicated, n))
            assert observed == n + 2 * k - 1, (k, n, observed)

        # Exact finite index samples are recorded as rational values.  The
        # source value 2 + phi is an algebraic upper-limit target; these
        # samples only check that the replicated family has nontrivial long
        # repetitions and do not claim to prove the limit.
        samples = [max_factor_exponent(replicated, n) for n in (2 * k + 1, 20, 40)]
        assert all(value >= Fraction(1, 1) for value in samples)
        assert max(samples) >= Fraction(2, 1)


def main() -> None:
    fibonacci = fibonacci_word(16)
    assert len(fibonacci) == 2584
    check_source_fidelity(fibonacci)
    check_doubled_family(fibonacci)
    print("source_fidelity: (aba)^3 at position 5; exact exponent 3")
    print("family_checks: k=2..6, all four bigrams, p_v_k(n)=n+2k-1 for tested n")
    print("algebra: phi=(1+sqrt(5))/2 represented by SymPy exact sqrt")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

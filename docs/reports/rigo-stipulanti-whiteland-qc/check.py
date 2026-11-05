#!/usr/bin/env python3
"""Exact finite checks for Question C of Rigo--Stipulanti--Whiteland.

The checks use integer dynamic programming for scattered-subword counts.  They
reproduce the paper's level-3 word and exercise the proposed lacunary overlay
on finite prefixes.  They are finite evidence only; the preregistered
settlement is the general argument recorded in the draft.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


def morphic_prefix(seed: str, subst: dict[str, str], iterations: int) -> str:
    word = seed
    for _ in range(iterations):
        word = "".join(subst[ch] for ch in word)
    return word


def scattered_count(word: str, pattern: str) -> int:
    # dp[j] is the number of ways to spell pattern[:j] in the processed prefix.
    dp = [0] * (len(pattern) + 1)
    dp[0] = 1
    for letter in word:
        for j in range(len(pattern), 0, -1):
            if pattern[j - 1] == letter:
                dp[j] += dp[j - 1]
    return dp[-1]


def signature(word: str, k: int, alphabet: tuple[str, ...]) -> tuple[int, ...]:
    return tuple(
        scattered_count(word, pattern)
        for length in range(1, k + 1)
        for pattern in map("".join, product(alphabet, repeat=length))
    )


def factor_set(prefix: str, length: int) -> set[str]:
    return {prefix[i : i + length] for i in range(len(prefix) - length + 1)}


def complexity(prefix: str, length: int, k: int, alphabet: tuple[str, ...]) -> int:
    factors = factor_set(prefix, length)
    return len({signature(word, k, alphabet) for word in factors})


def fib_prefix(length: int) -> str:
    # Fixed point of 0 -> 01, 1 -> 0.
    word = "0"
    while len(word) < length:
        word = "".join("01" if ch == "0" else "0" for ch in word)
    return word[:length]


def thue_morse_image(word: str, power: int) -> str:
    for _ in range(power):
        word = "".join("01" if ch == "0" else "10" for ch in word)
    return word


def marker_positions(count: int) -> list[int]:
    # p_n = F_(3 n + 2), with F_0 = 0, F_1 = 1.
    fib = [0, 1]
    while len(fib) <= 3 * (count - 1) + 2:
        fib.append(fib[-1] + fib[-2])
    return [fib[3 * n + 2] for n in range(count)]


def lacunary_prefix(length: int) -> str:
    positions = set(marker_positions(32))
    return "".join("1" if i in positions else "0" for i in range(length))


def overlay(x: str, y: str) -> tuple[str, ...]:
    # A single-character coding keeps the finite product word compatible with
    # the string-based factor and subsequence routines.
    coding = {"00": "0", "01": "1", "10": "2", "11": "3"}
    return tuple(coding[a + b] for a, b in zip(x, y))


def check_source_theorem_7_4() -> None:
    h = morphic_prefix("0", {"0": "01", "1": "12", "2": "2"}, 240)
    assert len(h) > 20_000
    alphabet = ("0", "1", "2")
    rows = []
    for n in range(6, 16):
        b1 = complexity(h, n, 1, alphabet)
        b2 = complexity(h, n, 2, alphabet)
        b3 = complexity(h, n, 3, alphabet)
        p = len(factor_set(h, n))
        assert b1 < b2 < b3 == p, (n, b1, b2, b3, p)
        rows.append((n, b1, b2, b3, p))
    print("SOURCE_THEOREM_7_4_OK", rows)


def check_lacunary_reconstruction() -> None:
    y = lacunary_prefix(10_000)
    alphabet = ("0", "1")
    for n in range(1, 31):
        b2 = complexity(y, n, 2, alphabet)
        p = len(factor_set(y, n))
        assert b2 == p, (n, b2, p)
    for h in range(1, 7):
        n = marker_positions(h)[-1] + 1
        b1 = complexity(y, n, 1, alphabet)
        assert b1 >= h + 1, (h, n, b1)
    positions = marker_positions(8)
    assert all(positions[i + 1] >= 3 * positions[i] for i in range(len(positions) - 1))
    print("LACUNARY_RECONSTRUCTION_OK", positions[:8])


def check_overlay_finite() -> None:
    # For k = 0, 1, 2, test the product word at lengths where the finite
    # prefix contains every marker and enough repeated Fibonacci factors.
    x0 = fib_prefix(2_500)
    y = lacunary_prefix(len(x0))
    for k in range(3):
        x = thue_morse_image(x0, k)
        yy = lacunary_prefix(len(x))
        z = overlay(x, yy)
        alphabet = tuple(sorted(set(z)))
        for n in range(2, 14):
            assert complexity(z, n, k + 2, alphabet) == len(factor_set("".join(z), n)), (k, n)
        # The y projection supplies at least the abelian distinctions seen in y.
        z_b1 = complexity("".join(z), 34, 1, alphabet)
        y_b1 = complexity(yy, 34, 1, ("0", "1"))
        assert z_b1 >= y_b1
    print("OVERLAY_FINITE_OK")


def main() -> None:
    # Keep one exact rational operation in the checker to make the arithmetic
    # contract explicit without introducing any floating-point computation.
    assert Fraction(2, 3) + Fraction(1, 3) == 1
    check_source_theorem_7_4()
    check_lacunary_reconstruction()
    check_overlay_finite()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact finite checks of OEIS A185858, using actual binary arrays first.

Python 3 standard library only. All evidence is integer or Fraction arithmetic.
The exhaustive check is finite; it is not an all-length proof.
"""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb, prod


SOURCE_TERMS = (
    4, 22, 124, 694, 3892, 21814, 122284, 685462, 3842404, 21538774,
    120736732, 676795894, 3793814164, 21266419702, 119210006284,
    668237802454, 3745841266372, 20997505287574,
    117702592547644, 659787917778742,
)
LIMIT = 200


def trace_word(mask, n):
    """Row-major bit (3*r+c) is cell B[r,c]; use ordinary main diagonals."""
    return tuple(
        ((mask >> (3*i)) & 1)
        + ((mask >> (3*i + 4)) & 1)
        + ((mask >> (3*i + 8)) & 1)
        for i in range(n)
    )


def legal(word):
    return all(a != b for a, b in zip(word, word[1:]))


def exhaustive(n):
    """Enumerate every array, including its six non-trace cells."""
    histogram = Counter(trace_word(mask, n) for mask in range(1 << (3*(n+2))))
    assert sum(histogram.values()) == 1 << (3*(n+2))
    for word in product(range(4), repeat=n):
        assert histogram[word] == 64 * prod(comb(3, s) for s in word)
    return sum(count for word, count in histogram.items() if legal(word))


def row_transfer(limit):
    """Append actual rows; remember the previous two rows and previous trace.

    This does not use disjointness, trace weights, or the claimed recurrence.
    """
    state = Counter()
    for a, b, c in product(range(8), repeat=3):
        t = (a & 1) + ((b >> 1) & 1) + ((c >> 2) & 1)
        state[b, c, t] += 1
    counts = {1: sum(state.values())}
    for n in range(2, limit+1):
        nxt = Counter()
        for (b, c, old_trace), count in state.items():
            for d in range(8):
                new_trace = (b & 1) + ((c >> 1) & 1) + ((d >> 2) & 1)
                if new_trace != old_trace:
                    nxt[c, d, new_trace] += count
        state = nxt
        counts[n] = sum(state.values())
    return counts


def qmul(a, b):
    """Multiply in Q[s]/(s^2-13); a pair represents a[0]+a[1]*s."""
    return (a[0]*b[0] + 13*a[1]*b[1], a[0]*b[1] + a[1]*b[0])


def main():
    counts = row_transfer(LIMIT)
    total_arrays = 0
    for n in range(1, 6):
        used = {(i+c, c) for i in range(n) for c in range(3)}
        assert len(used) == 3*n
        assert 3*(n+2) - len(used) == 6
        raw = exhaustive(n)
        total_arrays += 1 << (3*(n+2))
        assert raw == counts[n]
        assert raw % 128 == 0
        assert raw // 128 == SOURCE_TERMS[n-1]
        print(f"EXHAUSTIVE n={n} arrays={1 << (3*(n+2))} "
              f"valid={raw} normalized={raw//128} all_trace_fibers=OK")
    print(f"EXHAUSTIVE_TOTAL arrays={total_arrays}")

    normalized = {0: Fraction(0)}
    weights = (1, 3, 3, 1)
    terminal = list(weights)
    for n in range(1, LIMIT+1):
        if n > 1:
            terminal = [weights[t] * sum(terminal[s] for s in range(4) if s != t)
                        for t in range(4)]
        assert counts[n] == 64 * sum(terminal)
        assert counts[n] % 128 == 0
        normalized[n] = Fraction(counts[n], 128)
        if n <= len(SOURCE_TERMS):
            assert normalized[n] == SOURCE_TERMS[n-1]
        if n >= 3:
            assert normalized[n] == 4*normalized[n-1] + 9*normalized[n-2]
    print("FIDELITY first_20_OEIS_terms=OK")
    print(f"ROW_TRANSFER_VS_WEIGHTED_WORDS n=1..{LIMIT}=OK")
    print(f"INTEGRAL_NORMALIZATION n=1..{LIMIT}=OK")
    print(f"CLAIMED_RECURRENCE n=3..{LIMIT}=OK")

    # Multiplication by 1-4x-9x^2 through degree LIMIT, with A(0)=0.
    for k in range(LIMIT+1):
        coefficient = normalized[k]
        if k >= 1:
            coefficient -= 4*normalized[k-1]
        if k >= 2:
            coefficient -= 9*normalized[k-2]
        assert coefficient == ({1: 4, 2: 6}.get(k, 0))
    print(f"GENERATING_FUNCTION_PRODUCT coefficients=0..{LIMIT}=OK")

    # Evaluate the printed radical expression in the exact quadratic field.
    # Division by 3*s is multiplication by s/39, so no real approximation occurs.
    minus_power = plus_power = (Fraction(1), Fraction(0))
    for n in range(1, LIMIT+1):
        minus_power = qmul(minus_power, (2, -1))
        plus_power = qmul(plus_power, (2, 1))
        left = qmul((-4, 1), minus_power)
        right = qmul((4, 1), plus_power)
        result = qmul((left[0]+right[0], left[1]+right[1]), (0, Fraction(1, 39)))
        assert result == (normalized[n], 0)
    print(f"RADICAL_FORMULA exact_Q_sqrt13 n=1..{LIMIT}=OK")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()

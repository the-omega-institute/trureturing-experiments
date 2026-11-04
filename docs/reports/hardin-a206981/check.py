#!/usr/bin/env python3
"""Exact finite checks for Hardin's OEIS A206981, NAME and first FORMULA.

Program: SPDX-License-Identifier: MIT.
Quoted sequence data: R. H. Hardin, OEIS A206981, revision 17,
https://oeis.org/A206981, CC BY-SA 4.0.
No floating-point arithmetic or recurrence-generated oracle is used.
"""
from itertools import product

DIRECTIONS = ((0, 1), (1, 0), (1, 1), (1, -1))
FORBIDDEN = {(0, 1, 0), (1, 0, 1)}
SOURCE_DATA = (
    4, 16, 36, 100, 256, 676, 1764, 4624, 12100, 31684, 82944,
    217156, 568516, 1488400, 3896676, 10201636, 26708224, 69923044,
    183060900, 479259664, 1254718084, 3284894596, 8599965696,
    22515002500, 58945041796, 154320122896, 404015326884,
)


def segments(n, width=2):
    """Every contiguous triple in each stated straight direction, without wraparound."""
    return [
        tuple((r + t * dr) * width + c + t * dc for t in range(3))
        for r in range(n)
        for c in range(width)
        for dr, dc in DIRECTIONS
        if all(0 <= r + t * dr < n and 0 <= c + t * dc < width for t in range(3))
    ]


def literal_count(n):
    triples = segments(n)
    return sum(
        all(tuple((mask >> i) & 1 for i in triple) not in FORBIDDEN for triple in triples)
        for mask in range(1 << (2 * n))
    )


def row_dp(limit):
    """Count actual width-two arrays by their last two rows; no claimed recurrence."""
    rows = tuple(product((0, 1), repeat=2))
    out = {1: 4, 2: 16}
    states = {(r, s): 1 for r in rows for s in rows}
    for n in range(3, limit + 1):
        new = {(s, t): 0 for s in rows for t in rows}
        for (r, s), count in states.items():
            for t in rows:
                if all((r[c], s[c], t[c]) not in FORBIDDEN for c in range(2)):
                    new[(s, t)] += count
        states = new
        out[n] = sum(states.values())
    return out


def fib(n):
    u, v = 0, 1
    for _ in range(n):
        u, v = v, u + v
    return u


def polynomial_product(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def check_switches():
    checked = 0
    for n in range(1, 11):
        for w in product((0, 1), repeat=n):
            d = tuple(w[i] ^ w[i + 1] for i in range(n - 1))
            restored = [w[0]]
            for bit in d:
                restored.append(restored[-1] ^ bit)
            assert tuple(restored) == w
            column_valid = all(w[i:i + 3] not in FORBIDDEN for i in range(n - 2))
            switches_valid = all((d[i], d[i + 1]) != (1, 1) for i in range(n - 2))
            assert column_valid == switches_valid
            checked += 1
        for b in (0, 1):
            for d in product((0, 1), repeat=n - 1):
                restored = [b]
                for bit in d:
                    restored.append(restored[-1] ^ bit)
                assert tuple(restored[i] ^ restored[i + 1] for i in range(n - 1)) == d
    return checked


def main():
    assert segments(1) == [] and segments(2) == []
    literal = [literal_count(n) for n in range(1, 9)]
    assert tuple(literal) == SOURCE_DATA[:8]
    print("Literal four-direction counts n=1..8:", literal)
    print("Fidelity: reproduced the source's first 8 terms by exhaustive array enumeration.")
    counts = row_dp(210)
    assert all(counts[n] == literal[n - 1] for n in range(1, 9))
    assert tuple(counts[n] for n in range(1, 28)) == SOURCE_DATA
    print("Independent row-state DP agrees with enumeration and all 27 displayed source terms.")
    checked = check_switches()
    print(f"Switch encoding: both inverse checks and forbidden-pattern equivalence pass for {checked} words, lengths 1..10.")
    assert all(counts[n] == 4 * fib(n + 1) ** 2 for n in range(1, 211))
    assert all(counts[n] == 2 * counts[n - 1] + 2 * counts[n - 2] - counts[n - 3]
               for n in range(4, 211))
    print("Exact identity a(n)=4*F(n+1)^2 holds for n=1..210; recurrence holds for n=4..210.")
    denominator = [1, -2, -2, 1]
    numerator = [0, 4, 8, -4]
    assert polynomial_product([1, 1], [1, -3, 1]) == denominator
    # Offset one: the generating series has coefficient zero at degree zero.
    series = [0] + [counts[n] for n in range(1, 211)]
    for degree in range(211):
        actual = sum(coefficient * series[degree - j]
                     for j, coefficient in enumerate(denominator) if degree >= j)
        expected = numerator[degree] if degree < len(numerator) else 0
        assert actual == expected
    print("OGF product coefficients match 4*x+8*x^2-4*x^3 through degree 210; denominator factorization is exact.")
    print("Scope: finite checks only; the all-length array bridge and recurrence require a proof.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

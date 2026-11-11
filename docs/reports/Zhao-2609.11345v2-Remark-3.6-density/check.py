#!/usr/bin/env python3
"""Exact small-rank checks for Zhao, arXiv:2609.11345v2, Remark 3.6.

Coefficients are integer polynomials in q, represented sparsely by exponent.
The standard-basis multiplication is T_w T_s = T_ws on an ascent and
T_w T_s = (q-1)T_w + q T_ws on a descent.
"""

from itertools import permutations
from math import comb


Poly = dict[int, int]
Element = dict[tuple[int, ...], Poly]


def p_add(a: Poly, b: Poly, scale: int = 1) -> Poly:
    out = dict(a)
    for degree, value in b.items():
        out[degree] = out.get(degree, 0) + scale * value
        if out[degree] == 0:
            del out[degree]
    return out


def p_mul(a: Poly, b: Poly) -> Poly:
    out: Poly = {}
    for da, va in a.items():
        for db, vb in b.items():
            out[da + db] = out.get(da + db, 0) + va * vb
    return {degree: value for degree, value in out.items() if value}


def inv_count(permutation: tuple[int, ...]) -> int:
    return sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )


def right_simple(permutation: tuple[int, ...], simple: int) -> tuple[int, ...]:
    values = list(permutation)
    values[simple - 1], values[simple] = values[simple], values[simple - 1]
    return tuple(values)


def reduced_word(permutation: tuple[int, ...]) -> list[int]:
    """Return a reduced word whose product is the given permutation."""
    current = permutation
    word: list[int] = []
    while inv_count(current):
        for simple in range(1, len(current)):
            if current[simple - 1] > current[simple]:
                word.append(simple)
                current = right_simple(current, simple)
                break
    return list(reversed(word))


def multiply_simple(element: Element, simple: int) -> Element:
    out: Element = {}
    for permutation, coefficient in element.items():
        product = right_simple(permutation, simple)
        if inv_count(product) > inv_count(permutation):
            out[product] = p_add(out.get(product, {}), coefficient)
        else:
            # (q - 1) T_w + q T_{ws}.
            out[permutation] = p_add(out.get(permutation, {}), coefficient, -1)
            out[permutation] = p_add(
                out[permutation], {degree + 1: value for degree, value in coefficient.items()}
            )
            out[product] = p_add(
                out.get(product, {}), {degree + 1: value for degree, value in coefficient.items()}
            )
    return {permutation: coefficient for permutation, coefficient in out.items() if coefficient}


def multiply(left: Element, right: Element) -> Element:
    out: Element = {}
    for permutation, right_coefficient in right.items():
        term = left
        for simple in reduced_word(permutation):
            term = multiply_simple(term, simple)
        for product, coefficient in term.items():
            out[product] = p_add(
                out.get(product, {}), p_mul(coefficient, right_coefficient)
            )
    return {permutation: coefficient for permutation, coefficient in out.items() if coefficient}


def subtract(left: Element, right: Element) -> Element:
    out: Element = {}
    for permutation, coefficient in left.items():
        out[permutation] = p_add(out.get(permutation, {}), coefficient)
    for permutation, coefficient in right.items():
        out[permutation] = p_add(out.get(permutation, {}), coefficient, -1)
    return {permutation: coefficient for permutation, coefficient in out.items() if coefficient}


def shuffle_sums(n: int) -> list[Element]:
    sums: list[Element] = [{} for _ in range(comb(n, 2) + 1)]
    for permutation in permutations(range(1, n + 1)):
        sums[inv_count(permutation)][permutation] = {0: 1}
    return sums


def commutator(shuffles: list[Element], i: int, j: int) -> Element:
    return subtract(multiply(shuffles[i], shuffles[j]), multiply(shuffles[j], shuffles[i]))


def main() -> None:
    # Fidelity: the source's n=4 coefficient of T_2413 in [Y_1,Y_3] is q-1,
    # and [Y_1,Y_2] is nonzero for every n >= 4 (checked here through n=6).
    rank_four = shuffle_sums(4)
    coefficient = commutator(rank_four, 1, 3).get((2, 4, 1, 3))
    assert coefficient == {0: -1, 1: 1}, coefficient
    assert commutator(rank_four, 1, 2)

    for n in (4, 5, 6):
        shuffles = shuffle_sums(n)
        N = comb(n, 2)
        commuting: list[tuple[int, int]] = []
        for i in range(1, N + 1):
            for j in range(i + 1, N + 1):
                if not commutator(shuffles, i, j):
                    commuting.append((i, j))
        E = len(commuting)
        assert N - 1 <= E <= N - 1 + (4 * n - 14) ** 2

        if n >= 6:
            L = comb(n - 4, 2) + 3
            low = {(i, j) for i in range(1, L + 1) for j in range(i + 1, L + 1)}
            high_start = N - L
            high = {
                (i, j)
                for i in range(high_start, N)
                for j in range(i + 1, N)
            }
            assert not (low | high) & set(commuting), (n, commuting)

        print(f"n={n}: N={N}, commuting_pairs={E}, bound=[{N - 1}, {N - 1 + (4 * n - 14) ** 2}]")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

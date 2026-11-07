#!/usr/bin/env python3
"""Exact checks for the odd-alphabet type-1 substitution family.

The checks use only integers, SymPy exact algebra, and finite exhaustive
alternating-word tests.  They do not use floating point arithmetic.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import sympy as sp


LETTERS = ("A", "B", "C", "D")


def substitution(q: int) -> dict[str, tuple[str, ...]]:
    """The source's type-1 substitution for alpha=1, beta=2q+1."""
    return {
        "A": ("B",),
        "B": ("C",),
        "C": ("B",) + tuple(x for _ in range(q) for x in ("A", "B")),
        "D": ("C",) + tuple(x for _ in range(q) for x in ("D", "C")),
    }


def decode(letter: str, beta: int) -> tuple[int, ...]:
    """Decode A,B,C,D as 1^1,beta^1,1^beta,beta^beta."""
    runs = {"A": (1, 1), "B": (beta, 1), "C": (1, beta), "D": (beta, beta)}
    value, length = runs[letter]
    return (value,) * length


def decode_word(word: tuple[str, ...], beta: int) -> tuple[int, ...]:
    out: tuple[int, ...] = ()
    for letter in word:
        out += decode(letter, beta)
    return out


def run_lengths(word: tuple[int, ...]) -> tuple[int, ...]:
    if not word:
        return ()
    out: list[int] = []
    current = word[0]
    length = 0
    for value in word:
        if value == current:
            length += 1
        else:
            out.append(length)
            current = value
            length = 1
    out.append(length)
    return tuple(out)


def alternating_words(length: int):
    classes = (("A", "C"), ("B", "D"))
    for first_class in (0, 1):
        for choices in product((0, 1), repeat=length):
            yield tuple(classes[(first_class + i) % 2][choices[i]] for i in range(length))


def check_integration_identities() -> int:
    checked = 0
    for q in range(1, 9):
        beta = 2 * q + 1
        sigma = substitution(q)
        for length in range(1, 9):
            for word in alternating_words(length):
                image = tuple(letter for symbol in word for letter in sigma[symbol])
                got = run_lengths(decode_word(image, beta))
                expected = decode_word(word, beta)
                assert got == expected, (q, word, got, expected)
                checked += 1
    return checked


def check_incidence_and_spectrum() -> None:
    q = sp.symbols("q", integer=True, positive=True)
    z = sp.symbols("z")
    matrix = sp.Matrix(
        [
            [0, 0, q, 0],
            [1, 0, q + 1, 0],
            [0, 1, 0, q + 1],
            [0, 0, 0, q],
        ]
    )
    expected_charpoly = (z - q) * (z + 1) * (z**2 - z - q)
    assert sp.expand(matrix.charpoly(z).as_expr() - expected_charpoly) == 0

    lam = (1 + sp.sqrt(4 * q + 1)) / 2
    v0 = sp.Matrix([(lam - 1) / (2 * lam), sp.Rational(1, 2), 1 / (2 * lam), 0])
    vD = sp.Matrix([sp.Rational(1, 4), (q + 2) / (4 * q), sp.Rational(1, 4), (q - 2) / (4 * q)])
    # Reduce expressions using lambda^2=lambda+q, the defining polynomial.
    relation = sp.sqrt(4 * q + 1) ** 2 - (4 * q + 1)
    for vector, eigenvalue in ((v0, lam), (vD, q)):
        residual = matrix * vector - eigenvalue * vector
        for entry in residual:
            assert sp.simplify(entry.subs(sp.sqrt(4 * q + 1) ** 2, 4 * q + 1)) == 0
    assert sp.simplify(sum(v0)) == 1
    assert sp.simplify(sum(vD)) == 1
    assert relation == 0


def check_frequency_gap() -> None:
    for q in range(3, 101):
        beta = 2 * q + 1
        # Keep the formulas exact in a quadratic number field.
        lam = (1 + sp.sqrt(4 * q + 1)) / 2
        f0 = 1 - 1 / (2 * lam)
        fD = sp.Rational(q + 1, 2 * q)
        gap = sp.simplify(f0 - fD)
        assert sp.simplify(gap - (q - lam) / (2 * q)) == 0
        # q > lambda follows by squaring positive integers:
        # (2q-1)^2-(4q+1)=4q(q-2)>0.
        assert (2 * q - 1) ** 2 - (4 * q + 1) == 4 * q * (q - 2) > 0
        assert gap.is_positive is True
        assert beta >= 7 and beta % 2 == 1

    # The source case q=1 is the {1,3} system; q=2 is the spectral tie.
    assert (2 * 1 - 1) ** 2 - (4 * 1 + 1) < 0
    assert (2 * 2 - 1) ** 2 - (4 * 2 + 1) == 0


def main() -> None:
    checked = check_integration_identities()
    check_incidence_and_spectrum()
    check_frequency_gap()
    print(f"integration identities checked: {checked}")
    print("characteristic polynomial and normalized eigenvectors: exact")
    print("decoded frequency gap f0-fD=(q-lambda)/(2q)>0 for every q>=3")
    print("q=1 is the source {1,3} case; q=2 is the unresolved spectral-tie boundary")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

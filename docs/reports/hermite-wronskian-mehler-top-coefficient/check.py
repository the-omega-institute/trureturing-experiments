#!/usr/bin/env python3
"""Exact checks for Kuznetsov--Yuan, arXiv:2606.02405v1, Conjecture 2.

All arithmetic is SymPy exact polynomial/rational arithmetic.  The recurrence
implemented below is Proposition 1(iii), equation (7) in the source; its
coefficient identity is the finite coefficient form of equations (1)--(2).
"""

from __future__ import annotations

from itertools import combinations
from functools import lru_cache
from math import factorial

import sympy as sp


x, y = sp.symbols("x y")


@lru_cache(maxsize=None)
def hermite(n: int, var: sp.Symbol) -> sp.Expr:
    """Physicists' H_n, H_0=1 and H_{n+1}=2xH_n-H_n'."""
    if n < 0:
        return sp.Integer(0)
    return sp.hermite(n, var)


@lru_cache(maxsize=None)
def wronskian(indices: tuple[int, ...], var: sp.Symbol) -> sp.Expr:
    if not indices:
        return sp.Integer(1)
    cols = [hermite(n, var) for n in indices]
    matrix = sp.Matrix(
        [[sp.diff(col, var, row) for col in cols] for row in range(len(cols))]
    )
    return sp.expand(matrix.det())


def h_tuple(indices: tuple[int, ...], xv: sp.Symbol) -> sp.Expr:
    return wronskian(indices, xv)


@lru_cache(maxsize=None)
def ek_coeff(k: tuple[int, ...], n: int) -> sp.Expr:
    """Coefficient of z^n in source equation (1), with the omitted terms."""
    if n in k:
        return sp.Integer(0)
    denominator = 2 ** (n + len(k)) * factorial(n)
    for kj in k:
        denominator *= n - kj
    hk_n_x = wronskian(k + (n,), x)
    hk_n_y = wronskian(k + (n,), y)
    return sp.cancel(hk_n_x * hk_n_y / denominator)


@lru_cache(maxsize=None)
def e_coeff(n: int) -> sp.Expr:
    return sp.cancel(hermite(n, x) * hermite(n, y) / (2**n * factorial(n)))


@lru_cache(maxsize=None)
def q_coeffs(k: tuple[int, ...]) -> tuple[sp.Expr, ...]:
    """Build Q_0,...,Q_{k_l+1} using source equation (7)."""
    ell = len(k)
    out: list[sp.Expr] = []
    shifted = tuple(kj - 1 for kj in k)
    out.append(
        sp.expand((-2) ** ell * sp.prod(k) * h_tuple(shifted, x) * h_tuple(shifted, y))
    )
    for n in range(1, k[-1] + 2):
        term = ek_coeff(k, n)
        for m in range(1, n + 1):
            term -= e_coeff(m) * out[n - m]
        out.append(sp.expand(term))
    return tuple(out)


def falling(k: int, r: int) -> int:
    if r < 0:
        return 0
    if r > k:
        return 0
    result = 1
    for j in range(r):
        result *= k - j
    return result


def proposition3(k: int, n: int) -> sp.Expr:
    first = (
        2 ** (n - 1)
        * falling(k, n - 1)
        * hermite(k + 1 - n, x)
        * hermite(k + 1 - n, y)
        if n >= 1
        else 0
    )
    second = (
        2 ** (n + 1)
        * falling(k, n + 1)
        * hermite(k - 1 - n, x)
        * hermite(k - 1 - n, y)
    )
    return sp.expand(first - second)


@lru_cache(maxsize=None)
def finite_formula(k: tuple[int, ...], r: int) -> sp.Expr:
    """The finite coefficient formula obtained from the source recurrence."""
    ell = len(k)
    result = 0
    for j, kj in enumerate(k):
        prefix = k[:j] + k[j + 1 :]
        cj = 2 ** (kj + ell - 1) * factorial(kj)
        cj *= sp.prod(kj - ka for a, ka in enumerate(k) if a != j)
        def e(m: int) -> sp.Expr:
            return e_coeff(m) if m >= 0 else sp.Integer(0)
        result += cj * h_tuple(prefix, x) * h_tuple(prefix, y) * (e(kj + 1 - r) - e(kj - 1 - r))
    return sp.expand(result)


def assert_zero(expr: sp.Expr, label: str) -> None:
    reduced = sp.Poly(sp.expand(expr), x, y, domain=sp.QQ)
    if not reduced.is_zero:
        raise AssertionError(f"{label}: residual {reduced.as_expr()}")


def main() -> None:
    # Fidelity: Proposition 3, including the smallest allowed positive k.
    prop3_cases = 0
    for k in range(1, 7):
        qs = q_coeffs((k,))
        for n, q in enumerate(qs):
            assert_zero(q - proposition3(k, n), f"Proposition 3 k={k} n={n}")
            prop3_cases += 1

    tuples = [
        k for ell in (2, 3, 4) for k in combinations(range(1, 7), ell)
    ]
    highest_cases = 0
    recurrence_cases = 0
    finite_cases = 0
    for k in tuples:
        qs = q_coeffs(k)
        # Equation (2): E_k = E * sum z^r Q_r, coefficient by coefficient.
        for n, q in enumerate(qs):
            convolution = sum(e_coeff(m) * qs[n - m] for m in range(n + 1))
            assert_zero(ek_coeff(k, n) - convolution, f"equation (2) k={k} n={n}")
            recurrence_cases += 1
            if n >= 1:
                assert_zero(q - finite_formula(k, n), f"finite formula k={k} r={n}")
                finite_cases += 1
        ell = len(k)
        expected = (
            2 ** (k[-1] + ell - 1)
            * factorial(k[-1])
            * sp.prod(k[-1] - kj for kj in k[:-1])
            * h_tuple(k[:-1], x)
            * h_tuple(k[:-1], y)
        )
        assert_zero(qs[k[-1] + 1] - expected, f"Conjecture 2 k={k}")
        highest_cases += 1

    print(f"source Proposition 3 cases: {prop3_cases}")
    print(f"equation (2) coefficient cases: {recurrence_cases}")
    print(f"finite coefficient formula cases: {finite_cases}")
    print(f"Conjecture 2 highest-coefficient cases: {highest_cases}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact checks for the Stessin--Yang Problem 7.12(a) counterexample.

The numeric Clifford checks use sparse matrices over the Gaussian integers
(pairs of Python integers).  The small symbolic checks use SymPy only with
exact symbols and integer coefficients; no floating point arithmetic occurs.
"""

from itertools import combinations
from pathlib import Path

import sympy as sp


Gaussian = tuple[int, int]
Sparse = dict[tuple[int, int], Gaussian]
ZERO: Gaussian = (0, 0)


def gadd(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] + y[0], x[1] + y[1])


def gneg(x: Gaussian) -> Gaussian:
    return (-x[0], -x[1])


def gmul(x: Gaussian, y: Gaussian) -> Gaussian:
    return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])


def gscale(k: int, x: Gaussian) -> Gaussian:
    return (k * x[0], k * x[1])


def put(A: Sparse, key: tuple[int, int], value: Gaussian) -> None:
    if value == ZERO:
        return
    old = A.get(key, ZERO)
    new = gadd(old, value)
    if new == ZERO:
        A.pop(key, None)
    else:
        A[key] = new


def madd(A: Sparse, B: Sparse) -> Sparse:
    C = dict(A)
    for key, value in B.items():
        put(C, key, value)
    return C


def mscale(k: int, A: Sparse) -> Sparse:
    return {key: gscale(k, value) for key, value in A.items() if gscale(k, value) != ZERO}


def mmul(A: Sparse, B: Sparse) -> Sparse:
    by_row: dict[int, list[tuple[int, Gaussian]]] = {}
    for (row, col), value in B.items():
        by_row.setdefault(row, []).append((col, value))
    C: Sparse = {}
    for (row, mid), left in A.items():
        for col, right in by_row.get(mid, []):
            put(C, (row, col), gmul(left, right))
    return C


def identity(d: int) -> Sparse:
    return {(i, i): (1, 0) for i in range(d)}


def kron(A: Sparse, B: Sparse, da: int, db: int) -> Sparse:
    C: Sparse = {}
    for (i, j), x in A.items():
        for (k, l), y in B.items():
            put(C, (i * db + k, j * db + l), gmul(x, y))
    return C


def paulis() -> tuple[Sparse, Sparse, Sparse, Sparse]:
    one = {(0, 0): (1, 0), (1, 1): (1, 0)}
    x = {(0, 1): (1, 0), (1, 0): (1, 0)}
    y = {(0, 1): (0, -1), (1, 0): (0, 1)}
    z = {(0, 0): (1, 0), (1, 1): (-1, 0)}
    return one, x, y, z


def clifford_generators(m: int) -> list[Sparse]:
    one, x, y, z = paulis()
    result: list[Sparse] = []
    for j in range(m):
        for middle in (x, y):
            value = {(0, 0): (1, 0)}
            size = 1
            for k in range(m):
                factor = z if k < j else middle if k == j else one
                value = kron(value, factor, size, 2)
                size *= 2
            result.append(value)
    return result


def monomial(gammas: list[Sparse], subset: tuple[int, ...]) -> Sparse:
    value = identity(2 ** (len(gammas) // 2))
    for i in subset:
        value = mmul(value, gammas[i])
    return value


def is_scalar(A: Sparse, d: int, scalar: Gaussian) -> bool:
    expected = {(i, i): scalar for i in range(d)}
    return A == expected


def to_sympy(A: Sparse, d: int) -> sp.Matrix:
    return sp.Matrix(
        d,
        d,
        lambda i, j: A.get((i, j), ZERO)[0] + sp.I * A.get((i, j), ZERO)[1],
    )


def exact_pauli_fidelity() -> None:
    # This is the source's Example 2.2 (Pauli matrices and its Q_sigma).
    z0, z1, z2, z3 = sp.symbols("z0 z1 z2 z3")
    _, x, y, z = paulis()
    expr = sp.expand((z0 * sp.eye(2) + z1 * to_sympy(x, 2) + z2 * to_sympy(y, 2) + z3 * to_sympy(z, 2)).det())
    assert expr == z0**2 - z1**2 - z2**2 - z3**2
    assert to_sympy(mmul(x, x), 2) == sp.eye(2)
    assert to_sympy(madd(mmul(x, y), mmul(y, x)), 2) == sp.zeros(2)
    print("source Pauli characteristic polynomial and anticommutation: OK")


def check_family(m: int, symbolic: bool = False) -> None:
    gammas = clifford_generators(m)
    d = 2**m
    n = 2 * m
    unit = identity(d)
    for i, g in enumerate(gammas):
        assert mmul(g, g) == unit, (m, i, "involution")
        for j in range(i):
            assert madd(mmul(g, gammas[j]), mmul(gammas[j], g)) == {}

    subsets = [tuple(i for i in range(n) if mask & (1 << i)) for mask in range(2**n)]
    for subset in subsets[1:]:
        h = monomial(gammas, subset)
        assert sum(value[0] for (row, col), value in h.items() if row == col) == 0
    print(f"m={m}, d={d}: {len(subsets)} Clifford monomials, nonidentity traces zero")

    # Integer coefficients make every check exact.  The quartic is the one
    # in the proposed proof: (H^2-(a+b+c t^2)I)^2 - 4 a c t^2 I = 0.
    coefficients = list(range(1, n + 1))
    t = 3
    for subset in subsets:
        h = monomial(gammas, subset)
        size = len(subset)
        c = 1 if (size * (size - 1) // 2) % 2 == 0 else -1
        commuting: list[int] = []
        anticommuting: list[int] = []
        for i, g in enumerate(gammas):
            gh, hg = mmul(g, h), mmul(h, g)
            if gh == hg:
                commuting.append(i)
            else:
                assert madd(gh, hg) == {}
                anticommuting.append(i)
        H = mscale(t, h)
        U: Sparse = {}
        for i in range(n):
            H = madd(H, mscale(coefficients[i], gammas[i]))
        for i in commuting:
            U = madd(U, mscale(coefficients[i], gammas[i]))
        a = sum(coefficients[i] ** 2 for i in commuting)
        b = sum(coefficients[i] ** 2 for i in anticommuting)
        H2 = mmul(H, H)
        centered = madd(H2, mscale(-(a + b + c * t * t), unit))
        residual = madd(mmul(centered, centered), mscale(-4 * a * c * t * t, unit))
        assert residual == {}, (m, subset, residual)
        assert is_scalar(mmul(h, h), d, (c, 0))
    print(f"m={m}, d={d}: quartic residual zero for every {len(subsets)} monomial")

    if symbolic:
        # Representative subset sizes 0,...,6 in dimension 8.  The factors
        # have degree 2 with exponent 4 or degree 4 with exponent 2.
        z0, tvar, lam = sp.symbols("z0 t lam")
        zvars = sp.symbols("z1:7")
        s_gammas = [to_sympy(g, d) for g in gammas]
        for k in range(7):
            subset = tuple(range(k))
            h = to_sympy(monomial(gammas, subset), d)
            H = tvar * h
            for i, g in enumerate(s_gammas):
                H += zvars[i] * g
            # charpoly uses its generator as an indeterminate; substitute
            # after the call so that SymPy keeps the same z0 symbol.
            char = sp.expand(H.charpoly(lam).as_expr().subs(lam, -z0))
            factors = sp.factor_list(char)[1]
            expected_exp = 4 if k in (0, 1, 6) else 2
            expected_deg = 2 if expected_exp == 4 else 4
            assert len(factors) == 1
            factor, exponent = factors[0]
            assert exponent == expected_exp and sp.degree(sp.expand(factor), z0) == expected_deg
        print("m=3, d=8: symbolic characteristic factors match fourth-power/square pattern")


def main() -> None:
    exact_pauli_fidelity()
    check_family(3, symbolic=True)
    check_family(4, symbolic=False)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

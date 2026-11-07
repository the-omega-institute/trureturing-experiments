#!/usr/bin/env python3
"""Exact checks for the Mayer--Knill scaling transition candidate.

All numerical data are SymPy Integers/Rationals (with exact square roots for
Kraus amplitudes); no floating point values are used as evidence.
"""

from __future__ import annotations

import sympy as sp


def inner(x: sp.Matrix, y: sp.Matrix) -> sp.Expr:
    return (x.conjugate().T * y)[0]


def projector(x: sp.Matrix) -> sp.Matrix:
    return x * x.conjugate().T


def channel(kraus: list[sp.Matrix], rho: sp.Matrix) -> sp.Matrix:
    return sum((k * rho * k.conjugate().T for k in kraus), sp.zeros(rho.rows))


def process_fidelity(kraus: list[sp.Matrix], d: int) -> sp.Expr:
    # F(E) = d^(-2) sum_a |tr K_a|^2 for the normalized maximally entangled state.
    return sp.simplify(sum(abs(sp.trace(k)) ** 2 for k in kraus) / d**2)


def main() -> None:
    d = 2
    epsilon = sp.Rational(1, 10)
    p = d * epsilon / (d - 1)

    # The depolarizing witness E(rho)=(1-p)rho+p Tr(rho)I/d.
    q0 = 1 - 3 * p / 4
    q = p / 4
    I = sp.eye(2)
    X = sp.Matrix([[0, 1], [1, 0]])
    Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    Z = sp.diag(1, -1)
    kraus = [sp.sqrt(q0) * I, sp.sqrt(q) * X, sp.sqrt(q) * Y, sp.sqrt(q) * Z]
    tp = sp.simplify(sum((k.conjugate().T * k for k in kraus), sp.zeros(d)))
    assert tp == I

    e0 = sp.Matrix([1, 0])
    e1 = sp.Matrix([0, 1])
    plus = sp.Matrix([1, 1]) / sp.sqrt(2)
    probes = [e0, e1, plus]
    for psi in probes:
        survival = sp.simplify(inner(psi, channel(kraus, projector(psi)) * psi))
        assert sp.simplify(1 - survival - epsilon) == 0

    F = process_fidelity(kraus, d)
    assert sp.simplify(F - (1 - sp.Rational(3, 2) * epsilon)) == 0
    print("depolarizing_fidelity", F)
    print("depolarizing_input_error", epsilon)

    # A connected, spanning two-state family with a non-scalar obstruction.
    # Columns of V are the two unit probes; A=V diag(1,-1) V^{-1} preserves
    # each tested line, is traceless, and is therefore a square-root witness.
    v = sp.Matrix([sp.Rational(3, 5), sp.Rational(4, 5)])
    V = sp.Matrix.hstack(e0, v)
    A = sp.simplify(V * sp.diag(1, -1) * V.inv())
    assert sp.simplify(inner(e0, v)) == sp.Rational(3, 5)
    assert sp.simplify(A * e0 - e0) == sp.zeros(2, 1)
    assert sp.simplify(A * v + v) == sp.zeros(2, 1)
    assert sp.trace(A) == 0
    assert A != sp.zeros(2)
    print("connected_graph_overlap", sp.Rational(3, 5))
    print("traceless_line_preserving_A", A)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

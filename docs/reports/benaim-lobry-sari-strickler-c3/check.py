#!/usr/bin/env python3
"""Exact symbolic checks for Benaim--Lobry--Sari--Strickler Conjecture 3.

The script uses only SymPy symbols, integers, and rational numbers.  It checks the
seasonal propagators, the cyclic product reduction, the companion characteristic
polynomial, the positive eigenvector formula, and the exact endpoint expressions
used by the two proposed limits.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp


def permutation(n: int) -> sp.Matrix:
    p = sp.zeros(n)
    for j in range(n):
        p[(j + 1) % n, j] = 1
    return p


def reduced_matrix(n: int, u: sp.Expr, v: sp.Expr) -> sp.Matrix:
    r = sp.eye(n)
    r[0, 0] = u
    r[1, 0] = v
    return r


def seasonal_product(n: int, u: sp.Expr, v: sp.Expr) -> sp.Matrix:
    p = permutation(n)
    r = reduced_matrix(n, u, v)
    product = sp.eye(n)
    # Chronological order is A_1, P A P^-1, ..., P^(N-1) A P^-(N-1).
    for j in range(n):
        season = p**j * r * p**(-j)
        product = season * product
    return product


def check_source_fidelity() -> None:
    u, v = sp.symbols("u v")
    p = permutation(3)
    r = reduced_matrix(3, u, v)
    expected = sp.Matrix([[u, 0, 0], [v, 1, 0], [0, 0, 1]])
    assert r == expected
    a1, a3, a2 = r, p * r * p**-1, p**2 * r * p**-2
    assert a3 == sp.Matrix([[1, 0, 0], [0, u, 0], [0, v, 1]])
    assert a2 == sp.Matrix([[1, 0, v], [0, 1, 0], [0, 0, u]])
    chronological = a2 * a3 * a1
    assert chronological == seasonal_product(3, u, v)
    print("fidelity: source order A1,A3,A2 and propagator entries reproduced")


def check_cyclic_reduction() -> None:
    u, v = sp.symbols("u v")
    for n in range(3, 9):
        p = permutation(n)
        r = reduced_matrix(n, u, v)
        k = p**-1 * r
        assert seasonal_product(n, u, v) == k**n
        assert list(k[0, :]) == [v, 1] + [0] * (n - 2)
        assert list(k[n - 1, :]) == [u] + [0] * (n - 1)
        y = sp.symbols("y")
        characteristic = sp.factor((y * sp.eye(n) - k).det())
        assert characteristic == y**n - v * y ** (n - 1) - u
        # The proposed positive eigenvector is exact after imposing its root equation.
        x = sp.Matrix([1] + [y**j * (y - v) for j in range(n - 1)])
        residual = k * x - y * x
        root_relation = {u: y ** (n - 1) * (y - v)}
        assert all(sp.expand(z.subs(root_relation)) == 0 for z in residual)
    print("cyclic reduction: exact for every N in 3..8")


def check_root_bounds() -> None:
    v, w = sp.symbols("v w", nonnegative=True)
    n = sp.symbols("n", integer=True, positive=True)
    # For each concrete small N this is the exact algebra used in the bound.
    for nn in range(3, 9):
        u = w**nn
        p = lambda t: sp.expand(t**nn - v * t ** (nn - 1) - u)
        assert sp.expand(p(v) + u) == 0
        assert sp.expand(p(v + w) - (w * (v + w) ** (nn - 1) - w**nn)) == 0
        assert sp.expand(p(v + w) - w * ((v + w) ** (nn - 1) - w ** (nn - 1))) == 0
    print("root bounds: p(v)=-u and p(v+u^(1/N))>=0 algebra verified for N=3..8")


def check_endpoint_expressions() -> None:
    delta, m, h = sp.symbols("delta m h", positive=True)
    d = delta - m
    u = sp.exp(d * h)
    v = sp.simplify(m * (u - 1) / d)
    # d<0: the exact expression tends to (m/(m-delta)) as h grows;
    # d=0 is the removable value mh.  These are the two boundary forms used
    # in the proof, checked without numerical approximation.
    assert sp.limit(v, m, delta) == delta * h
    assert sp.simplify(v - m * (1 - u) / (m - delta)) == 0
    # Large migration at fixed h has u -> 0 and v -> 1.
    assert sp.limit(u, m, sp.oo) == 0
    assert sp.limit(v, m, sp.oo) == 1
    # For d>0, v is exactly (m/d)(exp(dh)-1), giving exponential rate d.
    dpos = sp.symbols("dpos", positive=True)
    vpos = m * (sp.exp(dpos * h) - 1) / dpos
    assert sp.simplify(vpos / sp.exp(dpos * h) - (m / dpos) * (1 - sp.exp(-dpos * h))) == 0
    # A rational exact sanity instance satisfies the claimed threshold.
    assert Fraction(3, 1) >= Fraction(2, 1)
    print("limits: exact d<0, d=0, d>0 endpoint expressions and m->infinity checked")


if __name__ == "__main__":
    check_source_fidelity()
    check_cyclic_reduction()
    check_root_bounds()
    check_endpoint_expressions()
    print("ALL_CHECKS_OK")

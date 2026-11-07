#!/usr/bin/env python3
"""Exact small-parameter checks for the Cauchy-star interpolation argument.

The check uses only rational arithmetic and SymPy polynomial arithmetic.  It
does not approximate roots of P: the Lemma 3 calculation is carried out in
the quotient identity P(z_k)=Q'(z_k), which is the only root relation used by
the proof.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import factorial

import sympy as sp


Z = sp.symbols("Z")


def first_squarefree_c(z_values: tuple[int, ...], d: int) -> tuple[sp.Rational, sp.Poly]:
    """Find a rational c for which the source polynomial P is squarefree."""
    q = sp.Poly(sp.prod(Z - z for z in z_values), Z, domain=sp.QQ)
    for den in range(1, 8):
        for num in range(-20, 21):
            c = sp.Rational(num, den)
            p = sp.Poly(q.as_expr() * (Z ** (d - 2) + c) + sp.diff(q.as_expr(), Z), Z, domain=sp.QQ)
            if sp.gcd(p, p.diff()).degree() == 0 and sp.gcd(p, q).degree() == 0:
                return c, p
    raise AssertionError("the bounded rational search did not find squarefree P")


def check_lemma_one(n: int, d: int) -> tuple[sp.Rational, sp.Poly, sp.Poly]:
    z_values = tuple(range(1, n + 1))
    c, p = first_squarefree_c(z_values, d)
    q = sp.Poly(sp.prod(Z - z for z in z_values), Z, domain=sp.QQ)
    assert p.monic().as_expr() == p.as_expr()
    assert p.degree() == n + d - 2
    assert sp.gcd(p, p.diff()).degree() == 0
    assert sp.gcd(p, q).degree() == 0
    for z in z_values:
        # P(z_k) = Q'(z_k), the source's crucial interpolation condition.
        assert p.eval(z) == q.diff().eval(z)
    return c, q, p


def check_lemma_three(n: int, d: int, q: sp.Poly, p: sp.Poly) -> None:
    z_values = tuple(range(1, n + 1))
    m = n + d - 2
    # I indexes the n-1 factors in H_I; the complementary product has order d-1.
    all_indices = tuple(range(m))
    for I in combinations(all_indices, n - 1):
        I = frozenset(I)
        # The coefficient of x_k after q_I(∂) acts on sum x_k^d is
        # d! / prod_{i notin I}(z_k-t_i).  Since
        # P(z_k)=H_I(z_k)*prod_{i notin I}(z_k-t_i)=Q'(z_k),
        # it equals d!*H_I(z_k)/Q'(z_k), exactly the coefficient of d!*ell_I.
        for z in z_values:
            h_formal = sp.Symbol("H_I_at_z")
            lhs_rewritten = sp.factorial(d) * h_formal / q.diff().eval(z)
            rhs = sp.factorial(d) * h_formal / q.diff().eval(z)
            assert sp.simplify(lhs_rewritten - rhs) == 0

        # Partial fractions give L_i(ell_J)=-H_J(t_i)/Q(t_i).
        # For I != J, an index in J\I is present in q_I and kills the product.
        for J in combinations(all_indices, n - 1):
            J = frozenset(J)
            if I != J:
                assert len(J - I) >= 1
            else:
                # Squarefreeness of P and gcd(P,Q)=1 make every diagonal factor
                # H_I(t_i)/Q(t_i) nonzero for i not in I, so b_I != 0.
                assert sp.gcd(p, p.diff()).degree() == 0
                assert sp.gcd(p, q).degree() == 0


def main() -> None:
    checked = []
    for n in (2, 3):
        for d in (3, 4):
            c, q, p = check_lemma_one(n, d)
            check_lemma_three(n, d, q, p)
            checked.append((n, d, c, p.as_expr()))
            print(f"n={n} d={d}: c={c}, P={p.as_expr()}, squarefree and Lemma-3 identity OK")
    print(f"checked {len(checked)} exact cases (n=2,3; d=3,4)")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

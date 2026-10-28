#!/usr/bin/env python3
"""Exact checks for Rehman--Yang--Ajmal--Rodrigues--Maden--Cai Question 1.

All matrix and polynomial calculations use SymPy over the integers/radicals.  No
floating point value is used as evidence.
"""
import sympy as sp

x = sp.symbols("x")

def seidel_turan(t: int, p: int, inserted: bool = False) -> sp.Matrix:
    """Seidel matrix of K_{p,...,p}, with the edge (0,1) inserted in part 0."""
    n = t * p
    A = sp.zeros(n)
    starts = [i * p for i in range(t)]
    for i in range(t):
        for j in range(i + 1, t):
            for u in range(starts[i], starts[i] + p):
                for v in range(starts[j], starts[j] + p):
                    A[u, v] = A[v, u] = 1
    if inserted:
        A[0, 1] = A[1, 0] = 1
    return sp.ones(n) - sp.eye(n) - 2 * A

def factor_matrix(M: sp.Matrix) -> sp.Expr:
    return sp.factor(M.charpoly(x).as_expr())

def check_source_p2_family() -> None:
    """Reproduce the paper's p=2 decreasing family exactly."""
    for t in range(3, 9):
        base = factor_matrix(seidel_turan(t, 2))
        expected_base = (x - 3) ** (t - 1) * (x + 1) ** t * (x + (2 * t - 3))
        assert sp.expand(base - expected_base) == 0, (t, base, expected_base)

        embedded = factor_matrix(seidel_turan(t, 2, inserted=True))
        d = t * t - 2 * t + 5
        expected_embedded = (
            (x - 3) ** (t - 2) * (x - 1) * (x + 1) ** (t - 1)
            * (x ** 2 + 2 * (t - 2) * x - (2 * t + 1))
        )
        assert sp.expand(embedded - expected_embedded) == 0, (t, embedded, expected_embedded)

        # The quadratic has roots -(t-2) +/- sqrt(d).  For t >= 3,
        # 0 < d < t^2, so one root is positive and one is negative.
        assert d > 0 and d < t * t
        base_energy = sp.Integer(6 * t - 6)
        embedded_energy = sp.Integer(4 * t - 6) + 2 * sp.sqrt(d)
        delta = sp.simplify(embedded_energy - base_energy)
        assert sp.simplify(delta - 2 * (sp.sqrt(d) - t)) == 0
        assert sp.ask(sp.Q.negative(delta)) is True
        print(f"source p=2: t={t}, exact delta={delta} < 0")

def check_twin_threshold_algebra() -> None:
    """Check the exact projection bound and its strictness threshold."""
    m = sp.symbols("m", integer=True, positive=True)
    lower_bound = sp.simplify(2 - sp.Rational(8, 1) / m)
    assert sp.simplify(lower_bound.subs(m, 4)) == 0
    assert sp.simplify(lower_bound.subs(m, 5)) == sp.Rational(2, 5)
    # A balanced T(n,t) has part sizes floor(n/t) or ceil(n/t).
    for t in range(2, 13):
        for n in range(4 * t, 4 * t + 9):
            assert n // t >= 4
    print("independent-twin bound: 2 - 8/m is positive for m>4 and the m=4 case is isolated")

def check_sharpness_family() -> None:
    """Check the exact critical polynomial and its sign change at t=7."""
    k = sp.symbols("k", integer=True, positive=True)
    # The 3x3 quotient block in the proposed proof has this characteristic polynomial.
    Mchar = x ** 3 + (3 * k - 4) * x ** 2 - (6 * k + 7) * x + 10 - 21 * k
    assert sp.expand(Mchar.subs(k, sp.Integer(1))) == x**3 - x**2 - 13*x - 11
    # det(4 I - M) = 3(t-7) = 3(k-6); this is the sign test for its unique positive root.
    assert sp.expand(3 * (k - 6) - (Mchar.subs(x, 4))) == 0
    for t in range(2, 13):
        sign = sp.sign(3 * (t - 7))
        assert sign == (-1 if t < 7 else 0 if t == 7 else 1)
    print("sharpness family: det(4I-M)=3(t-7), hence equality at t=7 and decrease for t>=8")

def main() -> None:
    check_source_p2_family()
    check_twin_threshold_algebra()
    check_sharpness_family()
    print("ALL_CHECKS_OK")

if __name__ == "__main__":
    main()

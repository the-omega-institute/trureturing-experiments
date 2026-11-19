"""Exact checks for arXiv:2410.18867v1, Sections 2 and 4."""

from math import factorial, prod

from sympy import Matrix, Poly, QQ, Rational, S, cancel, diff, fraction, integrate, symbols
from sympy.polys.matrices import DomainMatrix


t = symbols("t")


def wronskian(functions):
    """Differentiate every column, clear row denominators and take the determinant."""
    rows = []
    determinant_denominator = S.One
    derivatives = list(functions)
    for _ in functions:
        row_denominator = Poly(1, t, domain=QQ)
        for derivative in derivatives:
            denominator = fraction(derivative)[1]
            row_denominator = row_denominator.lcm(Poly(denominator, t, domain=QQ))
        factor = row_denominator.as_expr()
        rows.append([cancel(factor * derivative) for derivative in derivatives])
        determinant_denominator *= factor
        derivatives = [cancel(diff(derivative, t)) for derivative in derivatives]

    polynomial_ring = QQ.poly_ring(t)
    cleared_matrix = DomainMatrix.from_Matrix(Matrix(rows)).convert_to(polynomial_ring)
    numerator = polynomial_ring.to_sympy(cleared_matrix.det())
    return cancel(numerator / determinant_denominator)


def check_family(poles, dimensions):
    q = Poly(prod(t - pole for pole in poles), t, domain=QQ)
    h0 = integrate(q.as_expr() ** 2, t)
    bad_shifts = {-h0.subs(t, pole) for pole in poles}
    shift = next(value for value in range(len(poles) + 1) if value not in bad_shifts)
    h = h0 + shift
    assert Poly(diff(h, t), t, domain=QQ) == q**2
    assert all(h.subs(t, pole) != 0 for pole in poles)

    for n in dimensions:
        functions = [cancel(h**k / q.as_expr() ** (n - 1)) for k in range(n)]
        reduced = [fraction(function) for function in functions]
        denominators = [Poly(den, t, domain=QQ).monic() for _, den in reduced]
        assert all(denominator == q ** (n - 1) for denominator in denominators)
        degrees = {Poly(num, t, domain=QQ).degree() for num, _ in reduced}
        assert len(degrees) == n  # Normalize the common denominator to infer independence.
        expected = prod(factorial(k) for k in range(n))
        assert wronskian(functions) == expected
        print(f"poles={poles} shift={shift} n={n} pole_order={n - 1} W={expected}")


def main():
    example = (t**2 + t, 2 * t**2, t - 2)
    assert wronskian(example) == -8
    print("Section 2 polynomial example: W=-8")
    check_family((-1, 1), range(2, 6))
    check_family((0, 1), range(2, 5))
    check_family((Rational(-3, 2), Rational(1, 3), Rational(5, 2)), range(2, 5))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

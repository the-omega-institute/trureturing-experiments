"""Exact symbolic checks for Shen Shan's Problem 8.1.2."""

from sympy import (
    E, Function, Integer, Poly, diff, expand, exp, limit, oo, simplify, solve, symbols,
)


def require_equal(actual, expected):
    if simplify(actual - expected) != 0:
        raise AssertionError(f"expected {expected}, got {actual}")


def check_source_laplace_derivatives():
    """Thesis (2.2) for a positive measure supported at two positive points."""
    x, a, b, u, v = symbols("x a b u v", positive=True)
    n = symbols("n", integer=True, nonnegative=True)
    laplace = a * exp(-u * x) + b * exp(-v * x)
    signed_moment = a * u**n * exp(-u * x) + b * v**n * exp(-v * x)
    require_equal(signed_moment.subs(n, 0), laplace)
    require_equal(-diff(signed_moment, x), signed_moment.subs(n, n + 1))
    if signed_moment.is_positive is not True:
        raise AssertionError("the two positive atoms must give positive signed moments")
    print("Fidelity (thesis (2.2)): two-atom Laplace derivative base and recurrence agree.")


def check_inverse_derivatives():
    """Differentiate the formal identity H(R(y)) = y to third order."""
    y = symbols("y")
    forward = Function("H")
    inverse = Function("R")(y)
    first, second, third = (diff(forward(inverse), y, n) for n in range(1, 4))
    d1, d2, d3 = (diff(inverse, y, n) for n in range(1, 4))
    solved, = solve([first - 1, second, third], [d1, d2, d3], dict=True)
    h1, h2, h3 = (diff(forward(inverse), inverse, n) for n in range(1, 4))
    require_equal(solved[d1], 1 / h1)
    require_equal(solved[d2], -h2 / h1**3)
    require_equal(solved[d3], (3 * h2**2 - h1 * h3) / h1**5)
    print("Inverse derivative identities through order three agree.")


def check_completely_monotone_family():
    x, t = symbols("x T", positive=True)
    s = symbols("s", nonnegative=True)
    weight = exp(t - 1) / (4 * t)
    f = exp(-x / t) + weight * exp(-x)
    h = f / x
    f_zero = f.subs(x, 0)
    require_equal(f_zero, 1 + weight)
    if f_zero.is_positive is not True or diff(h, x).is_negative is not True:
        raise AssertionError("the quotient must have a positive numerator and decrease")
    if limit(h.subs(t, 30), x, 0, dir="+") != oo:
        raise AssertionError("wrong lower endpoint at T = 30")
    require_equal(limit(h.subs(t, 30), x, oo), 0)

    a = 2 + (t + 1) / (4 * t)
    b = 5 + (t**2 + 2 * t + 2) / (4 * t)
    c = 16 + (t**3 + 3 * t**2 + 6 * t + 6) / (4 * t)
    h1, h2, h3 = (simplify(diff(h, x, n).subs(x, t)) for n in range(1, 4))
    require_equal(h1, -a / (E * t**2))
    require_equal(h2, b / (E * t**3))
    require_equal(h3, -c / (E * t**4))
    polynomial = 6 * t**4 - 104 * t**3 - 831 * t**2 - 140 * t - 6
    require_equal(16 * t**2 * (a * c - 3 * b**2), polynomial)
    require_equal(3 * h2**2 - h1 * h3, -polynomial / (16 * E**2 * t**8))
    shifted = Poly(expand(polynomial.subs(t, s + 30)), s)
    coefficients = shifted.all_coeffs()
    if coefficients != [6, 616, 22209, 317200, 1299894]:
        raise AssertionError(f"incorrect shifted coefficients: {coefficients}")

    signed_jet = simplify((3 * h2**2 - h1 * h3) / h1**5)
    require_equal(signed_jet, E**3 * t**2 * polynomial / (16 * a**5))
    positive_jet = E**3 * (s + 30)**2 * shifted.as_expr() / (16 * a.subs(t, s + 30)**5)
    if positive_jet.is_positive is not True:
        raise AssertionError("the third inverse derivative must be positive for T >= 30")

    r, y = symbols("r y")
    inverse = Function("R")
    require_equal(diff(inverse(f_zero * y), y, 3),
                  f_zero**3 * diff(inverse(r), r, 3).subs(r, f_zero * y))
    require_equal(1 / f_zero + weight / f_zero, 1)
    if (1 / f_zero).is_positive is not True or (weight / f_zero).is_positive is not True:
        raise AssertionError("both normalized weights must be positive")
    print(f"CM family: shifted coefficients {coefficients}, so the third inverse derivative is positive.")
    print("CM normalization: the two weights sum to one and the derivative scales by a positive cube.")


def check_bernstein_family():
    x, a, b, c, d, epsilon = symbols("x a b c d epsilon", positive=True)
    g = a + d * x + b * (1 - exp(-c * x))
    j = x * g
    g_prime = d + b * c * exp(-c * x)
    require_equal(diff(g, x), g_prime)
    require_equal(diff(j, x), g + x * g_prime)
    require_equal(j.subs(x, 0), 0)
    require_equal(diff(j, x, 2), 2 * d + b * c * (2 - c * x) * exp(-c * x))
    second_at_witness = diff(j, x, 2).subs(x, 3 / c)
    require_equal(second_at_witness, 2 * d - b * c / E**3)

    # This substitution ranges over exactly 0 < d < b*c/(2*e**3).
    interior_drift = b * c / (2 * E**3 * (1 + epsilon))
    curvature_margin = b * c * epsilon / (E**3 * (1 + epsilon))
    require_equal(-second_at_witness.subs(d, interior_drift), curvature_margin)
    increasing_jet = (
        a + 6 * d / c + b * (1 - exp(-3)) + 3 * b * exp(-3)
    )
    require_equal(diff(j, x).subs(x, 3 / c), increasing_jet)
    if curvature_margin.is_positive is not True or increasing_jet.is_positive is not True:
        raise AssertionError("the inverse must have a positive second derivative")
    require_equal(
        second_at_witness.subs({b: 1, c: 1, d: Integer(1) / (4 * E**3)}),
        -1 / (2 * E**3),
    )
    print("BF family: the full strict drift interval gives negative product curvature and positive inverse curvature.")


if __name__ == "__main__":
    check_source_laplace_derivatives()
    check_inverse_derivatives()
    check_completely_monotone_family()
    check_bernstein_family()
    print("ALL_CHECKS_OK")

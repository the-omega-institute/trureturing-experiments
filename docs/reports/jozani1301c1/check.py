"""Exact checks for arXiv:1301.4292v1, Conjecture 1."""

from sympy import Rational, expand, integrate, log, oo, simplify, symbols


def collision_pair(density_quantile, rank_factors, quantile):
    """Integrate the squared marginals in the parent quantile coordinate."""
    return tuple(
        expand(integrate(expand(factor**2 * density_quantile), (quantile, 0, 1)))
        for factor in rank_factors
    )


def main():
    u, z = symbols("u z", real=True)
    first, second = 2 * (1 - u), 2 * u
    assert tuple(integrate(factor, (u, 0, 1)) for factor in (first, second)) == (1, 1)

    # For the source's unit-rate exponential, f(F^-1(u)) = 1 - u.
    exp_parent = integrate(1 - u, (u, 0, 1))
    exp_first, exp_second = collision_pair(1 - u, (first, second), u)
    assert (exp_parent, exp_first, exp_second) == (Rational(1, 2), 1, Rational(1, 3))
    # The three displayed individual formulas in Section 3 at alpha = 2.
    assert simplify(-log(exp_parent**2) - 2 * log(2)) == 0
    assert simplify(-log(exp_first) - (-log(2) + log(2))) == 0
    assert simplify(-log(exp_second) - (-2 * log(2) - log(Rational(1, 12)))) == 0
    print("Fidelity (unit-rate exponential, n=2, alpha=2): SRS product 1/4, RSS product 1/3")

    # For f(x) = (1+x)^-2 on the positive half-line, F = x/(1+x).
    x = symbols("x", nonnegative=True)
    f = 1 / (1 + x) ** 2
    cdf = integrate(f, x) + 1
    assert simplify(cdf - x / (1 + x)) == 0
    assert cdf.subs(x, 0) == 0 and cdf.limit(x, oo) == 1
    pareto_quantile_density = (1 - u) ** 2
    parent = integrate(pareto_quantile_density, (u, 0, 1))
    a, b = collision_pair(pareto_quantile_density, (first, second), u)
    cross = integrate(first * second * pareto_quantile_density, (u, 0, 1))
    assert (parent, a, b, cross) == (
        Rational(1, 3), Rational(4, 5), Rational(2, 15), Rational(1, 5)
    )
    print(f"Pareto integrals: parent {parent}, rank 1 {a}, rank 2 {b}, cross {cross}")

    p = (1 + z) / 2
    mixed = (p * first + (1 - p) * second, (1 - p) * first + p * second)
    assert tuple(expand(integrate(factor, (u, 0, 1))) for factor in mixed) == (1, 1)
    mixed_a, mixed_b = collision_pair(pareto_quantile_density, mixed, u)
    assert expand(15 * mixed_a - (2 * z**2 + 5 * z + 5)) == 0
    assert expand(15 * mixed_b - (2 * z**2 - 5 * z + 5)) == 0
    perfect_product = a * b
    mixed_product = expand(mixed_a * mixed_b)
    srs_product = parent**2
    assert expand(
        225 * (mixed_product - perfect_product) - (1 - z) * (1 + z) * (1 - 2 * z) * (1 + 2 * z)
    ) == 0
    assert expand(
        225 * (srs_product - mixed_product) - z**2 * (4 + (1 - 2 * z) * (1 + 2 * z))
    ) == 0
    # Every factor on the right of both identities is positive for 0 < z < 1/2.
    print("Exact positive factors certify perfect < imperfect < SRS in collision products")
    print("The interval is 0 < z < 1/2, equivalently 1/2 < p < 3/4")

    products_at_two_thirds = (perfect_product, mixed_product.subs(z, Rational(1, 3)), srs_product)
    assert tuple(18225 * product for product in products_at_two_thirds) == (1944, 1984, 2025)
    assert 0 < products_at_two_thirds[0] < products_at_two_thirds[1] < products_at_two_thirds[2]
    print("p=2/3: perfect, imperfect and SRS collision products are (1944, 1984, 2025)/18225")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

"""Check the Gaussian factor identity and four Taylor residues using exact algebra."""

from math import factorial

import sympy as sp


def check_gaussian_factor():
    """Check the differential identity and the normalization in (4.1)."""
    y = sp.Symbol("y", positive=True)
    # Put y = exp(z). The derivative of I(y) with respect to z is log(1 + y).
    log_g_derivative = 1 - sp.log(1 + y)
    reflected_derivative = -log_g_derivative.subs(y, 1 / y)
    log_product_derivative = sp.logcombine(
        log_g_derivative + reflected_derivative
    )
    assert sp.simplify(log_product_derivative + sp.log(y)) == 0

    # The first expression of (4.2) has zero integrand at z = 0.
    g_at_zero = (2 / sp.pi) ** sp.Rational(1, 4)
    phi_at_zero = 1 / sp.sqrt(2 * sp.pi)
    assert sp.simplify(g_at_zero**2 / 2 - phi_at_zero) == 0

    curvature = sp.simplify(y * sp.diff(log_g_derivative, y))
    assert curvature == -y / (1 + y)
    print("Fidelity (4.1): log derivative = -log(y), normalization agrees at z=0")
    print(f"Fidelity (4.2): log curvature = {curvature} for y=exp(z)>0")


def check_residue_obstruction():
    """Generate coefficients by F' = (log F)' F, then compare with the route."""
    a = sp.Symbol("a", positive=True)
    coefficients = [sp.S.One]
    for degree in range(1, 4):
        coefficient = sum(
            sp.Rational(1, power) * a * coefficients[degree - power]
            for power in range(1, degree + 1)
        ) / degree
        coefficients.append(sp.expand(coefficient))

    z = sp.Symbol("z")
    truncated_exponent = a * sum(
        z**power / power**2 for power in range(1, 4)
    )
    expanded = sp.Poly(
        sp.series(sp.exp(truncated_exponent), z, 0, 4).removeO(), z
    )
    assert coefficients == [expanded.nth(k) for k in range(4)]

    scaled = [factorial(k) * coefficient for k, coefficient in enumerate(coefficients)]
    b0, b1, b2, b3 = scaled
    obstruction = sp.factor(b0 * b2**3 - b1**3 * b3)
    assert obstruction == a**3 * (2 * a + 3) / 24
    assert obstruction.is_positive is True

    source_residues = [
        (-1) ** k * coefficient.subs(a, 1)
        for k, coefficient in enumerate(coefficients)
    ]
    assert source_residues == [1, -1, sp.Rational(3, 4), -sp.Rational(19, 36)]
    assert obstruction.subs(a, 1) == sp.Rational(5, 24)
    print(f"Normalized signed residues at a=b=1: {source_residues}")
    print(f"Factorial-scaled coefficients: {scaled}")
    print(f"PF-infinity identity defect: {obstruction} > 0 for a>0")
    print(f"Defect at a=1: {obstruction.subs(a, 1)}")


if __name__ == "__main__":
    check_gaussian_factor()
    check_residue_obstruction()
    print("ALL_CHECKS_OK")

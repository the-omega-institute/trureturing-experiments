#!/usr/bin/env python3
"""Exact finite checks for Jeffery's OEIS A052750 square-sum product.

The cyclotomic computation uses t = exp(pi*i/(2*n+1)) and hence
x(n,k) = t**(2*k-1) + t**(-(2*k-1)). All reductions are in QQ[t]
modulo the exact cyclotomic polynomial Phi_{2*(2*n+1)}. No floating
point evaluation or numerical root isolation is used.

A separate resultant computation multiplies the square-sum polynomial
over the roots of S_n - S_{n-1}. This is a cross-check of the finite
computations, not a proof for unbounded n.
"""
import sys

import sympy as sp


X = sp.Symbol("X")
T = sp.Symbol("T")


def s_polynomials(last):
    """The literal recurrence S_0=1, S_1=X, S_{r+2}=X*S_{r+1}-S_r."""
    values = [sp.Poly(1, X, domain=sp.ZZ), sp.Poly(X, X, domain=sp.ZZ)]
    for r in range(2, last + 1):
        values.append(sp.Poly(X * values[-1].as_expr() - values[-2].as_expr(), X))
    return values[:last + 1]


def fidelity():
    # G. C. Greubel's A052750 b-file, entries 0 through 12.
    published = [1, 1, 5, 49, 729, 14641, 371293, 11390625, 410338673,
                 16983563041, 794280046581, 41426511213649, 2384185791015625]
    computed = [1] + [(2*n + 1)**(n - 1) for n in range(1, 13)]
    assert computed == published
    # Lang's explicit n=4 example in the source entry, independent of the conjecture.
    s = s_polynomials(4)
    p4 = s[4] - s[3]
    assert p4 == sp.Poly(X**4 - X**3 - 3*X**2 + 2*X + 1, X)
    assert sp.discriminant(p4.as_expr(), X) == 729
    print("FIDELITY b052750[0..12]=" + ",".join(map(str, computed)))
    print("FIDELITY P4=X^4-X^3-3*X^2+2*X+1; discriminant=729")


def check_n(n):
    order = 2 * (2*n + 1)
    modulus = sp.Poly(sp.cyclotomic_poly(order, T), T, domain=sp.QQ)
    one = sp.Poly(1, T, domain=sp.QQ)
    zero = sp.Poly(0, T, domain=sp.QQ)
    product = one
    for k in range(1, n + 1):
        exponent = 2*k - 1
        xk = sp.Poly(T**exponent + T**(order - exponent), T, domain=sp.QQ) % modulus
        previous, current = one, xk
        factor = one  # m=0, S_0^2
        for r in range(2, 2*n - 1):
            previous, current = current, (xk*current - previous) % modulus
            if r % 2 == 0:
                factor = (factor + current*current) % modulus
        # Independent check of the proposed one-factor bridge, in this exact field.
        # 4*sin(theta)^2 = 4 - xk^2; thus no trigonometric approximation is needed.
        assert ((4*one - xk*xk)*factor - (2*n + 1)*one) % modulus == zero
        product = product*factor % modulus
    expected = (2*n + 1)**(n - 1)
    assert product == expected*one

    # Independent integer-polynomial route: roots of S_n-S_{n-1} are
    # exactly the source's n odd-angle nodes, and the polynomial is monic.
    s = s_polynomials(max(n, 2*n - 2))
    roots_polynomial = s[n] - s[n - 1]
    square_sum = sum((s[2*m]**2 for m in range(n)), sp.Poly(0, X, domain=sp.ZZ))
    resultant = sp.resultant(roots_polynomial.as_expr(), square_sum.as_expr(), X)
    assert resultant == expected
    print(f"n={n:2d} field_degree={modulus.degree():2d} "
          f"cyclotomic_product={expected} resultant={resultant} "
          f"individual_factors={n}:OK")


def main():
    print(f"Python {sys.version.split()[0]}; SymPy {sp.__version__}; arithmetic=exact ZZ/QQ")
    fidelity()
    for n in range(1, 21):
        check_n(n)
    print("CHECKED target n=1..20; no conclusion about untested n follows from these checks")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

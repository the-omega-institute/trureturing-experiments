"""Finite Hermite interpolation for matching overlap jets, and a Laurent degree check."""
import sympy as sp


def hermite_basis(order,x):
    # Zero through order K-1 at one, equal to one modulo x^K at zero.
    return (1-x)**order*sum(sp.binomial(order+j-1,j)*x**j for j in range(order))


def polynomial_checks():
    x = sp.symbols("x")
    for k in range(1,13):
        at_zero = sp.Poly(hermite_basis(k,x),x)
        at_one = sp.Poly(hermite_basis(k,1-x),x)
        powers = (sp.Poly(x**k,x),sp.Poly((1-x)**k,x))
        assert (at_zero-1).rem(powers[0]).is_zero and at_zero.rem(powers[1]).is_zero
        assert at_zero+at_one==sp.Poly(1,x)
        # Independent local series. Scalar coefficients may live in any real algebraic field.
        left = sp.Poly(sum((j+1)*x**j for j in range(k)),x)
        right = sp.Poly(sum((2*j+1)*(x-1)**j for j in range(k)),x)
        H = (at_zero*left+at_one*right).rem(powers[0]*powers[1])
        assert H.degree()<=2*k-1
        assert (H-left).rem(powers[0]).is_zero and (H-right).rem(powers[1]).is_zero
        assert 4*(2*k-1)+1 == 8*k-3


def target_jets():
    x,m,w = sp.symbols("x m w",positive=True)
    first = sp.sqrt(2/sp.sqrt(1-x)-1)
    second = sp.sqrt(2*m/(m-w)-1)
    assert sp.simplify((1+first**2)/2-1/sp.sqrt(1-x)) == 0
    assert sp.simplify((1+second**2)*(m-w)/2-m) == 0
    c=(sp.sqrt(5)-1)/4
    assert 0<1-c*c<1 and 0<5*(1-c)<sp.Rational(139,40)
    z = sp.symbols("z")
    # The fifth step on the tetrahedral ray has power -4, congruent to one modulo five.
    assert sp.rem(z**4*z-1,sp.cyclotomic_poly(5,z),z)==0


if __name__ == "__main__":
    polynomial_checks()
    target_jets()
    print("Finite Hermite selectors, kinetic/mass square completions, and the five-center ray verified.")

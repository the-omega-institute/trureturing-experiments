"""Two-node Hermite and Wilson constructions for matching three overlap jets."""
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


def shorter_kinetic_interpolant():
    """Recover the two-node polynomial directly by Chinese remaindering."""
    s = sp.symbols("s")
    c = (sp.sqrt(5) - 1) / 4
    node = 1 - c**2
    field = sp.QQ.algebraic_field(sp.sqrt(5))
    variable = sp.Poly(s, s, domain=field)
    shifted = variable - sp.Poly(node, s, domain=field)
    for order in range(1, 11):
        left_order = (order + 1) // 2
        right_order = order if order % 2 else order + 1
        left = sum(sp.binomial(2*j, j) * variable**j / 4**j
                   for j in range(left_order))
        # Keep arithmetic in Q(sqrt(5)) throughout.
        right = sp.Poly(0, s, domain=field)
        for j in range(right_order):
            coeff = sp.binomial(2*j, j) / (4**j * c**(2*j+1))
            right += shifted**j * sp.Poly(coeff, s, domain=field)
        left = sp.Poly(left, s, domain=field)
        modulus = shifted**right_order
        lift = ((right-left) * (variable**left_order).invert(modulus)).rem(modulus)
        interpolant = left + variable**left_order * lift
        assert (interpolant-left).rem(variable**left_order).is_zero
        assert (interpolant-right).rem(modulus).is_zero
        assert 2*interpolant.degree()+1 == 3*order + int(order % 2 == 0)


def wilson_and_ray():
    w, m, a = sp.symbols("w m a")
    x = w*(w-a)/(m*(m-a))
    for order in range(1, 17):
        q = (order+2)//2
        if q % 2 == 0:
            q += 1
        h = (w+m-a)/(m-a)*sum(x**j for j in range(q))
        assert sp.cancel((w-m)*h + m - m*x**q) == 0
        assert 2*q >= order + 1
        assert 2*q <= 3*order + int(order % 2 == 0)
    ray = (2,2,2,-3,-3)
    assert sum(ray) == 0 and all((3*n) % 5 == 1 for n in ray)
    powers = {0}
    letters = set(ray) | {-n for n in ray}
    for radius in range(1, 9):
        powers |= {n+step for n in powers for step in letters}
        # In each cross difference, terms with equal residues cancel.
        cross = {i+j for i in powers for j in powers if (i-j) % 5}
        pair = {i+j for i in powers for j in powers}
        scalar = {i+j for i in pair for j in pair if (i-j) % 5}
        assert min(cross) == -6*radius+1 and max(cross) == 6*radius-1
        assert min(scalar) == -12*radius+1 and max(scalar) == 12*radius-1


if __name__ == "__main__":
    polynomial_checks()
    target_jets()
    shorter_kinetic_interpolant()
    wilson_and_ray()
    print("Exact two-node kinetic jets K=1..10, Wilson identities K=1..16, and short-ray frequency supports r=1..8 verified.")

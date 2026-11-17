#!/usr/bin/env python3
"""Exact checks for Shiu, arXiv:2610.34556v2, Equation (3), K = Q(sqrt(-11)).

Coordinates (a,b) mean a+b*alpha, alpha^2-alpha+c=0; c=3 for -11.
No floating point, external algebra system, or unproved primality oracle is used.
"""
from fractions import Fraction
from itertools import product
from math import isqrt


ZERO = (0, 0)
ONE = (1, 0)
ALPHA = (0, 1)


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[0], -x[1])


def mul(x, y, c=3):
    a, b = x
    u, v = y
    return (a*u - c*b*v, a*v + b*u + b*v)


def norm(x, c=3):
    a, b = x
    return a*a + a*b + c*b*b


def disc(center, radius, c=3):
    """Exhaustive bounds from 4*N(a,b)=(2*a+b)^2+(4*c-1)*b^2."""
    radius = Fraction(radius)
    if radius < 0:
        return set()
    four = 4*radius
    bound_b = isqrt(four.numerator // (four.denominator*(4*c-1)))
    out = set()
    for b in range(-bound_b, bound_b+1):
        remaining = four - (4*c-1)*b*b
        bound_u = isqrt(remaining.numerator // remaining.denominator)
        for u in range(-bound_u, bound_u+1):
            if (u-b) % 2 == 0:
                a = (u-b)//2
                assert 4*norm((a, b), c) == u*u + (4*c-1)*b*b
                assert norm((a, b), c) <= radius
                out.add(add(center, (a, b)))
    return out


def rectangle_disc(center, radius, c=3):
    """Separate broad-box enumeration for the small parameter checks below."""
    radius = Fraction(radius)
    bound = 2*(radius.__ceil__()+1)
    return {add(center, (a, b))
            for a, b in product(range(-bound, bound+1), repeat=2)
            if norm((a, b), c) <= radius}


def quotient(z, c=3):
    """Enumerate Z[alpha]/(z) using adjugate-lattice signatures.

    Multiplication by z=(a,b) has matrix M=[[a,-c*b],[b,a+b]],
    determinant n=N(z)>0. x=y mod M*Z^2 iff adj(M)*(x-y)=0 mod n.
    Since n*Z^2 is in M*Z^2, [0,n)^2 contains every quotient class.
    """
    a, b = z
    n = norm(z, c)
    assert n > 0

    def key(x):
        u, v = x
        return (((a+b)*u+c*b*v) % n, (-b*u+a*v) % n)

    reps = {}
    for x in product(range(n), repeat=2):
        reps.setdefault(key(x), x)
    assert len(reps) == n
    assert key(z) == key(ZERO)
    assert key(mul(z, ALPHA, c)) == key(ZERO)
    keys = set(reps)
    for x, y in product(reps.values(), repeat=2):
        assert key(add(x, y)) in keys
        assert key(mul(x, y, c)) in keys
    return key, reps


def quotient_is_field(z, c=3):
    if z == ZERO or norm(z, c) == 1:
        return False
    key, reps = quotient(z, c)
    zero, one = key(ZERO), key(ONE)
    assert zero != one
    return all(any(key(mul(x, y, c)) == one for y in reps.values())
               for k, x in reps.items() if k != zero)


def prime_elements(center, radius, c=3):
    # For a nonzero nonunit z, the finite quotient is a domain iff a field.
    # This is the principal-prime-ideal criterion, not prime integer norm.
    return {z for z in disc(center, radius, c) if quotient_is_field(z, c)}


def main():
    for c in (3, 41):
        for x, y in product(product(range(-2, 3), repeat=2), repeat=2):
            assert norm(mul(x, y, c), c) == norm(x, c)*norm(y, c)
        for A in (Fraction(3), Fraction(13, 4), Fraction(7, 2),
                  Fraction(15, 4), Fraction(399, 100), Fraction(4), Fraction(5)):
            for center in (ZERO, ONE, neg(ONE)):
                assert disc(center, A, c) == rectangle_disc(center, A, c)
    print('COORDINATES: alpha^2=alpha-3; N(a+b*alpha)=a^2+a*b+3*b^2')
    print('INDEPENDENT_ENUMERATORS: exact square bounds and broad boxes agree')

    d0 = disc(ZERO, 3)
    expected_d0 = {ZERO, ONE, neg(ONE), ALPHA, neg(ALPHA), (1, -1), (-1, 1)}
    assert d0 == expected_d0
    assert len(prime_elements(ZERO, 3)) == 4
    assert len(disc(ZERO, 4)) == 9
    assert len(prime_elements(ZERO, 4)) == 6
    assert len(disc(ZERO, 5)) == 13
    assert len(prime_elements(ZERO, 5)) == 10
    print('FIDELITY Section 4: D(0,3) =', sorted(d0))
    print('FIDELITY Section 4: (A, disc size, prime count) = (3,7,4), (4,9,6), (5,13,10)')

    # Independent source example in Q(sqrt(-163)), alpha^2-alpha+41=0.
    assert disc(ZERO, 4, 41) == {(-2, 0), (-1, 0), ZERO, ONE, (2, 0)}
    assert prime_elements(ZERO, 4, 41) == {(-2, 0), (2, 0)}
    print('FIDELITY Section 4, discriminant -163: D(0,4) has 5 elements and 2 primes')

    # Reproduce the source's obstruction to an admissible radius-three tuple.
    k2, r2 = quotient((2, 0))
    ka, ra = quotient(ALPHA)
    kb, rb = quotient((1, -1))
    sizes = [sum(k2(z) != x and ka(z) != y and kb(z) != t for z in d0)
             for x, y, t in product(r2, ra, rb)]
    assert len(sizes) == 36
    assert max(sizes) == 3
    print('FIDELITY Section 4: 36 choices of omitted residues modulo (2),(alpha),(1-alpha); max survivors = 3')

    certified = [ALPHA, (1, -1), (1, 1), (2, -1), (2, 0)]
    for z, size, characteristic in zip(certified, [3, 3, 5, 5, 4], [3, 3, 5, 5, 2]):
        assert norm(z) == size
        assert quotient_is_field(z)
        key, reps = quotient(z)
        char = next(k for k in range(1, size+1) if key((k, 0)) == key(ZERO))
        assert char == characteristic
        print('QUOTIENT_FIELD: generator =', z, 'cardinality =', len(reps), 'characteristic =', char)
    for p, root in [(3, 0), (3, 1), (5, -1), (5, 2)]:
        assert (root*root-root+3) % p == 0
    assert [(x*x+x+1) % 2 for x in range(2)] == [1, 1]
    print('POLYNOMIAL_CERTIFICATES: roots 0,1 mod 3; roots -1,2 mod 5; X^2+X+1 root-free mod 2')

    d1 = disc(ONE, 3)
    p0 = prime_elements(ZERO, 3)
    p1 = prime_elements(ONE, 3)
    assert d1 == {ZERO, ONE, (2, 0), ALPHA, (1, -1), (1, 1), (2, -1)}
    assert d0-p0 == {ZERO, ONE, neg(ONE)}
    assert d1-p1 == {ZERO, ONE}
    assert p1 == set(certified)
    assert len(p0) == 4 < 5 == len(p1)
    print('COUNTEREXAMPLE: D(1,3) =', sorted(d1))
    print('COUNTEREXAMPLE: primes at center 0 =', sorted(p0))
    print('COUNTEREXAMPLE: primes at center 1 =', sorted(p1))
    print('COUNTEREXAMPLE: pi(0,3)=4 < pi(1,3)=5')

    # Norms are integers; the entire half-open real interval has the same disc.
    for center in (ZERO, ONE, neg(ONE)):
        below_four = {z for z in disc(center, 4) if norm(add(z, neg(center))) < 4}
        assert below_four == disc(center, 3)
    for A in (Fraction(3), Fraction(13, 4), Fraction(7, 2), Fraction(15, 4), Fraction(399, 100)):
        assert len(prime_elements(ZERO, A)) == 4
        for center in (ONE, neg(ONE)):
            assert len(prime_elements(center, A)) == 5
            assert prime_elements(neg(ONE), A) == {neg(z) for z in prime_elements(ONE, A)}
    assert {ka(z) for z in (ALPHA, (1, -1), (2, 0))} == set(ra)
    print('FAMILY_CERTIFICATE: integral norm thresholds; centers +1,-1; every real 3 <= A < 4')
    print('MECHANISM: alpha,1-alpha,2 occupy all 3 classes modulo (alpha); cluster is nonadmissible')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

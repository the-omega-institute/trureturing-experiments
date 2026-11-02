#!/usr/bin/env python3
"""Exact source-fidelity and saturation-counterexample checks; no Lean proof."""
from fractions import Fraction
from itertools import product
import sympy as sp


def zero(expr):
    result = sp.Poly(sp.expand(expr), *sorted(expr.free_symbols, key=str)) if expr.free_symbols else sp.expand(expr)
    assert result == 0 or getattr(result, 'is_zero', False), expr


def source_fidelity():
    # Lido--Mauri, arXiv:2605.03683v1, Example 3.2.3. The symbols g1,g2
    # stand for arbitrary series: a polynomial identity allows substitution.
    p, x, y, g1, g2 = sp.symbols('p x y g1 g2')
    F = (x*x+p*y+p*p*g1, x*y+p*p*g2, y*y+p*y*g1-p*x*g2)
    zero(y*F[0]-x*F[1]-p*F[2])
    syzygies = ((y, -x, -p),
                (y*y-p*x*g2, p*x*g1, -x*x-p*y),
                (-p*g2, y+p*g1, -x))
    for s in syzygies:
        zero(sum(c*f for c, f in zip(s, F)))
    assert [sp.expand(f.subs(p, 0)) for f in F] == [x*x, x*y, y*y]
    print('FIDELITY Example 3.2.3: one colon certificate, three lifted syzygies, reduction (x^2,xy,y^2): OK')


# B = Q[alpha,beta]/(alpha^2,beta^2), with coordinate order 1,alpha,beta,alpha*beta.
BASIS = ((0, 0), (1, 0), (0, 1), (1, 1))
ZERO = (Fraction(0),)*4
ONE = (Fraction(1), Fraction(0), Fraction(0), Fraction(0))
ALPHA = (Fraction(0), Fraction(1), Fraction(0), Fraction(0))
BETA = (Fraction(0), Fraction(0), Fraction(1), Fraction(0))
AB = (Fraction(0), Fraction(0), Fraction(0), Fraction(1))


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(c, a):
    return tuple(Fraction(c)*x for x in a)


def mul(a, b):
    out = [Fraction(0)]*4
    for i, (ia, ib) in enumerate(BASIS):
        for j, (ja, jb) in enumerate(BASIS):
            if ia+ja < 2 and ib+jb < 2:
                out[BASIS.index((ia+ja, ib+jb))] += a[i]*b[j]
    return tuple(out)


def power(a, n):
    out = ONE
    for _ in range(n):
        out = mul(out, a)
    return out


def evaluate(poly, variables, images):
    out = ZERO
    for monomial, coefficient in sp.Poly(poly, *variables, domain=sp.QQ).terms():
        term = scale(Fraction(coefficient), ONE)
        for image, exponent in zip(images, monomial):
            term = mul(term, power(image, exponent))
        out = add(out, term)
    return out


def counterexample():
    t, y, z, w = sp.symbols('t y z w')
    x, a = t-y, w-1
    f, g, h = x*x*a, y*y, x*z+y*w
    H = t*(w+2*z-1)+2*y
    identities = (
        t*H-(f+2*(t+y)*h+(2*z-3*w+1)*g),
        t*t*z*z-((h+2*y*(z-w))*h+(z-w)**2*g),
        t*t*w*a-((w-2*z)*f+2*a*x*h+a*w*g),
        t*t*z*a-(z*f+2*a*y*h+a*(z-2*w)*g),
        t**3*a-((t+2*y)*f+(3*t-2*y)*a*g),
    )
    for identity in identities:
        zero(identity)
    print('COUNTEREXAMPLE: five membership-certificate residuals = [0, 0, 0, 0, 0] over Z[t,y,z,w]')
    assert sp.expand(H.subs(t, 0)) == 2*y
    # These four elements are in I_2; at the origin their cotangent vectors
    # are independent. Together with I_2 contained in the maximal ideal,
    # this supplies the generator-bound obstruction after localization.
    first_order = sp.Matrix([[sp.diff(q, v).subs({t: 0, y: 0, z: 0, w: 0})
                              for v in (t, y, z, w)]
                             for q in (t*a, z*a, w*a, H)])
    assert first_order.det() in (2, -2)
    print('GENERATOR BOUND: four elements t*a,z*a,w*a,H in I_2 have cotangent determinant', first_order.det(), 'at the origin')
    J = sp.groebner((y, z*z, (w-1)*z, w*(w-1)), y, z, w, domain=sp.QQ)
    assert J.is_zero_dimensional
    # Verify the stated rewrite relations and all products of the proposed basis.
    assert [sp.expand(J.reduce(p)[1]) for p in (y, z*z, w*w-w, w*z-z)] == [0]*4
    standard = [sp.Integer(1), z, w]
    for b, c in product(standard, repeat=2):
        remainder = sp.Poly(J.reduce(b*c)[1], y, z, w, domain=sp.QQ)
        assert set(remainder.monoms()) <= {(0, 0, 0), (0, 1, 0), (0, 0, 1)}
    assert J.reduce(sp.Integer(1))[1] == 1
    print('SPECIAL FIBER: Q[y,z,w]/J has standard monomials [1,z,w]; J is proper and zero-dimensional')
    images = (add(ALPHA, BETA), BETA, ZERO, ZERO)
    phi = lambda p: evaluate(p, (t, y, z, w), images)
    assert mul(ALPHA, ALPHA) == mul(BETA, BETA) == ZERO
    assert mul(ALPHA, BETA) == AB != ZERO
    assert [phi(p) for p in (f, g, h)] == [ZERO]*3
    nil_generators = (ALPHA, BETA)
    assert all(mul(mul(u, v), q) == ZERO for u, v, q in product(nil_generators, repeat=3))
    assert phi(t*t) == scale(2, AB)
    assert phi(t*t*a) == scale(-2, AB) != ZERO
    assert phi(t**3*a) == ZERO
    # Multiplication by phi(t^2) retains exactly the constant coordinate.
    for i, e in enumerate((ONE, ALPHA, BETA, AB)):
        assert mul(phi(t*t), e) == (scale(2, AB) if i == 0 else ZERO)
    print('OBSTRUCTION: phi(f)=phi(g)=phi(h)=0; phi(t^2*(w-1))=(0,0,0,-2) != 0; phi(t^3*(w-1))=0')
    for p in (3, 5, 7, 11, 13):
        assert (-2) % p != 0
    print('FAMILY: identities are integral; obstruction coefficient -2 is nonzero in every field of characteristic != 2 (sample primes 3,5,7,11,13 checked)')
    assert [sp.expand(F.subs({y: 0, z: 0, w: 1})) for F in (f, g, h)] == [0]*3
    print('NONDEGENERACY: all three generators vanish at (y,z,w)=(0,0,1) over Q[[t]] and Q((t))')
    print('BOUNDARY: exact polynomial/finite-algebra checks support the prose argument; restricted-series and Krull-dimension bridges are not kernel-checked')


def main():
    source_fidelity()
    counterexample()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

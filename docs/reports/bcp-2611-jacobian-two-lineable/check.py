#!/usr/bin/env python3
"""Exact checks for arXiv:2611.03572v1, Question 2; no Lean proof.

Uses the coordinate formulas printed in the paper, not third-party program code.
Requires Python 3 and SymPy. All coefficients, points and comparisons are exact.
"""
import sys
import sympy as sp


def equal(left, right):
    assert sp.expand(left - right) == 0, (left, right)


def matrix_equal(left, right):
    assert left.shape == right.shape
    for x, y in zip(left, right):
        equal(x, y)


def at(vector, variables, point):
    return vector.subs(dict(zip(variables, point)), simultaneous=True)


def complex_structure(m):
    block = sp.Matrix([[0, -1], [1, 0]])
    return sp.diag(*([block] * m))


def main():
    print('Python', sys.version.split()[0], 'SymPy', sp.__version__)
    x, y, z = sp.symbols('x y z', real=True)
    a, b = sp.symbols('a b', real=True)
    F0 = sp.Matrix([
        (1 + x*y)**3*z + y**2*(1 + x*y)*(4 + 3*x*y),
        y + 3*x*(1 + x*y)**2*z + 3*x*y**2*(4 + 3*x*y),
        2*x - 3*x**2*y - x**3*z,
    ])
    assert not F0.atoms(sp.Float)
    DF0 = F0.jacobian((x, y, z))
    equal(DF0.det(method='domain-ge'), -2)
    print('FIDELITY: det(D F0) = -2 as an exact polynomial identity')

    # All three preimages printed in the source, independently of the proposed pair.
    paper_points = [
        (0, 0, -sp.Rational(1, 4)),
        (1, -sp.Rational(3, 2), sp.Rational(13, 2)),
        (-1, sp.Rational(3, 2), sp.Rational(13, 2)),
    ]
    target = sp.Matrix([-sp.Rational(1, 4), 0, 0])
    assert len(set(paper_points)) == 3
    for point in paper_points:
        matrix_equal(at(F0, (x, y, z), point), target)
    print('FIDELITY: the three source preimages all map to (-1/4, 0, 0)')

    u = (2, -sp.Rational(1, 2), sp.Rational(5, 4))
    v = (-2, 1, 2)
    assert u != v
    for point in (u, v):
        matrix_equal(at(F0, (x, y, z), point), sp.Matrix([0, 1, 0]))
    print('SETTLEMENT: (2,-1/2,5/4) and (-2,1,2) both map to (0,1,0)')

    rational_grid = [sp.Rational(-2), sp.Rational(-1), -sp.Rational(1, 2),
                     sp.Rational(0), sp.Rational(1, 2), sp.Rational(1), sp.Rational(2)]
    tested = 0
    for m in range(1, 7):
        N = 2*m
        I = sp.eye(N)
        J = complex_structure(m)
        A = a*I + b*J
        matrix_equal(J*J, -I)
        equal(A.det(method='domain-ge'), (a*a + b*b)**m)
        matrix_equal((a*I - b*J)*A, (a*a + b*b)*I)
        assert sp.Poly(A.det(), a, b).coeff_monomial(a**(2*m)) == 1
        for aa in rational_grid:
            for bb in rational_grid:
                if aa == 0 and bb == 0:
                    continue
                determinant = (aa*aa + bb*bb)**m
                assert determinant > 0
                tested += 1
        print(f'PENCIL: N={N}; J^2=-I; det(aI+bJ)=(a^2+b^2)^{m}; inverse numerator checked')
        if N < 4:
            continue

        # Identity padding uses actual variables; the collided points agree there.
        tail = sp.symbols('t0:'+str(N-3), real=True)
        variables = (x, y, z) + tail
        F = sp.Matrix(list(F0) + list(tail))
        DF = F.jacobian(variables)
        equal(DF.det(method='domain-ge'), -2)
        padded_u = u + (sp.Integer(0),)*(N-3)
        padded_v = v + (sp.Integer(0),)*(N-3)
        assert padded_u != padded_v
        fu = at(F, variables, padded_u)
        fv = at(F, variables, padded_v)
        matrix_equal(fu, fv)

        G = A*F
        DG = G.jacobian(variables)
        matrix_equal(DG, A*DF)
        # Determinant multiplicativity combines two exact determinant identities.
        equal(A.det(method='domain-ge')*DF.det(method='domain-ge'),
              -2*(a*a + b*b)**m)
        matrix_equal(at(G, variables, padded_u), at(G, variables, padded_v))
        # At the fixed origin DF is nonsingular: no coefficient pair can give
        # the zero polynomial map unless its matrix pencil is zero.
        derivative_at_origin = at(DF, variables, (sp.Integer(0),)*N)
        equal(derivative_at_origin.det(), -2)
        matrix_equal((A*derivative_at_origin)*derivative_at_origin.inv(), A)
        assert A[0, 0] == a and A[1, 0] == b
        if N == 4:
            # Independent direct determinant calculation of the transformed map.
            equal(DG.det(method='domain-ge'), -2*(a*a + b*b)**m)
        assert not G.atoms(sp.Float)
        print(f'FAMILY: N={N}; padded det=-2; DG=(aI+bJ)DF; symbolic collision preserved')
    print(f'POSITIVITY: {tested} exact rational nonzero coefficient pairs checked')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

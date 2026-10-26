#!/usr/bin/env python3
"""Exact diagnostics for Zimmermann--Stoye arXiv:2403.01879v3, Conjecture 1.

Equation (26) is a homogeneous curvature numerator on arbitrary tangent blocks;
for a Euclidean orthonormal pair it is the sectional curvature itself.
No floating point arithmetic is used. Symbolic checks in fixed dimensions do
not prove the dimension-independent conjecture or Ge's commutator inequality.
"""
from itertools import combinations
import sys
import sympy as S


def equal(label, lhs, rhs):
    residual = S.expand(lhs - rhs)
    assert residual == 0, (label, residual)


def frob2(m):
    return sum((x * x for x in m), S.Integer(0))


def curvature(a, c, b, d):
    return (frob2(b * c - d * a)
            + frob2(b * d.T - d * b.T) / 2
            - frob2(b.T * d - d.T * b) / 2
            + frob2(a * c - c * a - (d.T * b - b.T * d)) / 4)


def gauss(a, c, b, d):
    x, y = a.col_join(b), c.col_join(d)
    z = x.T * y
    sym = (z + z.T) / 2
    return S.trace(x.T * x * y.T * y) - frob2(sym)


def symbolic_skew(p, prefix):
    a = S.zeros(p)
    for i, j in combinations(range(p), 2):
        t = S.Symbol(f'{prefix}_{i}_{j}', real=True)
        a[i, j], a[j, i] = t, -t
    return a


def row(m, i):
    return m[i, :]


def dot(a, b):
    return (a * b.T)[0, 0]


def check_onb(a, c, b, d):
    assert a.T == -a and c.T == -c
    equal('unit X', frob2(a) + frob2(b), 1)
    equal('unit Y', frob2(c) + frob2(d), 1)
    equal('orthogonal XY', sum(a.multiply_elementwise(c)) + dot(b, d), 0)


def fidelity():
    # Theorem 10(5): exactly the displayed X,Y on St(3,2).
    a = S.zeros(2)
    c = S.Matrix([[0, -S.sqrt(2) / 2], [S.sqrt(2) / 2, 0]])
    b, d = S.Matrix([[-1, 0]]), S.zeros(1, 2)
    check_onb(a, c, b, d)
    equal('Theorem 10(5), eq.(26)', curvature(a, c, b, d), S.Rational(1, 2))
    equal('Theorem 10(5), Gauss', gauss(a, c, b, d), S.Rational(1, 2))
    # Equation (28), with zero A blocks and the displayed two-row B blocks.
    b = S.Matrix([[0, 1], [0, 0]])
    d = S.Matrix([[1, 0], [0, 0]])
    check_onb(S.zeros(2), S.zeros(2), b, d)
    equal('eq.(28) lower plane', curvature(S.zeros(2), S.zeros(2), b, d), -S.Rational(1, 2))
    print('FIDELITY: Theorem 10(5) K=1/2; eq.(28) K=-1/2 (exact)')


def symbolic_checks():
    # All entries are unrestricted symbols: these are identities, not samples.
    for p in range(2, 6):
        a, c = symbolic_skew(p, 'a'), symbolic_skew(p, 'c')
        b = S.Matrix(1, p, lambda i, j: S.Symbol(f'b_{j}', real=True))
        d = S.Matrix(1, p, lambda i, j: S.Symbol(f'd_{j}', real=True))
        equal(f'Gauss eq.(26), p={p}', curvature(a, c, b, d), gauss(a, c, b, d))
        r, s = S.symbols('r s', real=True)
        b, d = S.zeros(1, p), S.zeros(1, p)
        b[0, 0], d[0, 1] = r, s
        coupling = dot(row(c, 0), row(a, 1)) + dot(row(a, 0), row(c, 1))
        tail = sum((c[0, j] * a[1, j] + a[0, j] * c[1, j] for j in range(2, p)), S.Integer(0))
        equal(f'coupling cancellation, p={p}', coupling, tail)
        expansion = (frob2(a * c - c * a) / 4
                     + r * r * frob2(row(c, 0)) + s * s * frob2(row(a, 1))
                     - r * s * coupling - r * r * s * s / 2)
        equal(f'normalized row expansion, p={p}', curvature(a, c, b, d), expansion)
        av = [a[i, j] for j in range(2, p) for i in (0, 1)]
        cv = [c[i, j] for j in range(2, p) for i in (1, 0)]
        aa = sum((x * x for x in av), S.Integer(0))
        cc = sum((x * x for x in cv), S.Integer(0))
        lagrange = sum(((av[i] * cv[j] - av[j] * cv[i]) ** 2
                        for i, j in combinations(range(len(av)), 2)), S.Integer(0))
        equal(f'Cauchy sum-of-squares, p={p}', aa * cc - tail * tail, lagrange)
        # The half-Frobenius tail and individual-row budgets are sums of squares.
        for m, used in [(a, aa), (c, cc), (a, frob2(row(a, 1))), (c, frob2(row(c, 0)))]:
            poly = S.Poly(S.expand(frob2(m) / 2 - used), *sorted(m.free_symbols, key=str))
            assert all(coeff >= 0 for coeff in poly.coeffs())
        print(f'SYMBOLIC p={p}: Gauss identity, coupling expansion, Cauchy SOS, row budgets OK')
    alpha, gamma, r, s = S.symbols('alpha gamma r s', real=True)
    bound = (alpha**2 * gamma**2 / 4 + (r*r * gamma**2 + s*s * alpha**2) / 2
             + alpha * gamma * r * s / 2 - r*r*s*s / 2)
    equal('two-square deficit', (alpha**2+r*r)*(gamma**2+s*s)/2-bound,
          (r*s-alpha*gamma/4)**2 + 3*alpha**2*gamma**2/16)
    equal('normalization Cauchy', (alpha**2+r*r)*(gamma**2+s*s)-(alpha*gamma+r*s)**2,
          (alpha*s-gamma*r)**2)
    print('SCALAR: two-square upper deficit and normalization Cauchy identity OK')


def endpoint_and_finite_checks():
    total = 0
    for p in range(2, 9):
        z = S.zeros(p)
        e1, e2 = S.zeros(1, p), S.zeros(1, p)
        e1[0, 0], e2[0, 1] = 1, 1
        j = S.zeros(p)
        j[0, 1], j[1, 0] = S.sqrt(2)/2, -S.sqrt(2)/2
        zr = S.zeros(1, p)
        for a, c, b, d, k in [(z, j, e1, zr, S.Rational(1, 2)),
                               (z, z, e1, e2, -S.Rational(1, 2))]:
            check_onb(a, c, b, d)
            equal(f'endpoint p={p}', curvature(a, c, b, d), k)
        # Exhaust every coordinate 2-plane in the normalized block basis for p<=5.
        if p <= 5:
            basis = []
            for i, k in combinations(range(p), 2):
                m = S.zeros(p)
                m[i, k], m[k, i] = S.sqrt(2)/2, -S.sqrt(2)/2
                basis.append((m, zr))
            for i in range(p):
                e = S.zeros(1, p);e[0, i] = 1;basis.append((z, e))
            for (a, b), (c, d) in combinations(basis, 2):
                check_onb(a, c, b, d)
                k = S.simplify(curvature(a, c, b, d))
                assert k.is_Rational and -S.Rational(1, 2) <= k <= S.Rational(1, 2)
                total += 1
    print('ENDPOINTS: exact K=+/-1/2 for p=2,...,8 (n=p+1)')
    print(f'FINITE: all {total} coordinate tangent 2-planes for p=2,...,5 satisfy both bounds')


def main():
    print(f'Python {sys.version_info.major}.{sys.version_info.minor}; SymPy {S.__version__}; exact algebra only')
    fidelity()
    symbolic_checks()
    endpoint_and_finite_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for Zanardi arXiv:2609.28662v2, Conjecture 1.

These checks verify curvature identities and finite-dimensional certificates.
They do not decide density or verify the differential-geometric proof.
"""
from fractions import Fraction
from itertools import combinations
import sympy as sp


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def zero(matrix):
    return all(sp.cancel(entry) == 0 for entry in matrix)


def comm(a, b):
    return a * b - b * a


def xy(d, a, b):
    x, y = sp.zeros(d), sp.zeros(d)
    x[a, b], x[b, a] = 1, -1
    y[a, b] = y[b, a] = sp.I
    return x, y


def compositions(n):
    if n == 0:
        yield ()
    else:
        for first in range(1, n + 1):
            for tail in compositions(n - first):
                yield (first,) + tail


def band_labels(sizes):
    return [band for band, size in enumerate(sizes) for _ in range(size)]


def project(a, labels):
    return a.applyfunc(lambda entry: entry) if len(set(labels)) == 1 else sp.Matrix(
        a.rows, a.cols,
        lambda i, j: a[i, j] if labels[i] == labels[j] else 0,
    )


def real_column(a):
    column = []
    for z in a:
        re, im = sp.expand_complex(z).as_real_imag()
        require(re.is_Integer and im.is_Integer, 'nonintegral rank input')
        column.extend([int(re), int(im)])
    return column


def independent_indices(columns):
    """Exact rational column elimination; return a basis of the input span."""
    pivots, selected = {}, []
    for k, column in enumerate(columns):
        v = list(map(Fraction, column))
        for p, basis in sorted(pivots.items()):
            coefficient = v[p]
            if coefficient:
                v = [x - coefficient * y for x, y in zip(v, basis)]
        p = next((i for i, x in enumerate(v) if x), None)
        if p is not None:
            coefficient = v[p]
            pivots[p] = [x / coefficient for x in v]
            selected.append(k)
    return selected


def source_curvature():
    """Reproduce the source's block-curvature formula on an actual U(x,y)."""
    x, y = sp.symbols('x y', real=True)
    X, Y = xy(2, 0, 1)
    k = (x * X + y * Y) / 2
    unit = sp.eye(2)
    u = ((unit + k) * (unit - k).inv()).applyfunc(sp.cancel)
    require(zero(u.conjugate().T * u - unit), 'Cayley family not unitary')
    tx = (u.conjugate().T * u.diff(x)).applyfunc(sp.cancel)
    ty = (u.conjugate().T * u.diff(y)).applyfunc(sp.cancel)
    labels = [0, 1]
    ax, ay = project(tx, labels), project(ty, labels)
    curvature = (ay.diff(x) - ax.diff(y) + comm(ax, ay)).applyfunc(sp.cancel)
    offx, offy = tx - ax, ty - ay
    expected = -project(comm(offx, offy), labels)
    require(zero(curvature - expected), 'source block-curvature identity')
    require(sp.cancel(sp.trace(curvature)) == 0, 'source sum of band traces')
    require(zero(curvature.subs({x: 0, y: 0}) + comm(X, Y)), 'curvature jet at origin')
    print('FIDELITY: rational unitary two-band family; F=-P[theta_off,x,theta_off,y] and tr(F)=0 exactly')


def span_checks():
    count, pairs_checked, max_rank = 0, 0, 0
    representative_basis = None
    for d in range(2, 7):
        for sizes in compositions(d):
            if len(sizes) < 2:
                continue
            labels = band_labels(sizes)
            generators = []
            for a, b in combinations(range(d), 2):
                if labels[a] != labels[b]:
                    generators.extend(xy(d, a, b))
            curvatures = []
            for x, y in combinations(generators, 2):
                c = -project(comm(x, y), labels)
                require(zero(c + c.conjugate().T), 'curvature not skew Hermitian')
                require(sp.trace(c) == 0, 'curvature not traceless')
                require(zero(c - project(c, labels)), 'curvature not block diagonal')
                curvatures.append(c)
            basis_indices = independent_indices(list(map(real_column, curvatures)))
            expected = sum(size * size for size in sizes) - 1
            require(len(basis_indices) == expected, f'span failed for {sizes}')
            if sizes == (2, 1):
                representative_basis = [curvatures[j] for j in basis_indices]
            count += 1
            pairs_checked += len(curvatures)
            max_rank = max(max_rank, expected)
    require(count == 57, 'composition coverage count')
    print(f'SPAN: all {count} ordered positive band tuples, 2<=d<=6; {pairs_checked} projected commutators; rank=sum(d_mu^2)-1; maximum rank {max_rank}')
    return representative_basis


def generator_identities():
    # Indices 0 and 2 belong to one band, and 1 to another.
    x01, y01 = xy(3, 0, 1)
    x21, y21 = xy(3, 2, 1)
    x02, y02 = xy(3, 0, 2)
    diagonal = sp.diag(2 * sp.I, -2 * sp.I, 0)
    require(zero(comm(x01, y01) - diagonal), 'diagonal generator')
    require(zero(comm(x01, x21) + x02), 'within-band real generator')
    require(zero(comm(x01, y21) - y02), 'within-band imaginary generator')
    delta = sp.symbols('delta', real=True)
    require(zero(-comm(delta * x01, delta * y01) + delta**2 * diagonal), 'quadratic jet scaling')
    print('GENERATORS: three elementary commutator identities and delta^2 curvature scaling hold exactly')


def interpolation_certificate(basis):
    """Finite algebra model of the transported-curvature determinant step.

    Transports below are exact block-unitary Cayley curves, not computed WZ
    transports. The test checks the leading determinant and a rational small
    parameter. Analytic dependence of actual transports remains a proof duty.
    """
    delta = sp.symbols('delta', real=True)
    columns = []
    for j, c in enumerate(basis):
        # This block-diagonal skew-Hermitian matrix preserves the band Lie algebra.
        z = (j + 1) * sp.Matrix([[0, 1, 0], [-1, 0, 0], [0, 0, 0]])
        t = ((sp.eye(3) + delta * z / 2) * (sp.eye(3) - delta * z / 2).inv()).applyfunc(sp.cancel)
        transported = (t.inv() * c * t).applyfunc(sp.cancel)
        # Real coordinates on u(2) + u(1), with total trace zero.
        columns.append(sp.Matrix([sp.im(transported[0, 0]), sp.im(transported[1, 1]),
                                  sp.re(transported[0, 1]), sp.im(transported[0, 1])]))
    determinant = sp.factor(sp.Matrix.hstack(*columns).det())
    require(determinant.subs(delta, 0) != 0, 'zero leading coefficient')
    require(determinant.subs(delta, sp.Rational(1, 100)) != 0, 'rational transported basis failed')
    r = len(basis)
    require(r == 4, 'representative dimension')
    print(f'INTERPOLATION MODEL: bands (2,1), r={r}; scaled curvature determinant = delta^{2*r} * D(delta)')
    print('D(delta) =', determinant)
    print('D(0) =', determinant.subs(delta, 0))
    print('D(1/100) =', determinant.subs(delta, sp.Rational(1, 100)))


def main():
    source_curvature()
    generator_identities()
    interpolation_certificate(span_checks())
    print('LIMIT: finite exact algebra checks; no conclusion about C^2 density follows from these computations alone')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

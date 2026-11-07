#!/usr/bin/env python3
"""Exact small-case hull/volume checks; SciPy is an additional numerical diagnostic.

Source: Lee, Skipper, Speakman, arXiv:2605.01493v1, Sections 1--4.
No general proof is implemented. All asserted geometric evidence uses rationals.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import factorial, gcd, lcm, prod
from random import Random
import sympy as sp


def canonical(row):
    scale = lcm(*(x.denominator for x in row))
    ints = [int(x * scale) for x in row]
    divisor = gcd(*ints)
    assert divisor
    return tuple(x // divisor for x in ints)


def kernel_one(rows):
    """Exact one-dimensional right nullspace, or None at lower affine rank."""
    a = [list(map(F, r)) for r in rows]
    nr, nc = len(a), len(a[0])
    pivots = []
    r = 0
    for c in range(nc):
        k = next((i for i in range(r, nr) if a[i][c]), None)
        if k is None:
            continue
        a[r], a[k] = a[k], a[r]
        z = a[r][c]
        a[r] = [v / z for v in a[r]]
        for i in range(nr):
            if i != r and a[i][c]:
                z = a[i][c]
                a[i] = [v - z * w for v, w in zip(a[i], a[r])]
        pivots.append(c)
        r += 1
        if r == nr:
            break
    if r != nc - 1:
        return None
    free = next(c for c in range(nc) if c not in pivots)
    v = [F(0)] * nc
    v[free] = F(1)
    for i, c in enumerate(pivots):
        v[c] = -a[i][free]
    return canonical(v)


def evaluate(row, point):
    return sum(c * x for c, x in zip(row[:-1], point)) + row[-1]


def all_facets(vertices):
    """Exhaust every d-subset; a supporting plane is oriented as row >= 0."""
    d = len(vertices[0])
    checked = set()
    facets = {}
    for ids in combinations(range(len(vertices)), d):
        row = kernel_one([vertices[i] + (F(1),) for i in ids])
        if row is None:
            continue
        if row in checked or tuple(-x for x in row) in checked:
            continue
        checked.add(row)
        vals = [evaluate(row, v) for v in vertices]
        if min(vals) < 0 < max(vals):
            continue
        if max(vals) <= 0:
            row = tuple(-x for x in row)
            vals = [-x for x in vals]
        assert min(vals) == 0 < max(vals)
        facets[row] = frozenset(i for i, x in enumerate(vals) if x == 0)
    return facets


def exact_volume(vertices, facets):
    """Pulling triangulation using the complete exact facet lattice."""
    whole = frozenset(range(len(vertices)))
    ambient_facets = tuple(facets.values())

    @lru_cache(None)
    def triangulate(face):
        if len(face) == 1:
            return [(min(face),)]
        cuts = set(face & f for f in ambient_facets) - {face, frozenset()}
        maximal = [f for f in cuts if not any(f < g for g in cuts)]
        apex = min(face)
        return [(apex,) + s for f in maximal if apex not in f
                for s in triangulate(f)]

    simplices = triangulate(whole)
    d = len(vertices[0])
    assert simplices and all(len(s) == d + 1 for s in simplices)
    total = F(0)
    for simplex in simplices:
        v0 = vertices[simplex[0]]
        mat = sp.Matrix([[sp.Rational(x - y) for x, y in zip(vertices[i], v0)]
                         for i in simplex[1:]])
        total += F(abs(mat.det())) / factorial(d)
    return total, len(simplices)


def normalized_system(n, alpha, beta):
    """Affine inequalities in (z_1,...,z_m,u,v,y), each >= 0."""
    m = n - 2
    zero = tuple(F(0) for _ in range(n + 2))
    def var(i):
        return tuple(F(j == i) for j in range(n + 2))
    one = var(n + 1)
    z = [sp.Matrix(var(i)) for i in range(m)]
    u, v, y, c = map(sp.Matrix, [var(m), var(m + 1), var(n), one])
    h = sum(z, sp.Matrix(zero)) + (1 - m) * c
    rows = []
    for i in range(m):
        rows += [z[i], c - z[i]]
    rows += [u - alpha*c, c-u, v-beta*c, c-v]
    rows += [y, y-alpha*beta*h, y-beta*(h+u-c), y-alpha*(h+v-c),
             y-(alpha+beta-alpha*beta)*h-beta*(u-c)-alpha*(v-c),
             y-h-u-v+2*c]
    # The dummy z=1 is redundant for m>0 and handles the bilinear m=0 case.
    for w in z + [c]:
        rows += [w-y, alpha*w+u-alpha*c-y, beta*w+v-beta*c-y,
                 alpha*beta*w+u+alpha*v-(alpha+alpha*beta)*c-y,
                 alpha*beta*w+beta*u+v-(beta+alpha*beta)*c-y]
    return {canonical(tuple(map(F, row))) for row in rows}


def proposed_volume(n, lower, upper):
    alpha, beta = lower[-2]/upper[-2], lower[-1]/upper[-1]
    p = prod(upper)
    bracket = ((factorial(n)-1) + (factorial(n-1)-n)*(alpha+beta)
               + (2*factorial(n-2)-n*n+n+1)*alpha*beta)
    return (1-alpha)*(1-beta)*p*p*bracket/factorial(n+1)


def source_c1_volume(n, lower, upper):
    return ((upper[-1]-lower[-1])*prod(b*b for b in upper[:-1])
            * ((factorial(n)-1)*upper[-1]
               + (factorial(n-1)-n)*lower[-1])/factorial(n+1))


def source_c1_rows(n, alpha):
    """Section 2 Theorem 1, normalized, with nonnegative zero-lower variables."""
    # Obtain affine coefficients by evaluation of each printed source expression.
    def expressions(q):
        x, y = q[:-1], q[-1]
        return ([y-alpha*sum(x[:-1])+(n-2)*alpha,
                 y-sum(x)+(n-1)]
                + [-y+alpha*x[i]+x[-1]-alpha for i in range(n-1)]
                + [-y+x[i] for i in range(n-1)] + [y]
                + [1-a for a in x] + [x[-1]-alpha])
    base = expressions([F(0)]*(n+1))
    vals = [expressions([F(j==i) for j in range(n+1)]) for i in range(n+1)]
    return {canonical(tuple(vals[i][k]-v for i in range(n+1))+(v,))
            for k,v in enumerate(base)}


def numerical_diagnostic(vertices, exact_facets, expected):
    # Floating point is confined to this independently labelled diagnostic.
    import numpy as np
    from scipy.spatial import ConvexHull
    points = np.array([[float(x) for x in v] for v in vertices])
    hull = ConvexHull(points)
    recovered = set()
    for simplex in hull.simplices:
        row = kernel_one([vertices[i] + (F(1),) for i in simplex])
        if row is None:
            # Qhull may emit zero-volume simplices on nonsimplicial facets.
            continue
        if max(evaluate(row,v) for v in vertices) <= 0:
            row = tuple(-x for x in row)
        recovered.add(row)
    assert recovered == set(exact_facets)
    assert abs(hull.volume-float(expected)) <= 1e-8*max(1,abs(float(expected)))
    return len(hull.simplices)


def check_case(label, n, lower, upper, c1=False):
    alpha, beta = lower[-2]/upper[-2], lower[-1]/upper[-1]
    original = [tuple(x)+(prod(x),) for x in product(*zip(lower,upper))]
    normalized = [tuple(x[i]/upper[i] for i in range(n))+(x[-1]/prod(upper),)
                  for x in original]
    facets = all_facets(normalized)
    rows = normalized_system(n,alpha,beta)
    assert all(evaluate(row,v)>=0 for row in rows for v in normalized)
    assert set(facets) <= rows, ('missing facet',label,set(facets)-rows)
    volume, simplices = exact_volume(normalized,facets)
    volume *= prod(upper)**2
    expected = proposed_volume(n,lower,upper)
    assert volume == expected, (label,volume,expected)
    if c1:
        assert alpha == 0
        assert volume == source_c1_volume(n,lower,upper)
        paper_rows = source_c1_rows(n,beta)
        assert all(evaluate(row,v)>=0 for row in paper_rows for v in normalized)
        assert set(facets) <= paper_rows
    scales = list(upper)+[prod(upper)]
    original_facets = {canonical(tuple(F(row[i])/scales[i] for i in range(n+1))
                                 +(F(row[-1]),)): face
                       for row,face in facets.items()}
    scipy_simplices = numerical_diagnostic(original,original_facets,expected)
    print(f'{label}: n={n} alpha={alpha} beta={beta} exact_facets={len(facets)} '
          f'pulling_simplices={simplices} exact_volume={volume} '
          f'scipy_boundary_simplices={scipy_simplices} OK',flush=True)
    return volume


def main():
    count = 0
    for n in (2,3,4):
        for a,b in product((F(0),F(1,3),F(2,3)),repeat=2):
            lower = [F(0)]*(n-2)+[a,b]
            upper = [F(1)]*n
            check_case('grid',n,lower,upper,c1=(a==0))
            count += 1
    rng = Random(260501493)
    for n in (2,3,4):
        for trial in range(3):
            upper = [F(rng.randrange(2,12),rng.randrange(1,7)) for _ in range(n)]
            lower = [F(0)]*(n-2)+[upper[-2]*F(rng.randrange(1,9),10),
                                             upper[-1]*F(rng.randrange(1,9),10)]
            check_case(f'random-{trial}',n,lower,upper)
            count += 1
    example = check_case('search-seat-example',4,[F(0),F(0),F(1,2),F(1,2)],
                         [F(1)]*4)
    assert example == F(31,640)
    for n in (2,3,4):
        upper = [F(i+2,3) for i in range(n)]
        lower = [F(0)]*(n-1)+[upper[-1]*F(2,5)]
        check_case('paper-C1-fidelity',n,lower,upper,c1=True)
        count += 1
    # Symbolic consistency with the source's bilinear and C1 formulae.
    a,b = sp.symbols('a b')
    bracket = 1-a-b+a*b
    assert sp.expand(bracket-(1-a)*(1-b)) == 0
    al,be = sp.symbols('alpha beta')
    for n in range(2,21):
        A,B = 1-al,1-be
        integrated = (al*be*(sp.Rational(1,n-1)-sp.Rational(1,factorial(n-1)))
                      +(A*be+al*B)*(sp.Rational(1,n)-sp.Rational(1,factorial(n)))
                      +A*B*(sp.Rational(1,n+1)-sp.Rational(1,factorial(n+1))))
        closed = ((factorial(n)-1)+(factorial(n-1)-n)*(al+be)
                  +(2*factorial(n-2)-n*n+n+1)*al*be)/factorial(n+1)
        assert sp.expand(integrated-closed) == 0
    print('Symbolic envelope-integral/volume algebra: n=2..20 OK')
    print(f'Exact hull and volume cases: {count+1}; source C1 fidelity cases: 12')
    print('SciPy comparisons are numerical diagnostics; facet coverage and volumes above are exact.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

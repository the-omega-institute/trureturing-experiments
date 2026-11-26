"""Gauge-fixed Euclidean Regge Hessian on the permutation-symmetric momentum lines."""
from itertools import permutations
from pathlib import Path
from runpy import run_path
import sympy as sp

base = run_path(str(Path(__file__).with_name('ch39-regge-flat-cell.py')))


def gauge_minor(active):
    canonical = sum((v for _, _, v in base['canonical']()), sp.zeros(10))
    edges = base['EDGES']
    fixed = {1 << i for i in range(active)} | {1+(1 << i) for i in range(active,4)}
    kept = [mask for mask in range(1, 16) if mask not in fixed]
    assert len(kept) == 11
    poly = {i:sp.zeros(11) for i in range(-active,active+1)}
    for order in permutations(range(4)):
        vertices = [sum(1 << v for v in order[:k]) for k in range(5)]
        for (i,j),value in canonical.todok().items():
            a,b = [vertices[k] for k in edges[i]]
            c,d = [vertices[k] for k in edges[j]]
            if b-a in fixed or d-c in fixed:
                continue
            projection = sum(((c>>k)&1)-((a>>k)&1) for k in range(active))
            poly[projection][kept.index(b-a),kept.index(d-c)] += 240*value
    assert all(mat == poly[-p].T for p,mat in poly.items())
    return poly


def determinant(active):
    poly = gauge_minor(active)
    # The determinant is a symmetric Laurent polynomial, hence a polynomial in y=2-z-1/z.
    degree = 11*active
    points = []
    for a in range(1, degree+2):
        z = sp.Integer(a)
        num = sum((z**(p+active)*m for p,m in poly.items()), sp.zeros(11))
        points.append((2-z-1/z, num.det(method='domain-ge')/z**degree))
    # Newton interpolation in the falling-product basis.
    ys, ordinates = zip(*points)
    y = sp.Symbol('y')
    total, basis = sp.Poly(0,y), sp.Poly(1,y)
    column = list(ordinates)
    for k in range(len(points)):
        total += basis.mul_ground(column[0])
        basis *= sp.Poly(y-ys[k],y)
        column = [(b-a)/(ys[i+k+1]-ys[i]) for i,(a,b) in enumerate(zip(column,column[1:]))]
    return sp.factor(total.as_expr())


def verify():
    y = sp.Symbol('y')
    factors = (y+100, -(y-10)*(y-6)*(31*y**2+66*y-1800),
               (3*y**2-22*y+60)**2*(13*y**3-56*y**2-44*y+600),
               -(y-5)**2*(2*y**2-8*y+15)**3*(3*y**2-10*y-20))
    for r, factor in enumerate(factors, 1):
        result = determinant(r)
        c = sp.cancel(result/(y**6*factor))
        assert c.is_Integer and c > 0
        assert sp.Poly(factor,y).count_roots(0,4) == 0 and factor.subs(y,0)>0
        poly = gauge_minor(r)
        endpoint = sum(((-1)**abs(p)*mat for p, mat in poly.items()),sp.zeros(11))
        roots = endpoint.charpoly().as_poly().intervals(eps=sp.Rational(1,2))
        assert sum(n for (a,b),n in roots if b<0)==10 and sum(n for (a,b),n in roots if a>0)==1
        print(r, result, flush=True)


if __name__ == '__main__':
    verify()

#!/usr/bin/env python3
"""Exact finite shadow envelopes, with certified containment LPs; not a continuum proof."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
from hashlib import sha256
import sympy as sp


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def cross(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])


def hull(points):
    pts = sorted(set(points))
    lower, upper = [], []
    for chain, seq in ((lower, pts), (upper, reversed(pts))):
        for p in seq:
            while len(chain) >= 2 and cross(chain[-2], chain[-1], p) <= 0:
                chain.pop()
            chain.append(p)
    ans = tuple(lower[:-1]+upper[:-1])
    assert len(ans) >= 3
    assert all(cross(ans[i-1], ans[i], ans[(i+1) % len(ans)]) > 0 for i in range(len(ans)))
    return ans


def solve(a, b):
    """Exact square linear system; None for a singular matrix."""
    m = [list(map(Q, row))+[Q(rhs)] for row, rhs in zip(a, b)]
    n = len(m)
    for j in range(n):
        pivot = next((i for i in range(j, n) if m[i][j]), None)
        if pivot is None:
            return None
        m[j], m[pivot] = m[pivot], m[j]
        d = m[j][j]
        m[j] = [v/d for v in m[j]]
        for i in range(n):
            if i != j:
                d = m[i][j]
                m[i] = [v-d*w for v, w in zip(m[i], m[j])]
    return tuple(m[i][-1] for i in range(n))


lp_count = 0


def radius(a, b):
    """min lambda >= 0 with a subset t + lambda*b. Exact primal and dual witnesses.

    For each outward edge normal n of b, require n.t+lambda*h_b(n)>=h_a(n).
    A nonnegative dual combination has sum w*n=0, sum w*h_b=1 and gives
    lambda>=sum w*h_a, independently checking optimality of the primal point.
    """
    global lp_count
    a, b = hull(a), hull(b)
    rows, rhs = [], []
    for i, p in enumerate(b):
        q = b[(i+1) % len(b)]
        n = (q[1]-p[1], p[0]-q[0])
        hb = dot(n, p)
        assert all(dot(n, z) <= hb for z in b)
        rows.append((n[0], n[1], hb))
        rhs.append(max(dot(n, z) for z in a))
    rows.append((Q(0), Q(0), Q(1)))
    rhs.append(Q(0))
    feasible = []
    for ids in combinations(range(len(rows)), 3):
        v = solve([rows[i] for i in ids], [rhs[i] for i in ids])
        if v is not None and all(dot(row, v) >= h for row, h in zip(rows, rhs)):
            feasible.append(v)
    assert feasible
    v = min(feasible, key=lambda z: z[2])
    assert v[2] > 0
    active = [i for i in range(len(rows)-1) if dot(rows[i], v) == rhs[i]]
    dual_ok = False
    for ids in combinations(active, 3):
        w = solve([[rows[i][j] for i in ids] for j in range(3)], [0, 0, 1])
        if w is not None and all(x >= 0 for x in w):
            assert sum((w[k]*rhs[i] for k, i in enumerate(ids)), Q(0)) == v[2]
            dual_ok = True
            break
    assert dual_ok
    lp_count += 1
    return v[2]


def rotation(q):
    w, x, y, z = map(Q, q)
    n = w*w+x*x+y*y+z*z
    m = [[w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)],
         [2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
         [2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z]]
    return tuple(tuple(v/n for v in row) for row in m)


def finite_checks():
    frames = [rotation(q)[:2] for q in
              [(1, 0, 0, 0), (1, 1, 1, 0), (1, 2, 1, 1), (2, 1, 2, 2), (3, 1, 0, 2)]]
    frames.append(((Q(2, 3), Q(1, 3), Q(2, 3)), (Q(1, 3), Q(2, 3), Q(-2, 3))))
    frames.append(((Q(20, 29), Q(21, 29), Q(0)),
                   (Q(-21, 29), Q(20, 29), Q(0))))
    for frame in frames:
        assert [[dot(a, b) for b in frame] for a in frame] == [[1, 0], [0, 1]]
    octa = [tuple(Q(s if i == j else 0) for i in range(3)) for j in range(3) for s in (-1, 1)]
    cube = [tuple(map(Q, v)) for v in product((-1, 1), repeat=3)]
    tetra = [tuple(map(Q, v)) for v in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]]
    for name, vertices in [('octahedron', octa), ('cube', cube), ('tetrahedron', tetra)]:
        shadows = [hull([tuple(dot(row, v) for row in f) for v in vertices]) for f in frames]
        radii = [[radius(a, b) for b in shadows] for a in shadows]
        envelope = [min(row) for row in radii]
        m, big = min(envelope), max(envelope)
        nu_finite = max(1/r for row in radii for r in row)
        assert m > 0 and big == 1 and big/m == nu_finite
        if name == 'octahedron':
            assert nu_finite > 1
        # Proposition 8 for diameter, squared to avoid square roots.
        diam2 = [max(dot((p[0]-q[0], p[1]-q[1]), (p[0]-q[0], p[1]-q[1]))
                     for p in a for q in a) for a in shadows]
        assert nu_finite**2 <= max(diam2)/min(diam2)
        def f(a):
            return min(radius(a, b) for b in shadows)
        for a, fa in zip(shadows, envelope):
            t = (Q(7, 11), Q(-5, 13))
            assert f([tuple(p[i]+t[i] for i in range(2)) for p in a]) == fa
            assert f([tuple(Q(3, 2)*x for x in p) for p in a]) == Q(3, 2)*fa
            z = tuple(sum((p[i] for p in a), Q(0))/len(a) for i in range(2))
            contracted = [tuple(z[i]+Q(1, 2)*(p[i]-z[i]) for i in range(2)) for p in a]
            assert f(contracted) == fa/2 < fa
        print(f'{name}: frames={len(frames)}, all_pairs={len(frames)**2}, '
              f'min_F={m}, max_F={big}, nu_finite={nu_finite}, ratio_identity=OK')


def source_fidelity():
    """The source Theorem 1's displayed optimal pair, in Q(sqrt(2))."""
    r = sp.sqrt(2)
    m1 = sp.Matrix([[2, 1, 2], [1, 2, -2]])/3
    m2 = sp.Matrix([[1, 1, 0], [-1, 1, 0]])/r
    assert m1*m1.T == sp.eye(2) and sp.simplify(m2*m2.T) == sp.eye(2)
    c = 3*r/4
    assert c**2 == sp.Rational(9, 8)
    boundary = [sp.simplify(abs(c*v)) for v in m1]
    assert all(sp.simplify(r/2-x).is_nonnegative for x in boundary)
    assert any(x == r/2 for x in boundary)
    for mu in [sp.Rational(1), sp.Rational(21, 20), sp.Rational(53, 50)]:
        assert (c-mu).is_positive
        assert all(sp.simplify(r/2-abs(mu*v)).is_positive for v in m1)
    assert (sp.Rational(17, 16)-c).is_positive
    assert any(sp.simplify(abs(sp.Rational(17, 16)*v)-r/2).is_positive for v in m1)
    print('source Theorem 1 optimal pair: MM^T=I; strict-fit supremum=3*sqrt(2)/4; square=9/8')
    print('source value rational isolation: 53/50 < 3*sqrt(2)/4 < 17/16')
    c1 = (Q(152024884, 259375205), Q(0), Q(210152163, 259375205))
    c2 = tuple(Q(x, 10**10) for x in [6632738028, 6106948881, 3980949609])
    c3 = tuple(Q(x, 10**10) for x in [8193990033, 5298215096, 1230614493])
    assert dot(c1, c1) == 1
    assert all(Q(98, 100)**2 < dot(v, v) < Q(99, 100)**2 for v in (c2, c3))
    print('source Noperthedron generators: norm(C1)=1; 98/100 < norm(C2),norm(C3) < 99/100')


if __name__ == '__main__':
    source_fidelity()
    finite_checks()
    print(f'containment LPs with exact primal and dual certificates: {lp_count}')
    print('Scope: finite dictionaries only; no continuum or Noperthedron non-Rupert certificate')
    print('check.py sha256: '+sha256(Path(__file__).read_bytes()).hexdigest())
    print('ALL_CHECKS_OK')

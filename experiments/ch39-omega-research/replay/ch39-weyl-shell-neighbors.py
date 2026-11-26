"""Exact two-cell stiffness on edges of finite Behrend shells."""
from fractions import Fraction
from itertools import combinations, product
from math import isqrt, prod

import sympy as sp


def response10(weights):
    """Ten times the three numerators, with their common denominator suppressed."""
    q, u, d, l, e = weights
    return 6*(q+4*u+d+l+e), 90*q+10*l, 60*q+15*(u+d)


def multiply(*rows):
    return tuple(map(prod, zip(*rows, strict=True)))


def partners():
    f = Fraction
    x = (1, f(1, 2), 2, f(3, 4), f(21, 2))
    y = (1, f(19, 5), f(29, 25), f(111, 25), f(3, 5))
    for weights, common in [(x, f(39, 4)), (y, f(336, 25)),
                            (multiply(x, y), f(1233, 100))]:
        assert min(weights) > 0 and response10(weights) == (10*common,)*3
    r, s, v, z = sp.symbols('r s v z')
    l = 3*(r+s-2)/2
    h = 3*(v+z-2)/2
    a, b, c = response10(multiply((1, r, s, l, 12-3*r), (1, v, z, h, 12-3*v)))
    defect = (r-3)*(v-3)+2
    assert sp.expand(5*a-3*b-2*c-360*defect) == 0
    assert sp.expand((s+3*r-6)*(z+3*v-6)-12*(r+v)+68
                     -2*(b-c)/15-8*defect) == 0
    # A balanced pair followed by a cell with the first up ratio.
    # The hypercharge combination for three cells is positive on either branch.
    partner = 3-2/(r-3)
    triple_defect = (12-3*r)**2*(12-3*v)+3*r**2*v-12
    assert sp.cancel((triple_defect-12*v*(r-1)*(r-3)).subs(v, partner)) == 0
    print('Single- and double-cell common numerators:', f(39, 4), f(336, 25), f(1233, 100))


def sharp_cut():
    r, s, v, z = sp.symbols('r s v z')
    root = sp.sqrt(71)
    star = 6+2*root
    bound = (742-22*root)/445
    t = -sp.Rational(17, 178)+671*root/6319

    def row(r, s):
        return (1, r, s, 3*(r+s-2)/2, 12-3*r)

    a, b, c = (value/10 for value in response10(multiply(row(r, s), row(v, z))))
    h = sp.expand(t*a+(1-t)*b-bound*c)
    l = sp.Rational(9, 4)*(1-2*t/5)
    ss = 3*t/5+l-3*bound/2
    j = 39*t/5+l-3*bound/2-49*ss
    k = l-7*ss
    m = -108*t/5-2*l+7*ss*star
    sign1 = s+7*r-star
    sign2 = z+7*v-star
    correction = j*r*v+k*(r*z+v*s)+m*(r+v)
    assert sp.simplify(sp.expand(h-ss*sign1*sign2-correction)) == 0
    for positive in (t, 1-t, ss, k, m, j+3*m/4, j+m/2+k/4):
        assert positive.is_positive
    # With both up ratios in [2,4] and both down ratios in [0,1],
    # A - R B is multi-affine, so its 16 corners suffice.
    margins = [(a-bound*b).subs({r: x, v: y, s: d, z: e})
               for x, y, d, e in product((2, 4), (2, 4), (0, 1), (0, 1))]
    assert all(margin.is_positive for margin in margins)
    at_star = {r: 0, v: 0, s: star, z: star}
    assert sp.simplify((a-b).subs(at_star)) == sp.simplify((b-bound*c).subs(at_star)) == 0
    # The rational partners have opposite signs, giving a constructive sharp cut.
    assert (star-2-7*sp.Rational(1, 2)).is_positive
    assert (sp.Rational(29, 25)+7*sp.Rational(19, 5)-star).is_positive
    print('Sharp spread on an odd cycle:', bound, '=', bound.evalf(16))


def shell_cut(radius2, chord2, expected):
    pts = [p for p in product(range(isqrt(radius2)+1), repeat=3) if sum(x*x for x in p) == radius2]
    edges = [(i, j) for i, j in combinations(range(len(pts)), 2)
             if sum((x-y)**2 for x, y in zip(pts[i], pts[j])) <= chord2]
    adj = [0]*len(pts)
    for i, j in edges:
        adj[i] |= 1 << j
        adj[j] |= 1 << i
    # Fix the last vertex on one side. Gray order flips one vertex at a time.
    best = size = previous = witness = 0
    for i in range(1, 1 << (len(pts)-1)):
        mask = i ^ (i >> 1)
        flip = mask ^ previous
        neighbors = adj[flip.bit_length()-1]
        gain = neighbors.bit_count()-2*(neighbors & previous).bit_count()
        size += gain if flip & mask else -gain
        if size > best:
            best, witness = size, mask
        previous = mask
    assert sum(((witness >> i) ^ (witness >> j)) & 1 for i, j in edges) == best
    orbits = sorted({tuple(sorted(p)) for p in pts})
    orbit = [orbits.index(tuple(sorted(p))) for p in pts]
    invariant_best = max(sum(((mask >> orbit[i]) ^ (mask >> orbit[j])) & 1 for i, j in edges)
                         for mask in range(1 << (len(orbits)-1)))
    assert (len(pts), len(edges), best, invariant_best) == expected
    print('Shell, squared chord cutoff, edges, maximum equalities (free, permutation-invariant):',
          radius2, chord2, len(edges), best, invariant_best)


if __name__ == '__main__':
    partners()
    sharp_cut()
    shell_cut(17, 2, (9, 6, 5, 0))
    shell_cut(867, 804, (19, 135, 80, 78))

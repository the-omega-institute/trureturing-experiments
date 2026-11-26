"""Exact boundary curvature of two Weyl-weighted four-simplices sharing a tetrahedron."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations, product
from math import comb, prod
from pathlib import Path
from runpy import run_path

base = run_path(str(Path(__file__).with_name('ch39-weyl-four-simplex.py')))
CHARGES, color, parity = base['CHARGES'], base['color'], base['parity']
LINKS = tuple(e for e in combinations(range(1, 6), 2) if e != (4, 5))
SCALE = 24**(len(LINKS)-1)


def triangle(a, b, c):
    return ((a, b), (b, c), (c, a))


def cap(k):
    return (triangle(1, 2, k), triangle(1, 3, k), triangle(2, 3, k),
            ((1, k),), ((2, k),), ((3, k),))


FACES = (((1, 3),), ((2, 3),), triangle(1, 2, 3))+cap(4)+cap(5)
PERMS = {d: tuple((p, parity(p)) for p in permutations(range(d))) for d in (2, 3)}


def paired(a, b):
    if a[0] != b[0]:
        return ((a[1], b[1]), (a[2], b[2])), ()
    return (), ((a[1], b[1]), (a[2], b[2]))


def haar(dim, arms):
    """Integer coefficient, Kronecker pairs, alternating tensors, and a common denominator."""
    n = len(arms)
    if n == 0:
        return [(1, (), ())], 1
    if n == 2 and (dim == 2 or arms[0][0] != arms[1][0]):
        pairs, eps = paired(*arms)
        return [(1, pairs, eps)], dim
    if n == dim == 3 and len({sign for sign, _, _ in arms}) == 1:
        return [(1, (), (tuple(a for _, a, _ in arms), tuple(b for _, _, b in arms)))], 6
    if n == 4 and dim == 2:
        terms = []
        for j in (1, 2, 3):
            d1, e1 = paired(arms[0], arms[j])
            d2, e2 = paired(*(arms[k] for k in (1, 2, 3) if k != j))
            terms.append((1, d1+d2, e1+e2))
        return terms, 6
    if n == 4 and dim == 3:
        pos = [a for a in arms if a[0] > 0]
        neg = [a for a in arms if a[0] < 0]
        if len(pos) == len(neg) == 2:
            terms = []
            for s, t in product(((0, 1), (1, 0)), repeat=2):
                pairs = tuple((pos[i][end], neg[p[i]][end]) for end, p in ((1, s), (2, t)) for i in (0, 1))
                terms.append((3 if s == t else -1, pairs, ()))
            return terms, 24
    return [], 1


def index_sum(dim, size, boundary, terms):
    parent = list(range(size))

    def root(i):
        while i != parent[i]:
            i = parent[i]
        return i

    for _, pairs, _ in terms:
        for a, b in pairs:
            parent[root(a)] = root(b)
    for _, a, b in boundary:
        parent[root(a)] = root(b)
    eps = [tuple(root(i) for i in legs) for _, _, tensors in terms for legs in tensors]
    if any(len(set(legs)) != dim for legs in eps):
        return 0, 0
    constrained = {i for legs in eps for i in legs}
    free = sum(i == v and i not in constrained for i, v in enumerate(parent))
    charges = defaultdict(int)
    for sign, a, _ in boundary:
        charges[root(a)] += sign
    independent = sum(c*c for i, c in charges.items() if i not in constrained)
    fixed = [(i, c) for i, c in charges.items() if i in constrained]
    # Sum the traceless Cartan diag(1,-1[,0]) squared over all free index loops.
    variance = 2*independent*dim**(free-1) if independent else 0
    factor = dim**free
    labels = {}

    def visit(tensors):
        if not tensors:
            total = sum(c*((labels[i] == 0)-(labels[i] == 1)) for i, c in fixed)
            return factor, factor*total**2+variance
        legs = max(tensors, key=lambda e: sum(i in labels for i in e))
        remaining = list(tensors)
        remaining.remove(legs)
        new = [i for i in legs if i not in labels]
        mass = moment = 0
        for perm, sign in PERMS[dim]:
            if any(i in labels and labels[i] != v for i, v in zip(legs, perm)):
                continue
            labels.update(zip(legs, perm))
            z, h = visit(remaining)
            mass += sign*z
            moment += sign*h
            for i in new:
                del labels[i]
        return mass, moment

    return visit(eps)


@cache
def contract(dim, types, probe=(1, 2), faces=FACES):
    stars = defaultdict(list)
    size = 0
    for face, typ in zip(faces, types, strict=True):
        if typ:
            word = face if typ > 0 else tuple((b, a) for a, b in reversed(face))
            for j, (a, b) in enumerate(word):
                u, v = size+j, size+(j+1) % len(word)
                if a < b:
                    stars[a, b].append((1, u, v))
                else:
                    stars[b, a].append((-1, v, u))
            size += len(word)
    choices = []
    denominator = 1
    for edge in LINKS:
        if edge == probe:
            continue
        terms, scale = haar(dim, stars[edge])
        if not terms:
            return 0, 0
        choices.append(terms)
        denominator *= scale
    mass = moment = 0
    for terms in product(*choices):
        z, h = index_sum(dim, size, stars[probe], terms)
        sign = prod(k for k, _, _ in terms)
        mass += sign*z
        moment += sign*h
    assert SCALE % denominator == 0
    return mass*(SCALE//denominator), moment*(SCALE//denominator)


@cache
def allowed(labels):
    return sum(color(y) for y in labels) % 3 == 0 and sum(abs(y) in (1, 3) for y in labels) % 2 == 0


def caps():
    for a, b, c in product(CHARGES, repeat=3):
        u, v, w = a+b, c-a, -b-c
        if u in CHARGES and v in CHARGES and w in CHARGES:
            if allowed((u, -a, -b)) and allowed((v, a, -c)) and allowed((w, b, c)):
                yield a, b, c, u, v, w


def pack(labels):
    return (tuple(color(y) for y in labels), tuple(int(abs(y) in (1, 3)) for y in labels),
            sum(1 << i for i, y in enumerate(labels) if y))


def coefficients():
    parts = [(v, *pack(v)) for v in caps()]
    assert len(parts) == 563
    # Common units are SCALE**2 for mass and 60*SCALE**2 for each curvature.
    rows = [[0]*(1 << 15) for _ in range(4)]
    common = {vals: pack(vals) for vals in product(CHARGES, repeat=3)}
    seen = 0
    for j, (left, lc, lw, lm) in enumerate(parts):
        a, b, c = left[:3]
        for k, (right, rc, rw, rm) in enumerate(parts[:j+1]):
            d, e, f = right[:3]
            for g in CHARGES:
                u, v = g-b-e, -g-c-f
                core = common.get((u, v, g))
                if core is None or not allowed((u, -g, b, e)) or not allowed((v, g, c, f)):
                    continue
                seen += 1
                cc, cw, cm = core
                mz, mc = contract(3, cc+lc+rc)
                wz, mw = contract(2, cw+lw+rw)
                z = mz*wz
                values = (z, z*(g+a+d)**2, 15*mz*mw, 15*wz*mc)
                masks = [cm+(lm << 3)+(rm << 9)]
                if k < j:
                    masks.append(cm+(rm << 3)+(lm << 9))
                for idx in masks:
                    for row, value in zip(rows, values):
                        row[idx] += value
    assert seen == 889068
    return rows


def homogeneous_polynomials(rows, dimension, units):
    import sympy as sp
    t = sp.Symbol('t')
    polys = []
    for row, scale in zip(rows, units, strict=True):
        totals = [0]*(dimension+1)
        for i, value in enumerate(row):
            totals[i.bit_count()] += value
        polys.append(sp.Poly.from_list(list(reversed(totals)), t, domain='QQ')/sp.Integer(scale))
    return [sp.Poly(p, t, domain='QQ') for p in polys]


def positive_on_interval(poly):
    deg = poly.degree()
    coeff = [sum(poly.nth(i)*comb(deg-i, j-i)/12**i for i in range(j+1)) for j in range(deg+1)]
    return min(coeff) >= 0 and max(coeff) > 0


def verify(rows):
    import sympy as sp
    single = base['face_coefficients']()
    old_faces = [FACES.index(face) for face in base['FACES']]
    for subset in range(512):
        mask = sum(1 << j for i, j in enumerate(old_faces) if subset & (1 << i))
        for i, row in enumerate(rows):
            assert F(row[mask], SCALE**2*(60 if i else 1)) == single[i][subset]
    _, a, b, c = rows
    u = [x-y for x, y in zip(a, b)]
    v = [x-y for x, y in zip(b, c)]
    low, high = F(112, 5), F(1006712, 3025)
    for coeff in (v, [5*x-112*y for x, y in zip(u, v)],
                  [1006712*y-3025*x for x, y in zip(u, v)]):
        corners = [x*12**(15-i.bit_count()) for i, x in enumerate(coeff)]
        for j in range(15):
            bit = 1 << j
            for i in range(1 << 15):
                if i & bit:
                    corners[i] += corners[i ^ bit]
        assert min(corners) >= 0 and max(corners) > 0
    for active, ratio in [(0x7, low), (0x66f2, high)]:
        assert v[active] > 0 and u[active] == ratio*v[active]
        assert all(not (u[i] or v[i]) for i in range(len(v)) if i != active and i & active == i)
    t = sp.Symbol('t')
    homogeneous = homogeneous_polynomials(rows, 15, [SCALE**2]+[60*SCALE**2]*3)
    z, *h = homogeneous
    assert z.nth(3) == 14 and z.nth(4) == 18
    assert tuple(p.nth(3) for p in h) == (sp.Rational(71, 15), 1, sp.Rational(5, 6))
    kt = [4*z+sp.Poly(1+30*t,t)*p.exquo(sp.Poly(t,t)) for p in h]
    top, bottom = kt[0], kt[2]
    derivative = top.diff()*bottom-bottom.diff()*top
    assert positive_on_interval(derivative)
    lo, hi = sp.Rational('0.076393202250020'), sp.Rational('0.076393202250022')
    assert (3-10*hi)**2 < 5 < (3-10*lo)**2
    assert sp.Rational('1.01904359786172') < top.eval(lo)/bottom.eval(lo)
    assert top.eval(hi)/bottom.eval(hi) < sp.Rational('1.01904359786175')
    print('Complete shared-tetrahedron slope interval:', low, high)
    print('Equal strengths, upper endpoint of K1/K3 between 1.01904359786172 and 1.01904359786175.')


if __name__ == '__main__':
    verify(coefficients())

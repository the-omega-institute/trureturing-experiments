"""Exact Haar integral on the ten triangles of a four-simplex with one boundary loop."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations, permutations, product
from math import factorial

LINKS = tuple(combinations(range(1, 5), 2))
OUTER = tuple(combinations(range(1, 5), 3))
# After fixing the tree at vertex 0, the other five rooted triangles are monogons.
FACES = tuple((edge,) for edge in LINKS[1:]) + tuple(
    ((a, b), (b, c), (c, a)) for a, b, c in OUTER)
CHARGES = (0, 1, -1, 2, -2, 3, -3, 4, -4, 6, -6)


def incidence(face, edge):
    return int(edge in face)-int(edge[::-1] in face)


def color(y):
    if abs(y) in (1, 4):
        return 1 if y > 0 else -1
    if abs(y) == 2:
        return -1 if y > 0 else 1
    return 0


def parity(p):
    return (-1)**sum(a > b for a, b in combinations(p, 2))


@cache
def contraction(dim, types):
    """Return integral at identity and minus its boundary Cartan second derivative."""
    stars = defaultdict(list)
    size = 0
    for face, typ in zip(FACES, types, strict=True):
        if not typ:
            continue
        word = face if typ > 0 else tuple((b, a) for a, b in reversed(face))
        for j, (a, b) in enumerate(word):
            u, v = size+j, size+(j+1) % len(word)
            if a < b:
                stars[a, b].append((1, u, v))
            else:
                stars[b, a].append((-1, v, u))
        size += len(word)

    fixed = []
    scales = 1
    alternatives = []
    for link in LINKS[1:]:
        arms = stars[link]
        if not arms:
            continue
        if len(arms) == 2 and arms[0][0] != arms[1][0]:
            _, a, b = arms[0]
            _, c, d = arms[1]
            fixed.extend(((a, c), (b, d)))
            scales *= dim
        elif len(arms) == dim and len({s for s, _, _ in arms}) == 1:
            # Two equal SU(2) or three equal SU(3) factors. Expand the product of determinants.
            alternatives.append([(parity(p), tuple((arms[i][1], arms[j][2]) for i, j in enumerate(p)))
                                 for p in permutations(range(dim))])
            scales *= factorial(dim)
        else:
            return F(0), F(0)

    boundary = stars[LINKS[0]]
    fixed.extend((a, b) for _, a, b in boundary)
    mass = moment = 0
    for terms in product(*alternatives):
        parent = list(range(size))

        def root(i):
            while parent[i] != i:
                i = parent[i]
            return i

        def glue(pairs):
            for a, b in pairs:
                parent[root(a)] = root(b)

        glue(fixed)
        sign = 1
        for factor, pairs in terms:
            sign *= factor
            glue(pairs)
        loops = sum(i == v for i, v in enumerate(parent))
        mass += sign*dim**loops
        charges = defaultdict(int)
        for orientation, a, _ in boundary:
            charges[root(a)] += orientation
        if charges:
            moment += sign*2*dim**(loops-1)*sum(v*v for v in charges.values())
    return F(mass, scales), F(moment, 4*scales)


def face_coefficients():
    rows = [[F(0) for _ in range(512)] for _ in range(4)]
    outer_incidence = [tuple(incidence(face, e) for face in FACES[5:]) for e in LINKS]
    for outer in product(CHARGES, repeat=4):
        flux = [sum(s*y for s, y in zip(signs, outer, strict=True)) for signs in outer_incidence]
        local = tuple(-q for q in flux[1:])
        if any(y not in CHARGES for y in local):
            continue
        labels = local+outer
        c, hc = contraction(3, tuple(color(y) for y in labels))
        if not c and not hc:
            continue
        w, hw = contraction(2, tuple(int(abs(y) in (1, 3)) for y in labels))
        face_set = sum(1 << i for i, label in enumerate(labels) if label)
        for line, value in zip(rows, (c*w, c*w*F(flux[0]**2, 60), c*hw, w*hc), strict=True):
            line[face_set] += value
    return rows


def stiffness_ratios(rows):
    import sympy as sp
    from math import comb
    t = sp.Symbol('t')
    z, *hessians = [sp.Poly.from_list(list(reversed(row)), t, domain='QQ') for row in rows]
    strengths = [4*z + sp.Poly(1+30*t, t)*(h.exquo(sp.Poly(t, t))) for h in hessians]
    # Algebraic end of the positive local weight interval.
    low = sp.Rational('0.076393202250020')
    high = sp.Rational('0.076393202250022')
    assert (3-10*high)**2 < 5 < (3-10*low)**2 and high < sp.Rational(1, 12)
    for larger, smaller, limits in [
        (0, 2, ('1.01264026101293', '1.01264026101295')),
        (0, 1, ('1.01209551929217', '1.01209551929220')),
        (1, 2, ('1.00053823153088', '1.00053823153091')),
    ]:
        top, bottom = strengths[larger], strengths[smaller]
        derivative = top.diff()*bottom - top*bottom.diff()
        deg = derivative.degree()
        # Coefficients of (1+x)^deg N(x/[12(1+x)]) on x>0.
        transformed = [sum(derivative.nth(i)*comb(deg-i, j-i)/12**i for i in range(j+1))
                       for j in range(deg+1)]
        assert all(c >= 0 for c in transformed) and any(c > 0 for c in transformed)
        assert top.eval(0) == bottom.eval(0) == 4
        assert sp.Rational(limits[0]) < top.eval(low)/bottom.eval(low)
        assert top.eval(high)/bottom.eval(high) < sp.Rational(limits[1])
        print('Channels', larger+1, smaller+1, 'sharp upper endpoint in', *limits)
    # The two tetrahedral caps and the three closed tetrahedra.
    dims = (6, 3, 3, 2, 1)
    q, u, d, l, e = (F(1, dim)**3 for dim in dims)
    one_cap = (3*(q+4*u+d+l+e)/5, 9*q+l, 6*q+3*(u+d)/2)
    assert tuple(h.nth(3) for h in hessians) == tuple(4*v for v in one_cap)
    assert z.nth(4) == 6*sum(F(1, dim**2) for dim in dims) == 9


def heterogeneous_threshold(rows):
    _, hyper, weak, color_row = rows
    u = [h-w for h, w in zip(hyper, weak, strict=True)]
    v = [w-c for w, c in zip(weak, color_row, strict=True)]
    lower, upper = F(112, 5), F(776, 25)
    for coeff in (v, [a-lower*b for a, b in zip(u, v)], [upper*b-a for a, b in zip(u, v)]):
        scores = [value*F(1, 12)**i.bit_count() for i, value in enumerate(coeff)]
        for j in range(9):
            bit = 1 << j
            for i in range(512):
                if i & bit:
                    scores[i] += scores[i ^ bit]
        assert min(scores) >= 0 and max(scores) > 0
    # One three-face cap, or a five-face disk, can be left active to attain either slope.
    for active, target in [(0x25, lower), (0xb6, upper)]:
        assert v[active] > 0 and u[active] == target*v[active]
        assert all(not (u[i] or v[i]) for i in range(512) if i != active and i & active == i)
    sw, em = F('0.23129'), F('127.930')
    inputs = [F(3, 5)*(1-sw)*em, sw*em, 1/F('0.1180')]
    slopes = [F(41, 10), F(-19, 6), F(-7)]
    direction = (slopes[0]-slopes[1])/(slopes[1]-slopes[2])
    mismatch = inputs[0]-inputs[1]-direction*(inputs[1]-inputs[2])
    assert direction == F(218, 115) < lower and mismatch < 0
    crossing = (inputs[0]-inputs[2])/(slopes[0]-slopes[2])
    target = [a-b*crossing for a, b in zip(inputs, slopes)]
    assert target[1]/target[0] == F(4322773606174, 3962837583051)
    assert slopes[2]*inputs[1]-slopes[1]*inputs[2] < 0 < slopes[0]*inputs[1]-slopes[1]*inputs[0]
    print('Heterogeneous sharp slope interval:', lower, upper, 'one-loop separation:', mismatch)


if __name__ == '__main__':
    per_face = face_coefficients()
    rows = [[sum(v for mask, v in enumerate(line) if mask.bit_count() == degree) for degree in range(10)]
            for line in per_face]
    for label, coeffs in zip(('mass', 'hypercharge', 'weak', 'color'), rows, strict=True):
        print(label, ', '.join(map(str, coeffs)))
    stiffness_ratios(rows)
    heterogeneous_threshold(per_face)

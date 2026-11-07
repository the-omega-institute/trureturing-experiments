#!/usr/bin/env python3
"""Exact shadow checks for the section-3 scalar conjecture in arXiv:2609.10788v1.

Only integers, rationals and exact real algebraic expressions are used. The
square minimization below considers all orientations, not a mesh of angles.
"""
from itertools import combinations
from fractions import Fraction
from math import isqrt
import sympy as s

Q = s.Rational
sqrt = s.sqrt


def simp(x):
    return s.simplify(x)


def eq(a, b):
    assert simp(a - b) == 0, (a, b)


def algebraic_interval(x, bits):
    """Rational outward interval, including nested square roots and inverses."""
    if x.is_Rational:
        value = Fraction(int(x.p), int(x.q))
        return value, value
    if x.is_Add:
        bounds = [algebraic_interval(a, bits) for a in x.args]
        return sum(a for a, _ in bounds), sum(b for _, b in bounds)
    if x.is_Mul:
        lo = hi = Fraction(1)
        for arg in x.args:
            a, b = algebraic_interval(arg, bits)
            products = [lo*a, lo*b, hi*a, hi*b]
            lo, hi = min(products), max(products)
        return lo, hi
    if x.is_Pow:
        base, power = x.args
        lo, hi = algebraic_interval(base, bits)
        if power == Q(1, 2):
            if lo < 0:
                raise ArithmeticError("refine square-root domain")
            scale = 1 << bits
            a = isqrt((lo.numerator * scale**2) // lo.denominator)
            b = isqrt((hi.numerator * scale**2) // hi.denominator)
            return Fraction(a, scale), Fraction(b+1, scale)
        if power.is_Integer:
            exponent = int(power)
            if exponent < 0:
                if lo <= 0 <= hi:
                    raise ArithmeticError("refine reciprocal domain")
                lo, hi = 1/hi, 1/lo
                exponent = -exponent
            candidates = [lo**exponent, hi**exponent]
            if exponent % 2 == 0 and lo <= 0 <= hi:
                candidates.append(Fraction(0))
            return min(candidates), max(candidates)
    raise AssertionError(("unsupported exact interval expression", x))


def nonnegative(x):
    x = simp(s.sqrtdenest(s.radsimp(x)))
    if x == 0:
        return True
    for bits in [16, 32, 64, 128, 256, 512, 1024]:
        try:
            lo, hi = algebraic_interval(x, bits)
        except ArithmeticError:
            continue
        if lo >= 0:
            return True
        if hi < 0:
            return False
    raise AssertionError(("undecided exact sign", x))


def maximum(values):
    best = values[0]
    for value in values[1:]:
        if nonnegative(value - best):
            best = value
    return simp(best)


def minimum(values):
    return -maximum([-value for value in values])


def shadow_vertices(matrix):
    return [sign * matrix[:, i] for i in range(3) for sign in [-1, 1]]


def diameter_squared(vertices):
    # The convex hull has the same diameter as its finite vertex set:
    # distances between convex combinations are bounded by the largest
    # pairwise vertex distance, by the triangle inequality.
    return maximum([simp((a - b).dot(a - b))
                    for a in vertices for b in vertices])


def smallest_square(vertices):
    """Exact minimum side for a centrally symmetric finite planar set.

    For a fixed orthonormal frame (u, Ju), central symmetry makes the two
    coordinate ranges symmetric. Translation cannot decrease either range,
    so the best square has center zero and side 2 max_v(|u.v|, |Ju.v|).
    This is twice the support function of K = {v, -v, Jv, -Jv}.

    On every arc where a single k in K supports, h(theta)=k.u(theta)>0
    and h''=-h<0. Its minimum is therefore at an endpoint of the arc.
    Endpoints satisfy (k-l).u=0. Enumerating these exact algebraic normals
    (and quadrant endpoints) includes every minimum. Fourfold symmetry
    restricts the enumeration to the closed first quadrant.
    """
    J = s.Matrix([[0, -1], [1, 0]])
    supports = set()
    for v in vertices:
        for w in [v, -v, J * v, -J * v]:
            supports.add(tuple(map(simp, w)))
    supports = [s.Matrix(v) for v in sorted(supports, key=str)]
    normals = {(s.Integer(1), s.Integer(0)),
               (s.Integer(0), s.Integer(1))}
    for a, b in combinations(supports, 2):
        d = a - b
        n = s.Matrix([-d[1], d[0]])
        norm2 = simp(n.dot(n))
        if norm2 == 0:
            continue
        for sign in [-1, 1]:
            w = sign * n
            if nonnegative(w[0]) and nonnegative(w[1]):
                normals.add(tuple(simp(t / sqrt(norm2)) for t in w))
    sides = []
    for normal in normals:
        u = s.Matrix(normal)
        eq(u.dot(u), 1)
        sides.append(2 * maximum([simp(v.dot(u)) for v in supports]))
    return simp(minimum(sides)), len(normals)


H = s.Matrix([[1 / sqrt(2), -1 / sqrt(2), 0],
              [1 / sqrt(6), 1 / sqrt(6), -2 / sqrt(6)]])
I = Q(1, 3) * s.Matrix([[-1, 2, 2], [2, -1, 2]])
S = s.Matrix([[1, 0, 0], [0, 1, 0]])
R = s.Matrix([[0, 1 / sqrt(2), -1 / sqrt(2)],
              [-1 / sqrt(2), Q(1, 2), Q(1, 2)]])
# These are literally the two matrices in the source's proof of Theorem 1.
M1 = Q(1, 3) * s.Matrix([[2, 1, 2], [1, 2, -2]])
M2 = s.Matrix([[1, 1, 0], [-1, 1, 0]]) / sqrt(2)

expected = {
    "H": (H, 2 * sqrt(Q(2, 3)), 1 + 1 / sqrt(3)),
    "I": (I, 4 * sqrt(2) / 3, Q(4, 3)),
    "S": (S, s.Integer(2), sqrt(2)),
    "Q": (R, sqrt(3), sqrt(2)),
    "source_M2": (M2, s.Integer(2), sqrt(2)),
}
values = {}
for name, (matrix, diameter, side) in expected.items():
    assert (matrix * matrix.T - s.eye(2)).applyfunc(simp) == s.zeros(2)
    vertices = shadow_vertices(matrix)
    eq(diameter_squared(vertices), diameter**2)
    computed_side, count = smallest_square(vertices)
    eq(computed_side, side)
    values[name] = (diameter, side)
    print(f"{name}: MM^T=I_2; diameter={diameter}; square_side={side}; "
          f"exact_orientation_candidates={count}")

assert M1 * M1.T == s.eye(2)
assert set(map(tuple, shadow_vertices(M1))) == {
    (sign * a, sign * b) for sign in [-1, 1]
    for a, b in [(Q(2, 3), Q(1, 3)), (Q(1, 3), Q(2, 3)),
                 (Q(2, 3), Q(-2, 3))]}
print("FIDELITY: source M1 vertices and source M2 square side sqrt(2) reproduced")

c = Q(23, 10)
q = 3 * sqrt(2) / 4
phi = lambda name, weight: values[name][0] + weight * values[name][1]
delta = simp(phi("H", c) - q * phi("I", c))
eq(10 * sqrt(3) * delta, 23 + 20 * sqrt(2) + 3 * sqrt(3) - 23 * sqrt(6))
# Positive radicals, compared by exact rational squares.
assert Q(707, 500)**2 < 2
assert Q(433, 250)**2 < 3
assert Q(49, 20)**2 > 6
assert Q(7, 4)**2 > 3
numerator_lower = 23 + 20 * Q(707, 500) + 3 * Q(433, 250) - 23 * Q(49, 20)
eq(numerator_lower, Q(63, 500))
eq(numerator_lower / (10 * Q(7, 4)), Q(9, 1250))
assert nonnegative(delta - Q(9, 1250)) and simp(delta - Q(9, 1250)) != 0
print(f"COUNTEREXAMPLE: Delta={delta}; exact certificate Delta > 9/1250 > 0")

# Additional algebra for Triage; the settling target remains the fixed pair.
cminus = simp((4 * sqrt(2) / 3 - sqrt(3)) / (sqrt(2) - Q(4, 3)))
cplus = simp((2 - 2 * sqrt(Q(2, 3))) / (1 + 1 / sqrt(3) - sqrt(2)))
t = s.Symbol("t", real=True)
eq(phi("S", t), q * phi("I", t))
eq(phi("I", t) - phi("Q", t), (sqrt(2) - Q(4, 3)) * (cminus - t))
eq(phi("H", t) - phi("S", t), (1 + 1 / sqrt(3) - sqrt(2)) * (t - cplus))
assert nonnegative(sqrt(2) - Q(4, 3)) and sqrt(2) != Q(4, 3)
assert nonnegative(1 + 1 / sqrt(3) - sqrt(2))
assert nonnegative(c - cplus) and simp(c - cplus) != 0
print(f"TRIAGE: c_minus={cminus}; c_plus={cplus}; c=23/10 > c_plus")
# Conditional on the geometric 2-Lipschitz bounds stated in README, the
# following exact budget propagates the gap to nearby orthogonal matrices.
loss = (2 + 2 * c) / 1000
eq(loss, Q(33, 5000))
eq(Q(9, 1250) - loss, Q(3, 5000))
print("TRIAGE: 1/1000 perturbation budget leaves gap > 3/5000")
print("ALL_CHECKS_OK")

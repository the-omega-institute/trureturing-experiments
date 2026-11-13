#!/usr/bin/env python3
"""Exact checks of CLN (2021), Theorem 3.1, and the proposed rank-3 refutation.

No sampling result is a proof of the universal conjecture or its negation.
Q(sqrt(3)) comparisons are decided by rational signs and squaring.
SymPy is used only for exact algebraic projection identities.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import total_ordering
import sympy as S


@total_ordering
@dataclass(frozen=True)
class Q3:
    a: F = F(0)
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    @staticmethod
    def coerce(x):
        return x if isinstance(x, Q3) else Q3(F(x))

    def __add__(self, other):
        o = self.coerce(other)
        return Q3(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __neg__(self):
        return Q3(-self.a, -self.b)

    def __sub__(self, other):
        return self + -self.coerce(other)

    def __rsub__(self, other):
        return self.coerce(other) + -self

    def __mul__(self, other):
        o = self.coerce(other)
        return Q3(self.a*o.a + 3*self.b*o.b, self.a*o.b + self.b*o.a)
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = self.coerce(other)
        den = o.a*o.a - 3*o.b*o.b
        assert den != 0
        return self * Q3(o.a/den, -o.b/den)

    def sign(self):
        a, b = self.a, self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if a > 0 and b > 0:
            return 1
        if a < 0 and b < 0:
            return -1
        d = a*a - 3*b*b
        assert d != 0  # sqrt(3) is irrational; a,b here are rational.
        return ((d > 0) - (d < 0)) * ((a > 0) - (a < 0))

    def __eq__(self, other):
        o = self.coerce(other)
        return self.a == o.a and self.b == o.b

    def __lt__(self, other):
        return (self - other).sign() < 0


def add(x, y):
    return (x[0] + y[0], x[1] + y[1])


def neg(x):
    return (-x[0], -x[1])


def sub(x, y):
    return add(x, neg(y))


def smul(a, x):
    return (a*x[0], a*x[1])


def dot(x, y):
    return x[0]*y[0] + x[1]*y[1]


def cross(x, y):
    return x[0]*y[1] - x[1]*y[0]


R = Q3(0, 1)
C = [(Q3(1), Q3()), (Q3(F(-1, 2)), Q3(0, F(1, 2))),
     (Q3(F(-1, 2)), Q3(0, F(-1, 2)))]
O = (Q3(), Q3())
assert R*R == 3
for i in range(3):
    assert dot(C[i], C[i]) == 1
    for j in range(i):
        assert dot(sub(C[i], C[j]), sub(C[i], C[j])) == 3

# Stereographic rational parametrization of the unit circle; include the missing point.
directions = {(-F(1), F(0))}
for p in range(-12, 13):
    for q in range(1, 13):
        t = F(p, q)
        directions.add(((1-t*t)/(1+t*t), 2*t/(1+t*t)))
U = [(Q3(x), Q3(y)) for x, y in sorted(directions)]
assert all(dot(u, u) == 1 for u in U)

# Source fidelity: verify the displayed 2x2 Hermitian block characteristic polynomial
# symbolically, then Theorem 3.1's lambda_2 - lambda_5 = 2r at exact unit directions.
x, y, r, lam, alpha = S.symbols('x y r lam alpha', real=True)
H = S.Matrix([[alpha, r*(x-S.I*y)], [r*(x+S.I*y), alpha]])
assert S.expand((lam*S.eye(2)-H).det() - ((lam-alpha)**2-r*r*(x*x+y*y))) == 0
fidelity_count = 0
for radius in [R/2, Q3(1), R, Q3(2)]:
    for u in U:
        alphas = [dot(u, c) for c in C]
        vals = sorted([a+radius for a in alphas] + [a-radius for a in alphas], reverse=True)
        assert vals[1] - vals[4] == 2*radius
        assert min(a+radius for a in alphas) >= max(a-radius for a in alphas)
        fidelity_count += 1
print(f'FIDELITY: symbolic block characteristic polynomial; Theorem 3.1 in {fidelity_count} exact cases')


def inside(z):
    return all(dot(sub(z, c), sub(z, c)) <= 3 for c in C)


def support_certificate(u):
    """Exhibit an attaining point and certify a global support bound on all of K.

    Arc certificate: K is contained in its defining disk, so its support <= c.u+R.
    Vertex certificate: u is a nonnegative combination of the two active disk normals.
    Thus these checks do not infer a global maximum from a point cloud.
    """
    for c in C:
        z = add(c, smul(R, u))
        if inside(z):
            h = dot(u, z)
            assert h == dot(u, c) + R
            return h, z, 'arc'
    for i, v in enumerate(C):
        c1, c2 = [c for j, c in enumerate(C) if j != i]
        n1, n2 = sub(v, c1), sub(v, c2)
        det = cross(n1, n2)
        a, b = cross(u, n2)/det, cross(n1, u)/det
        if a >= 0 and b >= 0:
            assert add(smul(a, n1), smul(b, n2)) == u
            assert inside(v)
            assert dot(n1, n1) == 3 and dot(n2, n2) == 3
            return dot(u, v), v, 'vertex'
    raise AssertionError('No valid global support certificate')

kinds = {'arc': 0, 'vertex': 0}
for u in U:
    hi, zhi, kind = support_certificate(u)
    minus_lo, zlo, _ = support_certificate(neg(u))
    assert hi + minus_lo == R
    assert inside(zhi) and inside(zlo)
    kinds[kind] += 1
# Positive interior: B(0,1/2) lies in K, since |c|=1 and 1+1/2 < sqrt(3).
assert Q3(F(3, 2)) < R
# Any disk of width R would have radius R/2. Two vertices distance R apart force
# its centre to their midpoint, which is too far from the third vertex.
mid = smul(Q3(F(1, 2)), add(C[0], C[1]))
assert dot(sub(C[2], mid), sub(C[2], mid)) == Q3(F(9, 4))
assert Q3(F(9, 4)) > (R/2)*(R/2)
print(f'GEOMETRY: {len(U)} rational unit directions, width sqrt(3), global support certificates {kinds}; interior and non-disk obstruction exact')

# Family: k repeated disk blocks, then centroid scalar padding. All lower branches
# lie below all upper branches; the padding lies between them. Check ties as well.
family_count = 0
for k in range(3, 9):
    centers = C + [C[0]]*(k-3)
    for padding in range(4):
        n = 2*k + padding
        for u in U:
            upper = [dot(u, c)+R for c in centers]
            lower = [dot(u, c)-R for c in centers]
            assert max(lower) <= 0 <= min(upper)
            spectrum = sorted(upper+lower+[Q3()]*padding, reverse=True)
            assert len(spectrum) == n and spectrum[k-1] == min(upper)
            family_count += 1
print(f'FAMILY: kth spectral threshold for k=3..8, n=2k..2k+3, {family_count} exact cases')

# Construct genuine rank-k orthogonal projection compressions, independent of the
# support calculation. Do it at the centroid and the three vertices, with radicals.
rt = S.sqrt(3)
w = (-1+S.I*rt)/2
centers = [S.Integer(1), w, S.conjugate(w)]

def clean(a):
    return S.simplify(S.expand_complex(a))


def zero_matrix(M):
    assert all(clean(e) == 0 for e in M)


def projection(z, cs):
    blocks = []
    for c in cs:
        d = clean(z-c)
        norm2 = clean(S.conjugate(d)*d)
        a = S.sqrt((1+S.sqrt(1-norm2/3))/2)
        b = d/(2*rt*a)
        v = S.Matrix([a, b])
        assert clean((v.conjugate().T*v)[0]-1) == 0
        B = S.Matrix([[c, 2*rt], [0, c]])
        assert clean((v.conjugate().T*B*v)[0]-z) == 0
        blocks.append(v*v.conjugate().T)
    return S.diag(*blocks)

projection_count = 0
for z in [S.Integer(0), *centers]:
    P = projection(z, centers)
    A = S.diag(*[S.Matrix([[c, 2*rt], [0, c]]) for c in centers])
    zero_matrix(P*P-P)
    zero_matrix(P.conjugate().T-P)
    zero_matrix(P*A*P-z*P)
    assert clean(S.trace(P)-3) == 0  # idempotent over C: rank = trace
    # Scalar padding and additional disk block: rank 4 in dimension 9.
    P4 = S.diag(P, projection(z, [centers[0]]), S.Integer(0))
    A4 = S.diag(A, S.Matrix([[centers[0], 2*rt], [0, centers[0]]]), S.Integer(0))
    zero_matrix(P4*P4-P4)
    zero_matrix(P4.conjugate().T-P4)
    zero_matrix(P4*A4*P4-z*P4)
    assert clean(S.trace(P4)-4) == 0
    # Nonzero complex scaling and scalar translation preserve these compressions.
    zero_matrix(P*((2+S.I)*A+(3-S.I)*S.eye(6))*P-((2+S.I)*z+3-S.I)*P)
    projection_count += 2
print(f'COMPRESSION: {projection_count} explicit algebraic projections (rank 3/6 and rank 4/9); affine transport exact')
print('SCOPE: finite exact checks and symbolic identities; universal family proof remains a preregistered outline, with no Lean proof')
print('ALL_CHECKS_OK')

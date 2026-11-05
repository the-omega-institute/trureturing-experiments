#!/usr/bin/env python3
"""Exact checks for Bouafia's negative-type diversity-deletion question.

Source: P. Bouafia, Magnitude and diversity of trees, arXiv:2605.03681v1,
Proposition 5 (weighted-tree weights), Theorem 4, and the following open problem.
Run: python3 check.py. Requires SymPy; no floating point is used.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import sympy as s


def solve(a, b):
    """Independent Fraction Gauss-Jordan solver; does not call SymPy."""
    n = len(a)
    rows = [[F(x) for x in row] + [F(y)] for row, y in zip(a, b)]
    for j in range(n):
        k = next(k for k in range(j, n) if rows[k][j])
        rows[j], rows[k] = rows[k], rows[j]
        c = rows[j][j]
        rows[j] = [x / c for x in rows[j]]
        for k in range(n):
            if k != j:
                c = rows[k][j]
                rows[k] = [x - c*y for x, y in zip(rows[k], rows[j])]
    return [row[-1] for row in rows]


def mul(a, v):
    return [sum((x*y for x, y in zip(row, v)), F(0)) for row in a]


def energy(a, v):
    return sum((x*y for x, y in zip(v, mul(a, v))), F(0))


def deletion(a):
    active = list(range(len(a)))
    history = []
    for _ in range(len(a)):
        assert active
        w = solve([[a[i][j] for j in active] for i in active], [1]*len(active))
        full = [F(0)]*len(a)
        for i, value in zip(active, w):
            full[i] = value
        history.append((tuple(active), full))
        if all(x >= 0 for x in w):
            total = sum(w)
            assert total > 0
            return history, [x/total for x in full]
        active = [i for i, value in zip(active, w) if value > 0]
    raise AssertionError('Nontermination')


def zero(expr):
    assert s.cancel(expr) == 0, expr


def bernstein_positive(expr, q):
    """Exact interval certificate: all Bernstein coefficients on [19/20,1] > 0."""
    t = s.Symbol('t')
    a = s.Rational(19, 20)
    poly = s.Poly(s.expand(expr.subs(q, a+(1-a)*t)), t)
    degree = poly.degree()
    coeffs = [sum(poly.nth(j)*s.Rational(comb(k, j), comb(degree, j))
                  for j in range(k+1)) for k in range(degree+1)]
    assert all(c > 0 for c in coeffs), coeffs
    reconstruction = sum(coeffs[k]*comb(degree, k)*t**k*(1-t)**(degree-k)
                         for k in range(degree+1))
    zero(reconstruction-poly.as_expr())
    return min(coeffs)


q0 = F(19, 20)
# Fidelity: Proposition 5's weighted-tree formula, for a four-leaf equal-edge star.
# Each edge has length -log(19/20), so tree similarities are exact powers of q0.
star_dist = [[0 if i == j else 1 if i == 0 or j == 0 else 2
              for j in range(5)] for i in range(5)]
star = [[q0**d for d in row] for row in star_dist]
star_formula = [4/(1+q0)-3] + [1/(1+q0)]*4
assert solve(star, [1]*5) == star_formula
assert mul(star, star_formula) == [1]*5
star_history, star_p = deletion(star)
assert [a for a, _ in star_history] == [(0, 1, 2, 3, 4), (1, 2, 3, 4)]
assert star_p == [F(0)]+[F(1, 4)]*4
star_c = energy(star, star_p)
assert all(x >= star_c for x in mul(star, star_p))
assert all(s.Matrix(star)[:k, :k].det() > 0 for k in range(1, 6))
# Its magnitude is 1 + 4 tanh(L/2), using tanh(-log(q)/2) = (1-q)/(1+q).
assert sum(star_formula) == 1+4*(1-q0)/(1+q0)
print('FIDELITY: weighted-star formula and Theorem 4 output verified; energy =', star_c)

points = [(0, 0, 0, 0, 0), (1, 1, 0, 0, 1), (0, 0, 1, 0, 1),
          (0, 1, 0, 1, 1), (0, 0, 1, 1, 1)]
h = [[sum(x != y for x, y in zip(a, b)) for b in points] for a in points]
assert h == [[0, 3, 2, 3, 3], [3, 0, 3, 2, 4], [2, 3, 0, 3, 1],
             [3, 2, 3, 0, 2], [3, 4, 1, 2, 0]]
for i, j, k in product(range(5), repeat=3):
    assert h[i][i] == 0 and (h[i][j] == 0) == (i == j)
    assert h[i][j] == h[j][i] and h[i][j] <= h[i][k]+h[k][j]
    # The square-root embedding is sqrt(c) times these binary vectors.
    assert h[i][j] == sum((x-y)**2 for x, y in zip(points[i], points[j]))
print('GEOMETRY: all metric axioms and squared-distance/Hamming identities verified')

z0 = [[q0**d for d in row] for row in h]
history, p = deletion(z0)
w = [F(x, 558104619267) for x in
     [188194292000, 254397360000, -101751788000, -21051416000, 305206400000]]
v = [F(3192200, 11018319), F(4564000, 11018319), F(0), F(0), F(4564000, 11018319)]
assert history == [((0, 1, 2, 3, 4), w), ((0, 1, 4), v)]
assert list(s.Matrix(z0).inv()*s.ones(5, 1)) == [s.Rational(x.numerator, x.denominator) for x in w]
assert p == [F(15961, 61601), F(22820, 61601), F(0), F(0), F(22820, 61601)]
c0 = energy(z0, p)
assert c0 == F(11018319, 12320200)
assert mul(z0, p)[3] == F(438997299, 492808000)
delta0 = c0-mul(z0, p)[3]
assert delta0 == F(1735461, 492808000) > 0
pprime = [F(999, 1000)*x + (F(1, 1000) if i == 3 else 0) for i, x in enumerate(p)]
assert pprime == [F(x, 61601000) for x in [15945039, 22797180, 0, 61601, 22797180]]
assert sum(pprime) == 1 and all(x >= 0 for x in pprime)
improvement = energy(z0, p)-energy(z0, pprime)
assert improvement == F(1707687919, 246404000000000) > 0
print('WITNESS q=19/20: full weighting =', ', '.join(map(str, w)))
print('WITNESS: active sets =', [a for a, _ in history])
print('WITNESS: output =', ', '.join(map(str, p)))
print('WITNESS: C =', c0, '; C-(Zp)[3] =', delta0)
print('WITNESS: strict energy improvement =', improvement)

q = s.Symbol('q')
z = s.Matrix([[q**d for d in row] for row in h])
k = 2*q**8+2*q**6-3*q**4-3*q**2-1
ell = (1+q)*k
b = 1+q**2+2*q**4
numerators = s.Matrix([
    (1+q)*(q**5-q**4+2*q**3-2*q**2-1),
    (q**2+q+1)*(q**3-q**2-1),
    (q**2+1)*(q**4+2*q**3-q**2-1),
    2*q**6+2*q**4-q**2-q-1,
    -(q**2-q+1)*(2*q**2+q+1)])
wf = numerators/ell
for residual in z*wf-s.ones(5, 1):
    zero(residual)
zero(z.det()+(q-1)**4*(q+1)**4*k)
zero(k-(-2*q**2-1+q**2*(q**2-1)*(2*q**4+4*q**2+1)))
for i, sign in enumerate([-1, -1, 1, 1, -1]):
    lower = bernstein_positive(sign*numerators[i], q)
    print('FAMILY: numerator', i, 'signed Bernstein lower bound =', lower)
for i, sign in enumerate([1, 1, -1, -1, 1]):
    assert sign*wf[i].subs(q, s.Rational(19, 20)) > 0

vf = s.Matrix([1+q+q**2-q**3, 1+q+q**2, 0, 0, 1+q+q**2])/((1+q)*b)
retained = [0, 1, 4]
for i in retained:
    zero((z*vf)[i]-1)
    assert bernstein_positive(s.fraction(vf[i])[0], q) > 0
zero(z.extract(retained, retained).det()-(1-q**2)**2*(1+q**2)*b)
r3 = -(1-q)**2*(1+q)*(1+q+q**2)/b
zero((z*vf)[3]-1-r3)
print('FAMILY: determinants, two solve identities, signs and residual =', s.factor(r3))
removed = [2, 3]
cross = z.extract(removed, retained)
schur = z.extract(removed, removed)-cross*z.extract(retained, retained).inv()*cross.T
zero(schur[0, 1]-2*q**5*(q-1)*(q+1)/b)
print('MECHANISM: removed-block Schur off-diagonal =', s.factor(schur[0, 1]))

mass = sum(vf)
pf = vf/mass
cf = 1/mass
delta = cf-(z*pf)[3]
zero((pf.T*z*pf)[0]-cf)
zero(delta+r3/mass)
# Strict C<1 follows from positive mass on distinct points and off-diagonal q^H<1.
one_minus_c_num, one_minus_c_den = s.cancel((1-cf)/(1-q)).as_numer_denom()
if one_minus_c_den.subs(q, s.Rational(19, 20)) < 0:
    one_minus_c_num, one_minus_c_den = -one_minus_c_num, -one_minus_c_den
assert bernstein_positive(one_minus_c_num, q) > 0
assert bernstein_positive(one_minus_c_den, q) > 0
e = s.Symbol('epsilon')
e3 = s.eye(5)[:, 3]
pf_e = (1-e)*pf+e*e3
quad = (pf_e.T*z*pf_e)[0]
a = 1-cf+2*delta
zero(quad-(cf-2*e*delta+e**2*a))
zero((cf-quad).subs(e, delta/a)-delta**2/a)
print('FAMILY: exact energy expansion and epsilon=delta/(1-C+2delta) improvement verified')

# Scale zero: both the distance recursion and its terminal optimality certificate.
pzero = s.Matrix([s.Rational(1, 4), s.Rational(3, 8), 0, 0, s.Rational(3, 8)])
hm = s.Matrix(h)
assert hm*pzero == s.Matrix([s.Rational(9, 4), s.Rational(9, 4), 2,
                            s.Rational(9, 4), s.Rational(9, 4)])
history_zero, out_zero = deletion(h)
assert list(map(s.Rational, out_zero)) == list(pzero)
# Verify the conditional negative form symbolically for all zero-sum vectors.
a0, a1, a2, a3 = s.symbols('a0:4')
a_vec = s.Matrix([a0, a1, a2, a3, -a0-a1-a2-a3])
xmat = s.Matrix(points)
zero_sum_features = xmat.T*a_vec
zero((a_vec.T*hm*a_vec)[0]+2*(zero_sum_features.T*zero_sum_features)[0])
assert [s.limit(pf[i], q, 1, dir='-') for i in range(5)] == list(pzero)
zero(s.limit(r3/(1-q)**2, q, 1, dir='-')+s.Rational(3, 2))
print('SCALE ZERO: recursion is optimal here; positive-scale residual has coefficient -3/2')

# Confirm all 31 nonempty principal matrices are nonsingular at the exact witness.
for size in range(1, 6):
    for inds in combinations(range(5), size):
        assert s.Matrix(z0).extract(inds, inds).det() > 0
print('WITNESS: all 31 nonempty principal minors positive')
print('ALL_CHECKS_OK')

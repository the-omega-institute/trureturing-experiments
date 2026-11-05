#!/usr/bin/env python3
"""Exact checks for arXiv:2610.34388v1, Conjecture 11.1, item 6.

Requires Python 3 and SymPy. No floating-point arithmetic or randomized tests.
The tests support a preregistration; they are not a Lean proof of the conjecture.
"""
from itertools import product
import sympy as s


def eq(a, b):
    assert s.simplify(a - b) == 0, (a, b)


def matrix_eq(a, b):
    assert a.shape == b.shape
    for z in a - b:
        eq(z, 0)


def positive2(a):
    assert a == a.conjugate().T
    assert a[0, 0] > 0 and a.det() > 0


def sqrt2(a):
    """Principal square root of an exact positive definite Hermitian 2x2 matrix."""
    positive2(a)
    d = s.sqrt(a.det())
    r = (a + d * s.eye(2)) / s.sqrt(s.trace(a) + 2 * d)
    matrix_eq(r * r, a)
    positive2(r)
    return r


def geometric(a, b):
    h = sqrt2(a)
    hi = h.inv()
    return s.simplify(h * sqrt2(hi * b * hi) * h)


def spectral(a, b):
    h = sqrt2(geometric(a.inv(), b))
    return s.simplify(h * a * h)


# Fidelity: reproduce Proposition 5.43's displayed rank-one limiting quotient.
a, b, e = s.symbols('a b e', positive=True)
rho = s.diag(a, b)
v = s.Matrix([1, 1]) / s.sqrt(2)
x = rho.inv() * v
A = x * x.T + e * s.eye(2)
B0 = v * v.T
h = (v.T * A * v)[0]
P = B0 / s.sqrt(h)
matrix_eq(P * A * P, B0)
Phalf = B0 / h ** s.Rational(1, 4)
matrix_eq(Phalf * Phalf, P)
S0 = s.simplify(Phalf * A * Phalf)
matrix_eq(S0, s.sqrt(h) * B0)
ratio_sq = s.cancel(s.trace(rho * S0) ** 2 /
                    (s.trace(rho * A) * s.trace(rho * B0)))
eq(s.limit(ratio_sq, e, 0, dir='+'), (a + b) ** 2 / (4 * a * b))
print('FIDELITY: Proposition 5.43 rank-one squared-quotient limit = (a+b)^2/(4ab)')

# Evaluate the original means directly, without a Riccati-defined substitute.
A0 = s.Matrix([[5, 4], [4, 5]])
X0 = s.diag(1, 2)
B = X0 ** 2 * A0 * X0 ** 2
matrix_eq(geometric(A0.inv(), B), X0 ** 2)
matrix_eq(spectral(A0, B), X0 * A0 * X0)
print('MEAN FIDELITY: exact principal-root formulas verify a noncommuting Riccati pair')

# The algebraic identities in the proposed general upper-bound argument.
u, v0, t, x0, y0 = s.symbols('u v t x y', real=True)
F0 = u**2 + v0**2 + 2*t*u*v0
F2 = u**2*x0**2 + v0**2*y0**2 + 2*t*u*v0*x0*y0
F4 = u**2*x0**4 + v0**2*y0**4 + 2*t*u*v0*x0**2*y0**2
H = u**2*x0**2 + v0**2*y0**2 + t*u*v0*(x0**2+y0**2)
eq(F2, H-t*u*v0*(x0-y0)**2)
eq(F0*F4-H**2, (1-t**2)*u**2*v0**2*(x0**2-y0**2)**2)
r, j, aa, bb = s.symbols('r j a b', real=True)
cabs2 = r*r + j*j
delta_sq = ((aa-bb)**2 + 4*cabs2)/(aa+bb)**2
eq(delta_sq-cabs2/(aa*bb),
   (aa-bb)**2*(aa*bb-cabs2)/(aa*bb*(aa+bb)**2))
print('SYMBOLIC: both determinant identities and the off-diagonal/contrast identity hold')

# All small scalar cases, including zero components and repeated eigenvalues.
T = tuple(s.Rational(k, 4) for k in range(-3, 4))
UV = tuple(s.Integer(k) for k in range(4))
XY = (s.Rational(1, 4), s.Rational(1, 2), s.Integer(1), s.Integer(2), s.Integer(4))
count = 0
for tt, uu, vv, xx, yy in product(T, UV, UV, XY, XY):
    f0 = uu**2+vv**2+2*tt*uu*vv
    f2 = uu**2*xx**2+vv**2*yy**2+2*tt*uu*vv*xx*yy
    f4 = uu**2*xx**4+vv**2*yy**4+2*tt*uu*vv*xx**2*yy**2
    assert f0 >= 0 and f2 >= 0 and f4 >= 0
    assert (1-tt**2)*f2**2 <= f0*f4
    count += 1
print('EXACT SCALAR GRID:', count, 'cases satisfy (1-t^2) F2^2 <= F0 F4')

# Complex vectors and genuinely complex faithful densities, in an X eigenbasis.
vectors = [s.Matrix(z) for z in [(0,0),(1,0),(0,1),(1,1),(1,s.I),(s.I,1),(1,-1)]]
count = 0
for aa, bb, cc, xx, yy, z in product((1,2), (1,2), (0,s.Rational(1,2),s.I/2),
                                    XY[1:4], XY[1:4], vectors):
    R = s.Matrix([[aa,cc],[s.conjugate(cc),bb]])
    positive2(R)
    X = s.diag(xx, yy)
    q = lambda w: s.expand((w.conjugate().T*R*w)[0])
    ds = s.cancel(((aa-bb)**2+4*cc*s.conjugate(cc))/(aa+bb)**2)
    assert q(z) >= 0 and q(X*z) >= 0 and q(X**2*z) >= 0
    assert (1-ds)*q(X*z)**2 <= q(z)*q(X**2*z)
    count += 1
print('EXACT COMPLEX GRID:', count, 'vector/density/eigenvalue cases satisfy the squared bound')

# A non-diagonal complex X, obtained by an exact unitary conjugation.
U = s.Matrix([[3,4*s.I],[4*s.I,3]])/5
matrix_eq(U.conjugate().T*U, s.eye(2))
R = s.Matrix([[2,s.I/2],[-s.I/2,1]])
X = s.diag(s.Rational(1,2),2)
z = s.Matrix([1,s.I])
q = lambda R, z: s.expand((z.conjugate().T*R*z)[0])
for k in [0,1,2]:
    eq(q(U*R*U.conjugate().T, (U*X*U.conjugate().T)**k*(U*z)), q(R,X**k*z))
print('UNITARY CHECK: quadratic quantities agree for a non-diagonal complex X')

# The strictly positive sharpness family: symbolic traces, determinants and limit.
d = s.symbols('d', real=True)
R = s.Matrix([[1,-d],[-d,1]])/2
A = s.Matrix([[1+e,d],[d,d*d+e]])
X = s.diag(1,e)
B = X**2*A*X**2
S = X*A*X
eq(A.det(), e*(1+d*d+e))
eq(B.det(), e**4*A.det())
f0 = 1-d*d+2*e
f2 = 1+e+(d*d+e)*e**2-2*d*d*e
f4 = 1+e+(d*d+e)*e**4-2*d*d*e**2
for M, f in [(A,f0),(S,f2),(B,f4)]:
    eq(2*s.trace(R*M), f)
matrix_eq(X**2*A*X**2, B)
eq(s.limit(f2**2/(f0*f4),e,0,dir='+'),1/(1-d*d))
print('SHARPNESS: exact trace polynomials and squared-quotient limit = 1/(1-d^2)')
print('SCOPE: exact algebra and finite checks; no kernel verification or exhaustive proof of all matrices')
print('ALL_CHECKS_OK')

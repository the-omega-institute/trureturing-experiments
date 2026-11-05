#!/usr/bin/env python3
"""Exact algebraic checks for arXiv:2609.40355v2, Section VII's question.

Requires Python 3 and SymPy 1.14.0. No floating-point input or inference.
Finite identities and symbolic certificates support the accompanying argument;
this program does not enumerate all mixed-unitary channels or certify an exponent.
"""
import itertools
import sympy as s

R = s.Rational
I = s.I
kron = s.kronecker_product
E = s.eye(3)
E9 = s.eye(9)
e = [E[:, i] for i in range(3)]
P = s.diag(1, 1, 0)
Q = s.diag(0, 0, 1)

def projector(v):
    return v * v.conjugate().T

def equal(A, B):
    assert A.shape == B.shape
    assert all(s.expand(x) == 0 for x in A - B)

def unit(i, j):
    return e[i] * e[j].T

def wh(X):
    return (s.trace(X) * E - X.T) / 2

def two_use(images, X):
    out = s.zeros(9)
    for i, j, k, l in itertools.product(range(3), repeat=4):
        coefficient = X[3*i+k, 3*j+l]
        if coefficient != 0:
            out += coefficient * kron(images[i, j], images[k, l])
    return out

def pt_left(X):
    return s.Matrix(9, 9, lambda a, b:
                    X[3*(b//3)+a%3, 3*(a//3)+b%3])

def marginal_left(X):
    return s.Matrix(3, 3, lambda i, j: sum(X[3*i+k, 3*j+k] for k in range(3)))

def marginal_right(X):
    return s.Matrix(3, 3, lambda k, l: sum(X[3*i+k, 3*i+l] for i in range(3)))

psi = (kron(e[0], e[1]) - kron(e[1], e[0])) / s.sqrt(2)
rho = projector(psi)
wh_images = {(i, j): wh(unit(i, j)) for i, j in itertools.product(range(3), repeat=2)}

# Channel fidelity: the three antisymmetric Kraus operators realize N exactly.
kraus = [(unit(i, j) - unit(j, i))/s.sqrt(2) for i in range(3) for j in range(i+1, 3)]
equal(sum((K.conjugate().T*K for K in kraus), s.zeros(3)), E)
for i, j in itertools.product(range(3), repeat=2):
    X = unit(i, j)
    equal(sum((K*X*K.conjugate().T for K in kraus), s.zeros(3)), wh(X))
print('Fidelity: antisymmetric Kraus realization equals N on all 9 matrix units.')

# Source Theorem 6, its displayed output and the eq:entcomp compression.
output = two_use(wh_images, rho)
source_output = (E9-kron(E, P/2)-kron(P/2, E)+rho)/4
equal(output, source_output)
PP = kron(P, P)
equal(PP*output*PP, rho/4)
equal(marginal_left(rho), P/2)
equal(marginal_left(output), P/4+Q/2)
equal(marginal_right(rho), P/2)
equal(marginal_right(output), P/4+Q/2)
assert s.trace(output) == 1
assert output.eigenvals() == {R(1, 4): 2, R(1, 8): 4, s.Integer(0): 3}
print('Fidelity: Theorem 6 output, both marginals, and P tensor P compression reproduced exactly.')

# Proposed separable probe, with an explicit eleven-term decomposition.
sigma = (E9+2*rho)/11
local = [(e[0]+e[1])/s.sqrt(2), (e[0]-e[1])/s.sqrt(2),
         (e[0]+I*e[1])/s.sqrt(2), (e[0]-I*e[1])/s.sqrt(2), e[0], e[1]]
vectors = [kron(local[a], local[a^1]) for a in range(6)]
vectors += [kron(e[i], e[j]) for i, j in [(0, 2), (1, 2), (2, 0), (2, 1), (2, 2)]]
assert len(vectors) == 11
for v in vectors:
    equal(v.conjugate().T*v, s.ones(1))
equal(sum((projector(v) for v in vectors), s.zeros(9))/11, sigma)
assert list(11*sigma.diagonal()) == [1, 2, 1, 2, 1, 1, 1, 1, 1]
assert sigma[1, 3] == sigma[3, 1] == -R(1, 11)
assert s.trace(sigma) == 1
assert sigma.eigenvals() == {R(3, 11): 1, R(1, 11): 8}
print('Settlement certificate: eleven unit product projectors average to (I9+2 rho)/11.')
print('Probe spectrum: 3/11 (multiplicity 1), 1/11 (multiplicity 8).')

# Generic product-vector obstruction, verified as polynomial identities.
a = s.Matrix(s.symbols('a0:2'))
b = s.Matrix(s.symbols('b0:2'))
c = s.Matrix(s.symbols('c0:2'))
d = s.Matrix(s.symbols('d0:2'))
A, B = a*b.T, c*d.T
psi2 = s.Matrix([0, 1, -1, 0])/s.sqrt(2)
phi = kron(A, B)*psi2
coefficient_matrix = s.Matrix(2, 2, list(phi))
assert s.expand(coefficient_matrix.det()) == 0
assert s.Matrix(2, 2, list(psi2)).det() == R(1, 2)
assert min(pt_left(rho/4).eigenvals()) == -R(1, 8)
x, y, u, v = s.symbols('x y u v', real=True)
z = s.Symbol('z')
bb = s.Matrix([x+I*y, u+I*v])
C = s.eye(2)-bb*bb.conjugate().T
normsq = x*x+y*y+u*u+v*v
assert s.expand(C.det()-(1-normsq)) == 0
assert s.expand(s.trace(C)-(2-normsq)) == 0
print('Obstruction certificates: rank-one local images have zero coefficient determinant; psi determinant = 1/2.')
print('Compressed target partial-transpose minimum eigenvalue: -1/8.')
print('Unitary-block certificates: det(I2-b b*)=1-||b||^2 and trace=2-||b||^2.')

# Entire parameter family: symbolic, not interpolation from sample points.
t = s.symbols('t', real=True)
sigma_t = (1-t)*E9/9+t*rho
alpha = 11*t/2
equal(sigma_t, alpha*sigma+(1-alpha)*E9/9)
assert s.solve_univariate_inequality(alpha >= 0, t, relational=False) == s.Interval(0, s.oo)
assert s.solve_univariate_inequality(alpha <= 1, t, relational=False) == s.Interval(-s.oo, R(2, 11))
lambda_min = (2-11*t)/18
assert s.expand(pt_left(sigma_t).charpoly(z).as_expr() -
    (z-lambda_min)*(z-((1-t)/9+t/2))**3*(z-(1-t)/9)**5) == 0
print('Family: sigma_t=(11t/2) sigma+(1-11t/2) I9/9 for 0 <= t <= 2/11.')
print('Partial-transpose eigenvalues: (2-11t)/18 once, (2+7t)/18 three times, (1-t)/9 five times.')

# Arbitrary unital complex-linear map: eight free 3x3 images, the ninth fixed
# by L(I3)=I3. All 72 scalar variables remain independent.
images = {}
for i, j in itertools.product(range(3), repeat=2):
    if (i, j) != (2, 2):
        images[i, j] = s.Matrix(3, 3, s.symbols('f%d%d_0:9' % (i, j)))
images[2, 2] = E-images[0, 0]-images[1, 1]
equal(sum((images[i, i] for i in range(3)), s.zeros(3)), E)
gap_rho = two_use(images, rho)-output
gap_t = two_use(images, sigma_t)-two_use(wh_images, sigma_t)
equal(gap_t, t*gap_rho)
equal(two_use(images, sigma)-two_use(wh_images, sigma), R(2, 11)*gap_rho)
print('Symbolic dilution identity: (L tensor L-N tensor N)(sigma_t)=t (L tensor L-N tensor N)(rho), for every unital L.')

# Strict full-rank lower bound for the constructed endpoint's target output.
target_sigma = two_use(wh_images, sigma)
assert min(target_sigma.eigenvals()) >= R(1, 11)
print('Target endpoint output is full rank, with minimum eigenvalue 1/11.')
print('Scope: exact algebraic certificates; compactness, Pinsker and the operational Stein bridge are analytic proof obligations.')
print('ALL_CHECKS_OK')

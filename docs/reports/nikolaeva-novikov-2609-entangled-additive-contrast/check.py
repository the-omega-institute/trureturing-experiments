#!/usr/bin/env python3
"""Exact source-fidelity and symbolic-family checks; no numerical approximations."""
import itertools
import sys
import sympy as s

checks = 0

def equal(name, lhs, rhs):
    global checks
    delta = lhs - rhs
    entries = list(delta) if isinstance(delta, s.MatrixBase) else [delta]
    assert all(s.simplify(x) == 0 for x in entries), (name, delta)
    assert not any(x.has(s.Float) for x in entries), name
    checks += 1
    print('OK:', name)


def psd_certificate(name, M, parameter=None):
    """All principal minors nonnegative; for t>=0 use rational nonnegative coefficients."""
    global checks
    equal(name + ' Hermitian', M, M.H)
    for count in range(1, M.rows + 1):
        for indices in itertools.combinations(range(M.rows), count):
            minor = s.factor(M.extract(indices, indices).det())
            if parameter is None:
                assert minor.is_nonnegative is True, (name, indices, minor)
            else:
                poly = s.Poly(s.expand(minor), parameter, domain=s.QQ)
                assert all(c >= 0 for c in poly.all_coeffs()), (name, indices, poly)
            assert not minor.has(s.Float)
    checks += 1
    print('OK:', name, '(all principal minors exact)')


def pt(M):
    return s.Matrix(4, 4, lambda i, j: M[2*(i//2) + j%2, 2*(j//2) + i%2])

# Fidelity: Theorem 4.4 and Section 5.9 of arXiv:2609.06997v2.
X = s.diag(4, 3, 1, 1)
Y = s.Matrix([[4, 0, 0, 0], [0, 2, 1, 0], [0, 1, 2, 0], [0, 0, 0, 1]])
Z = X*Y
z = s.symbols('z')
expected_charpoly = (z-4)*(z-3)*(z-1)**2
equal('source X spectrum (4,3,1,1)', X.charpoly(z).as_expr(), expected_charpoly)
equal('source Y spectrum (4,3,1,1)', Y.charpoly(z).as_expr(), expected_charpoly)
equal('source product entries', Z, s.Matrix([[16,0,0,0],[0,6,3,0],[0,1,2,0],[0,0,0,1]]))
psd_certificate('source product singular values >=1', Z.H*Z-s.eye(4))
psd_certificate('source product singular values <=16', 256*s.eye(4)-Z.H*Z)
equal('source attained singular extrema', (Z.H*Z).extract([0,3],[0,3]), s.diag(256,1))
equal('source rho_Y partial transpose', pt(Y/9), s.Matrix([[4,0,0,1],[0,2,0,0],[0,0,2,0],[1,0,0,1]])/9)
pt_spectrum = ((z-s.Rational(2,9))**2 * (z-(5+s.sqrt(13))/18) * (z-(5-s.sqrt(13))/18))
equal('source rho_Y partial-transpose spectrum', pt(Y/9).charpoly(z).as_expr(), pt_spectrum)
psd_certificate('source rho_Y is PPT', pt(Y/9))
assert X*Y != Y*X
print('SOURCE_FIDELITY_OK: kappa(X)=kappa(Y)=4; kappa(XY)=16; rho_Y is PPT')

# Family: k = 3+t, with t>0; the exact certificates below hold throughout this domain.
t = s.symbols('t', positive=True)
k = 3+t
b = (k-1)/2
T = s.Rational(3,2)*(k+1)
e = [s.eye(4)[:,i] for i in range(4)]
p, q = (e[0]+e[3])/s.sqrt(2), (e[0]-e[3])/s.sqrt(2)
u, v = (e[1]+e[2])/s.sqrt(2), (e[1]-e[2])/s.sqrt(2)
P = lambda x: x*x.H
U = s.Matrix.hstack(p,q,u,v)
V = s.Matrix.hstack(p,q,e[1],e[2])
D = s.diag(k,1,b,1)
A = k*P(p)+P(q)+b*P(u)+P(v)
B = k*P(p)+P(q)+b*P(e[1])+P(e[2])
rho, sigma = A/T, B/T
for name, basis in [('A', U), ('B', V)]:
    equal(name + ' orthonormal eigenbasis', basis.H*basis, s.eye(4))
for name, M, basis in [('A', A, U), ('B', B, V)]:
    equal(name + ' diagonalization', basis.H*M*basis, D)
    equal(name + ' Hermitian', M.H, M)
    equal(name + ' trace normalization', s.trace(M), T)
    equal(name + ' top eigenvector', M*p, k*p)
    equal(name + ' bottom eigenvector', M*q, q)
for expression in [k, b, k-b, b-1, T]:
    assert expression.is_positive is True, expression
print('OK: k>b>1 and T>0 for every t>0')
for name, M in [('rho',rho),('sigma',sigma)]:
    equal(name + ' trace one', s.trace(M), 1)
    equal(name + ' negative PT expectation', (v.H*pt(M)*v)[0], -t/(6*(t+4)))
assert (-t/(6*(t+4))).is_negative is True
commutator = (u.H*(rho*sigma-sigma*rho)*v)[0]
equal('nonzero commutator coefficient', commutator, t**2/(18*(t+4)**2))
assert (t**2/(18*(t+4)**2)).is_positive is True
M = (s.diag(b,1)*s.Matrix([[b+1,b-1],[b-1,b+1]])/2)
expected_block = s.zeros(4)
expected_block[0,0], expected_block[1,1] = k*k, 1
expected_block[2:4,2:4] = M
equal('AB block decomposition in Bell basis', U.H*A*B*U, expected_block)
G = s.simplify(M.H*M)
psd_certificate('family middle singular values >=1', G-s.eye(2), t)
psd_certificate('family middle singular values <=k^2', k**4*s.eye(2)-G, t)
# Together with the two scalar blocks, the singular extrema of AB are exactly k^2 and 1.
# Positive rescaling leaves their ratio unchanged, so kappa(rho sigma)=k^2.
equal('positive normalization cancels in condition number', (k/T)/(1/T), k)
equal('positive product normalization cancels in condition number', (k**2/T**2)/(1/T**2), k**2)
equal('OLC additivity on the positive domain', s.log(k**2)/2, s.log(k)/2+s.log(k)/2)
assert all(not x.has(s.Float) for mat in [rho,sigma,A,B,G] for x in mat)
print('FAMILY_OK: every real k>3; both states NPT; noncommuting; kappa=k; kappa(product)=k^2')
print('EXACT_CHECK_COUNT:', checks)
print('PYTHON:', sys.version.split()[0], 'SYMPY:', s.__version__)
print('ALL_CHECKS_OK')

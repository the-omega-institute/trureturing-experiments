#!/usr/bin/env python3
"""Exact checks for Kanachova--Liu--Tu, arXiv:2609.11151v1, Appendix 7.

SymPy proves polynomial identities by exact expansion. Exhaustive finite checks
use Gaussian integers as integer pairs; 2*M is integral. Dividing a nonzero
input by sqrt(N) divides the checked expectation and Gram determinant by N**2.
No floating point, numerical diagonalization, or random sampling is used.
"""
from fractions import Fraction
from itertools import product
import sympy as sp


def pt(matrix):
    """Qubit-major indexing, transpose only the qutrit factor."""
    return sp.Matrix(6, 6, lambda i, j: matrix[3*(i//3)+j%3, 3*(j//3)+i%3])


U = sp.Matrix([[0, 0, -1], [0, 1, 0], [-1, 0, 0]])
K = sp.kronecker_product(sp.eye(2), U)
assert U.T == U and U.H*U == sp.eye(3)


def reconstruct(psi):
    rho = psi*psi.H
    return (rho + K*pt(rho)*K.H)/2


# Fidelity: Proposition 1 in its explicitly specified meson Schmidt basis.
# The characteristic polynomial exposes -a*b/2 and a positive 2x2 block.
a, b, t = sp.symbols('a b t', positive=True)
psi = sp.Matrix([a, 0, 0, 0, b, 0])
M = reconstruct(psi)
expected_charpoly = t*(t-a*a/2)*(t*t-a*a*b*b/4)*(t*t-(a*a/2+b*b)*t+a*a*b*b/4)
assert sp.expand(pt(M).charpoly(t).as_expr().subs(sp.Symbol('t'), t)-expected_charpoly) == 0
# The last block has positive trace and determinant for a,b>0, so both roots
# are positive. The only negative eigenvalue is therefore -a*b/2.
for aa, bb in [(sp.Rational(3, 5), sp.Rational(4, 5)),
               (sp.sqrt(2)/2, sp.sqrt(2)/2)]:
    p = psi.subs({a: aa, b: bb})
    assert (p.H*p)[0] == 1
    eigs = pt(reconstruct(p)).eigenvals()
    negative = [x for x in eigs if x.is_negative]
    assert negative == [-aa*bb/2]
    print(f'FIDELITY a={aa}, b={bb}: lambda_min(measured_PT)={negative[0]}')
print('FIDELITY symbolic Proposition 1 characteristic polynomial: OK')

# Arbitrary complex coordinates are represented by twelve independent reals.
x = sp.symbols('x0:12', real=True)
p = sp.Matrix([x[2*j] + sp.I*x[2*j+1] for j in range(6)])
u, v = p[:3, :], p[3:, :]
r, s = U*u.conjugate(), U*v.conjugate()
eta = s.col_join(-r)
N = (p.H*p)[0]
delta = (u.H*u)[0]*(v.H*v)[0] - (u.H*v)[0]*(v.H*u)[0]
minor_sum = sum(sp.expand_complex(u[j]*v[k]-u[k]*v[j]) *
                sp.conjugate(sp.expand_complex(u[j]*v[k]-u[k]*v[j]))
                for j in range(3) for k in range(j+1, 3))
assert sp.expand((p.H*eta)[0]) == 0
assert sp.expand((eta.H*eta)[0]-N) == 0
assert sp.expand((eta.H*reconstruct(p)*eta)[0]+delta) == 0
assert sp.expand(delta-minor_sum) == 0
assert sp.simplify(reconstruct(p)-reconstruct(p).H) == sp.zeros(6)
print('SYMBOLIC orthogonality, norm, expectation=-Gram determinant, minors, Hermiticity: OK')

# Independent implementation: exact Gaussian-integer matrix operations.
ZERO = (0, 0)
def add(z, w): return (z[0]+w[0], z[1]+w[1])
def neg(z): return (-z[0], -z[1])
def sub(z, w): return add(z, neg(w))
def conj(z): return (z[0], -z[1])
def mul(z, w): return (z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0])
def abs2(z): return z[0]*z[0]+z[1]*z[1]
def gsum(zs):
    result = ZERO
    for z in zs: result = add(result, z)
    return result
def inner(z, w): return gsum(mul(conj(a), b) for a, b in zip(z, w))

alphabet = [ZERO, (1, 0), (-1, 0), (0, 1), (0, -1)]
permutation, signs = [2, 1, 0, 5, 4, 3], [-1, 1, -1, -1, 1, -1]
counts = {'zero': 0, 'product': 0, 'entangled': 0}
for p in product(alphabet, repeat=6):
    n = sum(map(abs2, p))
    if n == 0:
        counts['zero'] += 1
        continue
    u, v = p[:3], p[3:]
    r = tuple(mul((signs[j], 0), conj(u[permutation[j]])) for j in range(3))
    s = tuple(mul((signs[j], 0), conj(v[permutation[j]])) for j in range(3))
    e = s+tuple(map(neg, r))
    rho = [[mul(p[i], conj(p[j])) for j in range(6)] for i in range(6)]
    rho_pt = [[rho[3*(i//3)+j%3][3*(j//3)+i%3] for j in range(6)] for i in range(6)]
    twice_M = [[add(rho[i][j], mul((signs[i]*signs[j], 0),
                      rho_pt[permutation[i]][permutation[j]]))
                for j in range(6)] for i in range(6)]
    q = gsum(mul(conj(e[i]), mul(twice_M[i][j], e[j]))
             for i in range(6) for j in range(6))
    d = sum(map(abs2, u))*sum(map(abs2, v))-abs2(inner(u, v))
    minors = [sub(mul(u[j], v[k]), mul(u[k], v[j]))
              for j in range(3) for k in range(j+1, 3)]
    assert d == sum(map(abs2, minors)) >= 0
    assert q == (-2*d, 0)
    assert inner(p, e) == ZERO and sum(map(abs2, e)) == n
    # Decide product structure independently by exact cross multiplication.
    pivot = next((j for j in range(3) if u[j] != ZERO), None)
    is_product = pivot is None or all(mul(u[pivot], v[j]) == mul(v[pivot], u[j])
                                          for j in range(3))
    assert is_product == (d == 0)
    if is_product:
        counts['product'] += 1
    else:
        counts['entangled'] += 1
        assert Fraction(q[0], 2*n*n) == -Fraction(d, n*n) < 0
print(f'EXHAUSTIVE alphabet={{0,1,-1,i,-i}}, six coordinates: {counts}')
assert sum(counts.values()) == 5**6
print('ALL_CHECKS_OK')

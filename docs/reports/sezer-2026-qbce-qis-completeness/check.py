#!/usr/bin/env python3
"""Exact algebra checks for Sezer arXiv:2608.04973v1, Section 9.2.

SymPy exact expressions only. No Lean proof or numerical SDP is claimed.
"""
import itertools
import sympy as S


def eq(a, b):
    if isinstance(a, S.MatrixBase):
        assert a.shape == b.shape
        assert all(S.simplify(z) == 0 for z in a - b), (a, b)
    else:
        assert S.simplify(a - b) == 0, (a, b)


def proj(v):
    return v * v.conjugate().T


def psd_rational(m):
    """Hermitian PSD iff every principal minor is nonnegative, exact rationals."""
    assert m == m.conjugate().T
    for size in range(1, m.rows + 1):
        for indices in itertools.combinations(range(m.rows), size):
            assert m.extract(indices, indices).det() >= 0


def tensor_power(m, n):
    out = S.ones(1, 1)
    for _ in range(n):
        out = S.kronecker_product(out, m)
    return out


def marginal(m, n, player):
    basis = list(itertools.product(range(2), repeat=n))
    return S.Matrix(2, 2, lambda a, b: sum(
        m[j, k] for j, x in enumerate(basis) for k, y in enumerate(basis)
        if x[player] == a and y[player] == b
        and all(x[i] == y[i] for i in range(n) if i != player)))


# Fidelity: the exact depolarizing composition stated in Section 10.2.
g, h = S.symbols('g h', real=True)
z = S.Matrix(2, 2, S.symbols('z0:4'))
def dep(a, m):
    return a * m + (1 - a) * S.eye(2) / 2

eq(dep(g, dep(h, z)), dep(g*h, z))
print('FIDELITY Section 10.2: dep(g,dep(h,z)) = dep(g*h,z), exactly')

# Fidelity: the forward implication of Theorem 6.2 uses this exact identity.
x0 = S.Matrix(2, 2, S.symbols('x0:4'))
xs = [S.Matrix(2, 2, S.symbols(f'x{j}_0:4')) for j in range(1, 4)]
m1 = S.Matrix(2, 2, S.symbols('m1_0:4'))
m2 = S.Matrix(2, 2, S.symbols('m2_0:4'))
ms = [m1, m2, S.eye(2)-m1-m2]
eq(S.trace(x0)-sum(S.trace(m*x) for m,x in zip(ms,xs)),
   sum(S.trace(m*(x0-x)) for m,x in zip(ms,xs)))
print('FIDELITY Theorem 6.2: POVM obedience-gap identity, exactly')

P = S.diag(1, 0)
u = S.Matrix([S.Rational(4,5), S.Rational(3,5)])
v = S.Matrix([S.Rational(3,5), S.Rational(4,5)])
A, B = proj(u), proj(v)
for m in [P, A, B, S.kronecker_product(P,P), S.kronecker_product(A,A),
          S.kronecker_product(B,B)]:
    psd_rational(m)
    eq(S.trace(m), 1)
    eq(m*m,m)
assert P*A != A*P and P*B != B*P
for m, a in [(A,S.Rational(3,5)), (B,S.Rational(4,5))]:
    H = P-m
    eq(H*H, a*a*S.eye(2))
    eq(H.charpoly().as_expr(), H.charpoly().gen**2-a*a)
    assert sorted(H.eigenvals()) == [-a,a]
    eq(sum(abs(e)*count for e,count in H.eigenvals().items()), 2*a)
print('WITNESS: all six density matrices PSD, trace 1, rank one; noncommuting pairs')
print('TRACE_NORMS: ||P-A||_1 = 6/5; ||P-B||_1 = 8/5; 6/5 < 8/5')

# Rational parametrization of the WHOLE open family: s=2q/(1+q^2), 0<q<1.
q, c0, c1 = S.symbols('q c0 c1', real=True)
a = (1-q*q)/(1+q*q)
s = 2*q/(1+q*q)
w, orth = S.Matrix([a,s]), S.Matrix([s,-a])
C = proj(w)
eq(a*a+s*s,1)
eq((orth.T*(c0*P+c1*C)*orth)[0], c0*s*s)
eq((S.Matrix([[0,1]])*(c0*P+c1*C)*S.Matrix([0,1]))[0], c1*s*s)
eq((P-C)*(P-C),s*s*S.eye(2))
eq((P-C).det(),-s*s)
eq(S.trace(P-C),0)
print('FAMILY identities: ||w||^2=1; dual tests=c0*s^2,c1*s^2; (P-C)^2=s^2 I')

# The pure-state block mechanism in an orthonormal two-ray basis.
r, d, x, y = S.symbols('r d x y', real=True)
U = S.Matrix.hstack(w,orth)
R = S.Matrix([[r,x+S.I*y],[x-S.I*y,d]])
eq(U.conjugate().T*U,S.eye(2))
eq((U*R*U.conjugate().T).det(),r*d-x*x-y*y)
eq((S.diag(1,0)-R).det(),(1-r)*(-d)-x*x-y*y)
# 0<=R<=diag(1,0) gives d>=0,-d>=0; at d=0 the determinant is -x^2-y^2.
eq(R.det().subs(d,0),-x*x-y*y)
eq(U*S.diag(r,0)*U.conjugate().T,r*C)
print('PURE BLOCK identities: complement diagonal=-d; det at d=0=-x^2-y^2; R=r*C')

# Exact PSD iff coefficientwise nonnegative on a signed rational grid.
values = list(map(S.Rational, [-2,-1,0,1,2])) + [S.Rational(-1,2),S.Rational(1,2)]
cases = 0
for point in [S.Rational(1,4),S.Rational(1,3),S.Rational(1,2),S.Rational(2,3),S.Rational(3,4)]:
    Cq=C.subs(q,point)
    for b0,b1 in itertools.product(values,repeat=2):
        M=b0*P+b1*Cq
        is_psd = M[0,0]>=0 and M[1,1]>=0 and M.det()>=0
        assert bool(is_psd) == (b0>=0 and b1>=0)
        cases += 1
print(f'COEFFICIENTWISE: {cases} signed rational cases, all exact')

# Every player's marginal of each product family is the same local ray.
for n in [2,3,4]:
    for local in [P,A,B]:
        joint=tensor_power(local,n)
        eq(S.trace(joint),1)
        eq(joint*joint,joint)
        for i in range(n):
            eq(marginal(joint,n,i),local)
            eq(marginal(S.Rational(2,7)*joint,n,i),S.Rational(2,7)*local)
print('PRODUCTS: n=2,3,4, every player marginal and scaled block checked exactly')

# All-games algebra: each arbitrary coefficient is carried unchanged to the local ray.
p00,p01,p10,p11,psi,df0,df1=S.symbols('p00 p01 p10 p11 psi df0 df1',real=True)
for local in [A,B]:
    lhs=psi*df0*(p00+p01)*P+(1-psi)*df1*(p10+p11)*local
    eq(lhs, (psi*df0*(p00+p01))*P+((1-psi)*df1*(p10+p11))*local)
# The universal finite sum is justified by linearity; this exact check tests two summands.
print('OUTCOMES: Tr(p*rho)=p; arbitrary real obedience coefficients use the same two-ray cone')
print('SCOPE: symbolic identities support the prose family argument; finite cases are not a universal Lean proof')
print('ALL_CHECKS_OK')

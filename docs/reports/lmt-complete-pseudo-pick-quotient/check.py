#!/usr/bin/env python3
"""Exact algebra checks for arXiv:2607.25088v2 Q9.3; Python 3 + SymPy.
Universal interpolation is proved in the preregistration outline, not by sampling.
"""
from fractions import Fraction as F
from itertools import combinations
import sympy as S


def zero(a):
    assert all(S.cancel(x) == 0 for x in a), a


def psd(a):
    """Exact symmetric rational elimination, including zero pivots."""
    a = [[F(x) for x in row] for row in a]
    n = len(a)
    assert all(a[i][j] == a[j][i] for i in range(n) for j in range(n))
    for i in range(n):
        pivot = a[i][i]
        if pivot < 0:
            return False
        if pivot == 0:
            if any(a[i][j] for j in range(i + 1, n)):
                return False
            continue
        for j in range(i + 1, n):
            for k in range(j, n):
                a[j][k] -= a[j][i]*a[i][k]/pivot
                a[k][j] = a[j][k]
    return True


# Source fidelity: Example 8.3, label eg:recursivel. Only the exact displayed
# decomposition and a0,a1,a2 are reproduced; no approximate roots are evidence.
u = S.symbols('u')
Q = 6+2*u+3*u**3+2*u**4
kinv = 2*(3-2*u-u**2)+u**3*(3-u-2*u**2)
assert S.expand(kinv-(1-u)*Q) == 0
H = (7+5*u+5*u**2+2*u**3)/(13*Q)
assert S.cancel(1/kinv-1/(13*(1-u))-H) == 0
series = S.series(1/kinv,u,0,3).removeO().expand()
coef = [series.coeff(u,i) for i in range(3)]
assert coef == [S.Rational(1,6),S.Rational(1,9),S.Rational(7,54)]
print('FIDELITY Example 8.3: exact decomposition; a0,a1,a2 =',coef)

# Family and negative form as symbolic identities, not sampled decimals.
e = S.symbols('epsilon',real=True)
x = S.Matrix([e,1,2])
one = S.ones(3,1)
K = one*one.T+x*x.T
L = S.ones(3)+S.diag(e**2,1,4)
R = S.Matrix(3,3,lambda i,j:L[i,j]/K[i,j])
v = S.Matrix([-2,1,1])
q = S.factor((v.T*R*v)[0])
assert S.cancel(q-(S.Rational(20,3)-4/(1+e)-4/(1+2*e))) == 0
assert S.cancel(q-4*(10*e**2+6*e-1)/(3*(1+e)*(1+2*e))) == 0
upper = (S.sqrt(19)-3)/10
assert S.simplify(10*upper**2+6*upper-1) == 0
assert S.Rational(1,10)<upper
assert q.subs(e,S.Rational(1,10)) == -S.Rational(10,33)
print('QUOTIENT q(epsilon) =',q)
print('QUOTIENT q(1/10) =',q.subs(e,S.Rational(1,10)))
print('INTERVAL: denominators positive; numerator polynomial strictly increasing on epsilon>0; positive root =',upper)

# Two relations force a common value entrywise for arbitrary matrix size.
t = S.Matrix([-1,2-e,e-1])
assert (t.T*one)[0] == 0
assert S.expand((t.T*x)[0]) == 0
relations = S.Matrix([list(t),list(S.matrix_multiply_elementwise(t,x))])
assert relations.rank() == 2
assert relations.nullspace() == [one]
for j,k in combinations(range(3),2):
    assert S.factor(K.extract([j,k],[j,k]).det()-(x[j]-x[k])**2) == 0
print('RIGIDITY: [t;t*x] has rank 2 and nullspace span(1,1,1)')

# Symbolic complex block identity. The written argument uses it for every m.
for m in [1,2]:
    W = [S.Matrix(m,m,lambda a,b:S.Symbol(f'w{j}_{a}_{b}')) for j in range(3)]
    P = S.BlockMatrix([[K[i,j]*(S.eye(m)-W[i]*W[j].conjugate().T)
                       for j in range(3)] for i in range(3)]).as_explicit()
    T = S.kronecker_product(t,S.eye(m))
    A = sum((t[i]*W[i] for i in range(3)),S.zeros(m))
    B = sum((t[i]*x[i]*W[i] for i in range(3)),S.zeros(m))
    zero(T.T*P*T+A*A.conjugate().T+B*B.conjugate().T)
print('NULL FORM: T* P T = -A A* - B B* exact for unrestricted complex symbolic m=1,2')

# Anchor compression and reconstruction, scalar coefficient identities which
# tensor with any matrix size when the three core data agree.
lam = S.Matrix([[1/(1-e),-e/(1-e),0]])
assert S.cancel((lam*one)[0]-1) == 0
assert S.cancel((lam*x)[0]) == 0
for r in range(7):
    features = S.zeros(3+r,2+r)
    for i in range(3+r):
        features[i,0] = 1
        if i<3:
            features[i,1] = x[i]
        else:
            features[i,i-1] = 1
    fullK = features*features.T
    C = S.zeros(1+r,3+r)
    C[0,:3] = lam
    for j in range(r):
        C[1+j,3+j] = 1
    anchorK = S.ones(1+r)+S.diag(0,*([1]*r))
    zero(C*fullK*C.T-anchorK)
    duplicate = S.zeros(3+r,1+r)
    for i in range(3):
        duplicate[i,0] = 1
    for j in range(r):
        duplicate[3+j,1+j] = 1
    d = [e**2,S.Integer(1),S.Integer(4)]+[1+2**(j+1) for j in range(1,r+1)]
    fullL = S.ones(3+r)+S.diag(*d)
    correction = S.diag(e**2,1,4,*[2**(j+1) for j in range(1,r+1)])
    zero(fullL-duplicate*anchorK*duplicate.T-correction)
print('ANCHOR: compression and diagonal reconstruction exact for 0..6 extra nodes')

# Exact finite checks include nonconsecutive extra labels; positivity here
# uses Fraction elimination, not floating eigenvalues or numerical tolerances.
count = 0
for eps in [F(1,100),F(1,20),F(1,10),F(1,8)]:
    core = [eps,F(1),F(2)]
    extras = [1,2,4,7,10]
    d = [a*a for a in core]+[F(1+2**(j+1)) for j in extras]
    N = 3+len(extras)
    fk = [[F(1)+(core[i]*core[j] if i<3 and j<3 else F(i==j and i>=3))
           for j in range(N)] for i in range(N)]
    fl = [[F(1)+(d[i] if i==j else F(0)) for j in range(N)] for i in range(N)]
    assert psd(fk) and psd(fl)
    assert S.Matrix(fk).rank() == 2+len(extras)
    assert S.Matrix(fl).det()>0
    fq = [[fl[i][j]/fk[i][j] for j in range(N)] for i in range(N)]
    assert not psd([row[:3] for row in fq[:3]])
    cpkernel = [[F(1)-F(1)/fl[i][j] for j in range(N)] for i in range(N)]
    assert psd(cpkernel)
    assert all(0<cpkernel[i][i]<1 for i in range(N))
    assert sum((F(1,2**j) for j in extras),F(0))<=1
    for size in range(1,N+1):
        for inds in combinations(range(N),size):
            if sum(i<3 for i in inds)>2:
                continue
            assert psd([[fq[i][j] for j in inds] for i in inds])
            count += 1
    for a,b in combinations(core,2):
        rho = F(1)/(1+a*b)
        sigma = 2/(1+rho)
        assert 0<rho<1 and 1<sigma<2
print('FINITE CHECKS:',count,'quotient submatrices with <=2 core nodes PSD; three-core block non-PSD')
print('KERNELS: exact finite rank, pairwise independence, target strict positivity, complete-Pick representation checked')
print('SCOPE: universal interpolation is a separate analytic argument, not inferred from finite checks')
print('ALL_CHECKS_OK')

#!/usr/bin/env python3
"""Exact source-fidelity and proposed GLT nonclosure certificates. Requires SymPy."""
import sympy as S
from sympy import Matrix, Rational as Q

def equal(A, B):
    assert A.shape == B.shape
    assert all(S.simplify(x) == 0 for x in A-B), (A, B)

def pd2(M):
    equal(M, M.T)
    assert M[0, 0] > 0 and S.simplify(M.det()) > 0

def mean2(A, B):
    """Native positive 2x2 mean (ALM Corollary 2.2), all radicals exact."""
    a, b = S.sqrt(A.det()), S.sqrt(B.det())
    N = A/a + B/b
    return S.sqrt(a*b) * N / S.sqrt(S.simplify(N.det()))

# Fidelity: the exact overlap matrix in source Section 4.4.2,
# 'Candidate symbols predicted by Conjecture 1', Case 1, Example 2.
A0 = Matrix([[2, 1], [1, 2]])
B0 = Matrix([[3, 1], [1, 1]])
N0 = Matrix([[2*S.sqrt(2)+3*S.sqrt(3), S.sqrt(2)+S.sqrt(3)],
             [S.sqrt(2)+S.sqrt(3), 2*S.sqrt(2)+S.sqrt(3)]])
C0 = N0 / (6**Q(1,4)*S.sqrt(2+S.sqrt(6)))
pd2(A0); pd2(B0); pd2(C0)
equal(C0*A0.inv()*C0, B0)
equal(mean2(A0, B0), C0)
print('FIDELITY_OK: source overlap matrix C is exactly A0 # B0')

# One symbolic certificate covers all t>0, s>=0, rather than only small cases.
t = S.symbols('t', positive=True)
s = S.symbols('s', nonnegative=True)
P = S.diag(1, 0)
A = S.diag(1, t*t)
B = Matrix([[1,s*t],[s*t,(1+s*s)*t*t]])
H = (A+B)/S.sqrt(4+s*s)
pd2(A); pd2(B); pd2(H)
equal(H*A.inv()*H, B)
assert S.simplify(A.det()-t*t) == 0
assert S.simplify(B.det()-t*t) == 0
assert S.simplify(H.det()-t*t) == 0
x,y = S.symbols('x y', real=True)
assert S.expand((Matrix([[x,y]])*B*Matrix([x,y]))[0] -
                ((x+s*t*y)**2+t*t*y*y)) == 0
equal(A.applyfunc(lambda z:S.limit(z,t,0,dir='+')), P)
equal(B.applyfunc(lambda z:S.limit(z,t,0,dir='+')), P)
equal(H.applyfunc(lambda z:S.limit(z,t,0,dir='+')), 2*P/S.sqrt(4+s*s))
print('FAMILY_OK: positivity, determinants, Riccati identity, all boundary limits')

C = S.diag(1,4)
Z = Matrix([[1,1],[1,2]])
D = Z*C.inv()*Z
pd2(C);pd2(D);pd2(Z)
equal(D,Matrix([[Q(5,4),Q(3,2)],[Q(3,2),2]]))
equal(mean2(C,D),Z)
equal(C*D-D*C,Matrix([[0,-Q(9,2)],[Q(9,2),0]]))
kappa,xi = S.diag(P,C),S.diag(P,D)
assert kappa.rank() == xi.rank() == 3
assert kappa*xi != xi*kappa
print('AUGMENTATION_OK: fixed noncommuting rank-three symbols; C # D = Z')

n=S.symbols('n',integer=True,positive=True)
limits=[]
for shear,parity in [(Q(3,2),'even'),(Q(8,3),'odd')]:
    An=A.subs(t,1/n)
    Bn=B.subs({t:1/n,s:shear})
    Hn=S.simplify(H.subs({t:1/n,s:shear}))
    if parity=='even':
        target=Matrix([[Q(4,5),Q(3,5)/n],[Q(3,5)/n,Q(17,10)/n**2]])
        tr=Q(19,20)+Q(17,40)/n**2
    else:
        target=Matrix([[Q(3,5),Q(4,5)/n],[Q(4,5)/n,Q(41,15)/n**2]])
        tr=Q(9,10)+Q(41,60)/n**2
    equal(Hn,target)
    assert S.simplify((S.trace(Hn)+S.trace(Z))/4-tr)==0
    lim=S.limit(tr,n,S.oo);limits.append(lim)
    # Exact uniform spectral bounds: PD plus trace bounds all eigenvalues.
    assert S.trace(Hn).subs(n,1) <= Q(10,3)
    assert S.simplify(S.diff(S.trace(Hn),n)) <= 0
    # Frobenius bounds imply the required input operator-norm convergence.
    EA,EB=An-P,Bn-P
    a_sq=S.simplify(S.trace(EA.T*EA));b_sq=S.simplify(S.trace(EB.T*EB))
    assert a_sq == 1/n**4
    assert S.simplify(b_sq-(2*shear**2/n**2+(1+shear**2)**2/n**4))==0
    assert S.limit(a_sq,n,S.oo)==S.limit(b_sq,n,S.oo)==0
    assert S.trace(An).subs(n,1) <= 4
    assert S.trace(Bn).subs(n,1) < 10
    assert S.trace(C) <= 5 and S.trace(D) < 10
    print(f'{parity.upper()}_EXACT: normalized_trace={tr}; subsequential_limit={lim}; perturbation_Frob_squared={b_sq}')
    # Independent direct-sum and tensor Riccati identities, small actual sizes.
    for k in range(1,9):
        a,b,h = [M.subs(n,k) for M in (An,Bn,Hn)]
        AA=S.kronecker_product(S.eye(k),S.diag(a,C))
        BB=S.kronecker_product(S.eye(k),S.diag(b,D))
        HH=S.kronecker_product(S.eye(k),S.diag(h,Z))
        equal(HH*AA.inv()*HH,BB)
        assert HH.trace()/(4*k)==tr.subs(n,k)
assert limits==[Q(19,20),Q(9,10)] and limits[0]-limits[1]==Q(1,20)
assert S.trace(Z)==3 and S.trace(Z)+S.trace(P)==4
e=S.symbols('epsilon',positive=True)
Pe=P+e*S.eye(2)
equal(S.simplify(mean2(Pe,Pe)),Pe)
CDregularized=mean2(C+e*S.eye(2),D+e*S.eye(2))
equal(CDregularized.applyfunc(lambda z:S.simplify(S.limit(z,e,0,dir='+'))),Z)
# An explicit continuous compactly supported complex spectral test:
# F(z)=z*phi(|z|), phi(u)=1 for u<=4, 5-u for 4<=u<=5, 0 for u>=5.
# On every certified output spectrum F is identity; its averages above do not converge.
print('SPECTRAL_TEST_OK: output spectra in [0,10/3] subset [0,4]; compact cutoff equals identity')
print('PROPOSED_SYMBOL_TRACE: tr(P direct_sum Z)/4 = 1')
print('NONCONVERGENCE_CERTIFICATE: even limit 19/20; odd limit 9/10; gap 1/20')
print('SCOPE: exact algebraic and limit certificates; GLT generator/distribution theorems remain Lean obligations')
print('ALL_CHECKS_OK')

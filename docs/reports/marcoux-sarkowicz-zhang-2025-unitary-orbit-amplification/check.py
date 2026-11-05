#!/usr/bin/env python3
"""Exact source-fidelity and conditional padding checks; no concrete gap seed.

ALL_CHECKS_OK refers only to assertions executed below. The published existence
seed is not instantiated, so this program does not certify a counterexample.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement, permutations
from math import gcd
import sympy as s


def bottleneck(xs, ys):
    """Independent threshold perfect-matching calculation, retaining multiplicity."""
    assert len(xs) == len(ys) and xs
    for threshold in sorted({abs(x-y) for x in xs for y in ys}):
        owner = {}
        def augment(i, seen):
            for j, y in enumerate(ys):
                if j in seen or abs(xs[i]-y) > threshold:
                    continue
                seen.add(j)
                if j not in owner or augment(owner[j], seen):
                    owner[j] = i
                    return True
            return False
        if all(augment(i, set()) for i in range(len(xs))):
            return threshold
    raise AssertionError('perfect matching absent')


def fidelity():
    # Section 3.1: ordered spectral matching, and its repetition corollary.
    lists = list(combinations_with_replacement([F(-2), F(0), F(3)], 3))
    pairs = [(x, y) for x in lists for y in lists]
    pairs.append(((F(-3,2), F(1,3), F(7,4)), (F(-1), F(2,3), F(5,2))))
    checks = 0
    for x, y in pairs:
        direct = min(max(abs(a-b) for a,b in zip(x,p)) for p in permutations(y))
        assert direct == bottleneck(x,y) == max(abs(a-b) for a,b in zip(x,y))
        for m in range(1,5):
            assert bottleneck(x*m,y*m) == direct
            checks += 1
    print(f'FIDELITY: {checks} exact spectral-repetition checks, including repeated eigenvalues')
    print('FIDELITY_SCOPE: Section 3.1 matching proposition/corollary; orbit interpretation uses the cited self-adjoint theorem')


def pad(A, lam, r):
    return s.diag(A, lam*s.eye(r))


def regroup(m, n, r):
    # Original indices (copy, base-or-pad) -> all base copies followed by pad copies.
    order = [a*(n+r)+i for a in range(m) for i in range(n)]
    order += [a*(n+r)+n+j for a in range(m) for j in range(r)]
    return s.eye(m*(n+r))[order,:]


def amplification_blocks():
    A = s.Matrix([[0,1],[-1,1]])
    B = s.Matrix([[1,2],[0,-1]])
    lam = s.Rational(19,2)
    checks=0
    for m in range(1,5):
        for r in range(1,5):
            P=regroup(m,2,r)
            assert P*P.T == s.eye(m*(2+r))
            for M in [A,B]:
                assert P*s.kronecker_product(s.eye(m),pad(M,lam,r))*P.T == pad(s.kronecker_product(s.eye(m),M),lam,m*r)
                checks+=1
    print(f'AMPLIFICATION: {checks} exact padded regroupings, m=1..4 and r=1..4')


def near_unitary():
    # A square compression with rational cosine/sine singular values.
    for t in [s.Rational(1,50),s.Rational(1,20),s.Rational(1,9)]:
        c=(1-t*t)/(1+t*t)
        z=2*t/(1+t*t)
        U=s.Matrix([[c,0,-z,0],[0,c,0,-z],[z,0,c,0],[0,z,0,c]])
        assert U.T*U == s.eye(4)
        X=U[:2,:2];Z=U[2:,:2]
        V=s.eye(2)
        assert X.T*X == s.eye(2)-Z.T*Z
        assert s.simplify(s.sqrt(1-z*z)-c)==0
        assert 0 <= 1-c <= z*z < 1
        assert X == V*(c*s.eye(2))
        # Exact block residuals, with genuinely nonnormal 2-by-2 inputs.
        S=s.Matrix([[0,1],[0,0]]);T=s.Matrix([[0,-1],[0,0]])
        lam=s.Integer(9)
        R=pad(S,lam,2)*U-U*pad(T,lam,2)
        Y=U[:2,2:];W=U[2:,2:]
        assert R[:2,:2] == S*X-X*T
        assert R[:2,2:] == (S-lam*s.eye(2))*Y
        assert R[2:,:2] == Z*(lam*s.eye(2)-T)
        assert R[2:,2:] == s.zeros(2)
        # Neumann resolvent bound for both exact nilpotents (M=1).
        for A in [S,T]:
            K=(lam*s.eye(2)-A).inv()
            Q=s.eye(2)/(lam-1)**2-K.T*K
            assert Q[0,0]>=0 and Q[1,1]>=0 and Q.det()>=0
    print('POLAR_BLOCKS: 3 exact unitary compressions, polar estimates, residual identities, and Neumann bounds')


def nonnegative_for_q_ge_one(expr,q):
    # Certify a rational function by nonnegative numerator/denominator coefficients
    # after q=1+t, t>=0. A strictly positive constant certifies strict positivity.
    t=s.Symbol('t',nonnegative=True)
    num,den=s.fraction(s.cancel(expr.subs(q,t+1)))
    pn,pd=s.Poly(num,t),s.Poly(den,t)
    assert all(c>=0 for c in pn.all_coeffs()) and pn.eval(0)>0
    assert all(c>=0 for c in pd.all_coeffs()) and pd.eval(0)>0


def conditional_family_budget():
    q=s.Symbol('q',positive=True)
    # g is a rational lower bound for an actual normalized seed gap, if supplied.
    # These parameters do not assert existence of matrices realizing g.
    g=1/q;lam=1+8*q;error=8/(lam-1)**2;perturb=g/32
    nonnegative_for_q_ge_one(lam-3,q)
    nonnegative_for_q_ge_one(g/2-error,q)
    retained=s.factor(g-error-4*perturb)
    nonnegative_for_q_ge_one(retained-g/4,q)
    gamma=g/(4*(lam+1))
    nonnegative_for_q_ge_one(gamma,q)
    assert s.simplify(error-1/(8*q*q))==0
    assert s.simplify(retained-(7*q-1)/(8*q*q))==0
    for N in [2,4,8,16,32]:
        ks=[k for k in range(N,N+60) if k%3==0]
        assert len(ks)==20
        for k in ks:
            for n in [3,6,9,12]:
                d=gcd(n,k)
                assert d>=3 and d%3==0
    print('CONDITIONAL_FAMILY: all q>=1, lambda=1+8q, error=1/(8q^2), perturb=1/(32q), retained=(7q-1)/(8q^2)>g/4')
    print('GCD: dimensions N=2,4,8,16,32 with the next 20 positive multiples of 3 each')


def irreducible_examples():
    # Simple H, every off-diagonal of Hermitian K nonzero: the joint commutant
    # must be scalar. Compute it independently as a linear-system nullspace.
    for k in range(2,9):
        H=s.diag(*range(k));K=s.ones(k)-s.eye(k)
        columns=[]
        for a in range(k):
            for b in range(k):
                E=s.zeros(k);E[a,b]=1
                columns.append(s.Matrix(list(H*E-E*H)+list(K*E-E*K)))
        system=s.Matrix.hstack(*columns)
        assert len(system.nullspace())==1
    print('IRREDUCIBILITY: exact scalar joint commutants in dimensions 2..8')


if __name__=='__main__':
    fidelity()
    amplification_blocks()
    near_unitary()
    conditional_family_budget()
    irreducible_examples()
    print('SEED_NOT_COMPUTED: no concrete N,S,T,delta_1,delta_3 or positive gap is certified')
    print('CHECK_SCOPE: source fidelity and conditional family mechanism only; no settlement certificate')
    print('ALL_CHECKS_OK')

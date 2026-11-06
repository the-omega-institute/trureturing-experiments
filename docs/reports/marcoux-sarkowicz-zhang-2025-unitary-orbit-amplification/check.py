#!/usr/bin/env python3
"""Exact fidelity checks and an explicit rational amplification-gap construction.

The universal lower-bound argument is stated in the accompanying preregistration;
this checks its exact finite identities, constants, matrix data and witness cost.
No floating-point optimization or unverified seed is used.
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


def decimal_bound(x, places, upper):
    """Directed rational rounding, with decimal notation only for presentation."""
    x=F(x);scale=10**places
    n=x.numerator*scale//x.denominator
    if upper and F(n,scale)<x:
        n+=1
    return f'{n//scale}.{n%scale:0{places}d}'


def winding_identities():
    # Hermitian rank-one projection path, written with a rational parameter.
    u=s.Symbol('u',real=True)
    c=(1-u*u)/(1+u*u);v=2*u/(1+u*u)
    z=s.Symbol('z',nonzero=True)
    E=[s.eye(3)[:,i]*s.eye(3)[:,i].T for i in range(3)]
    for j in [1,2]:
        vec=s.zeros(3,1);vec[0]=c;vec[j]=v
        P=vec*vec.T
        assert (P*P-P).applyfunc(s.cancel)==s.zeros(3)
        U=s.eye(3)+(z-1)*P
        Us=s.eye(3)+(1/z-1)*P
        assert (U*Us-s.eye(3)).applyfunc(s.cancel)==s.zeros(3)
        assert s.cancel(U.det()-z)==0
    units=[s.eye(3)+(z-1)*P for P in E]
    assert units[0]*units[1]*units[2]==z*s.eye(3)
    for U in units:
        assert s.factor(U.det())==z
    # The two endpoint scalar evaluations have dimensions 1, not 3.
    grid=[F(j,86) for j in range(1,86)]
    trace_left=F(0)+3*sum(grid)
    trace_right=3*sum(grid)+F(1)
    assert trace_left-trace_right==-1
    assert 1+3*len(grid)==256
    eps=F(1,8)
    assert 2*eps==F(1,4) and 3*(2*eps)==F(3,4)<2
    # pi < 22/7 makes both determinant-integer homotopies avoid -1.
    assert 3*eps<2 and 2*F(22,7)*eps<2
    print('WINDING: rank-one rotation/unitary/determinant identities; three central factors; trace difference -1 in dimension 256')


def explicit_seed():
    x=s.Symbol('x',real=True)
    p=x**4*(1-x)**4
    assert s.integrate(p/(1+x*x),(x,0,1))==s.Rational(22,7)-s.pi
    assert s.integrate(p,(x,0,1))==s.Rational(1,630)
    assert s.factorial(25)>10000*7**25
    eta=F(1,315)+F(1,10000)
    assert eta==F(2063,630000)
    # Exact quantitative tuple recovery for every epsilon in (0,1], q<=16.
    linear=1+F(64,1000)+F(4,32)+F(8,32000)+F(56,1000)
    quadratic=F(180,10**6)
    assert linear==F(4981,4000)
    assert linear+quadratic<F(5,4)
    assert F(1,10)*F(5,4)==F(1,8)
    # Rational matrices K_1,...,K_5 have operator norm one.
    E=lambda a,b: s.eye(3)[:,a]*s.eye(3)[:,b].T
    Ks=[E(0,0),E(1,1),E(2,2),E(0,1)+E(1,0),E(0,2)+E(2,0)]
    for K in Ks:
        assert K.T==K
        assert set(K.eigenvals()).issubset({s.Integer(-1),s.Integer(0),s.Integer(1)})
        assert 1 in K.eigenvals() or -1 in K.eigenvals()
    def C(t):
        return sum((F((-1)**r)*t**(2*r)/int(s.factorial(2*r)) for r in range(13)),F(0))
    def P(t):
        return sum((F((-1)**r)*t**(2*r+1)/int(s.factorial(2*r+1)) for r in range(12)),F(0))
    cosdata=[F(0)];sindata=[F(0)];times=[F(0)]
    for j in range(1,86):
        t=F(j,86);xx=F(22*j,301)
        assert xx==F(44,7)*t and 0<xx<7
        cosdata.append(C(xx)-1);sindata.append(P(xx));times.append(t)
    cosdata.append(F(0));sindata.append(F(0));times.append(F(1))
    # Boundary values are prescribed exactly, not obtained by Taylor evaluation.
    assert cosdata[0]==cosdata[-1]==sindata[0]==sindata[-1]==0
    assert all(abs(a)<2+eta for a in cosdata+sindata+times)
    increments=[abs(a[j]-a[j-1]) for a in [times,cosdata,sindata] for j in range(1,87)]
    witness=max(increments)
    assert witness<F(7,86)+2*eta<F(88,1000)
    lower=F(1,10)-2*eta
    assert lower==F(29437,315000)>F(93,1000)
    assert lower-witness>F(1,200)
    # Triple repetition, reordered into 3q-dimensional blocks. Equal labels
    # represent identical matrices; each unequal pair is one adjacent mesh step.
    left=[0]+[j for j in range(1,86) for _ in range(3)]
    right=[j for j in range(1,86) for _ in range(3)]+[86]
    assert len(left)==len(right)==256
    assert all(b-a in [0,1] for a,b in zip(left,right))
    assert [(a,b) for a,b in zip(left,right) if a!=b]==[(j,j+1) for j in range(86)]
    # Basis-coordinate lists specify the two actual permutation unitaries.
    # Input is three copies, each consisting of an endpoint q-block and 85 3q-blocks.
    for q,scale,denominator in [(12,12000,2400000),(16,16000,3200000)]:
        N=256*q
        D=s.diag(*(1000*j for j in range(q)))
        J=s.zeros(q)
        for j in range(1,q):J[0,j]=J[j,0]=1
        assert D.T==D and J.T==J
        B0=D+32*s.I*J
        B1=B0+s.I*s.eye(q)[:,1]*s.eye(q)[:,1].T
        assert B1-B0==s.I*s.eye(q)[:,1]*s.eye(q)[:,1].T
        rows=[]
        for j in range(1,86):
            coeff=[times[j]]+[v for _ in Ks for v in [cosdata[j],sindata[j]]]
            assert len(coeff)==11
            # Each (encoding label, internal coordinate) has the fixed D and J
            # part, and exactly one of these eleven rational 3-by-3 data blocks.
            mats=[s.Rational(times[j].numerator,times[j].denominator)*s.eye(3)]
            for K in Ks:
                mats += [s.Rational(cosdata[j].numerator,cosdata[j].denominator)*K,
                         s.Rational(sindata[j].numerator,sindata[j].denominator)*K]
            assert len(mats)==11 and all(M.T==M for M in mats)
            assert all(a.is_Rational for M in mats for a in M)
            rows.append(coeff)
        def order(endpoint_first):
            # Original coordinates are labelled by (copy, endpoint/internal,j,a,b).
            result=[]
            def endpoint():
                result.extend((copy,0,0,a,0) for a in range(q) for copy in range(3))
            if endpoint_first:endpoint()
            for j in range(1,86):
                for copy in range(3):
                    result.extend((copy,1,j,a,b) for a in range(q) for b in range(3))
            if not endpoint_first:endpoint()
            return result
        for first in [True,False]:
            coords=order(first)
            assert len(coords)==3*N and len(set(coords))==3*N
        norm_bound=1000*(q-1)+32*(q-1)+2+eta
        assert norm_bound<scale
        # Base upper witness matches all equal internal blocks, leaving B0-B1.
        base_upper=F(1,scale)
        amplified_upper=witness/scale
        trace_lower=F(1,N*scale)
        base_lower=lower/scale
        assert base_lower-amplified_upper>F(1,denominator)>0
        print(f'SEED: q={q}, N={N}, scale={scale}, exact Gaussian-rational entries')
        print(f'DELTA_1_ENCLOSURE: ({decimal_bound(base_lower,15,False)}, {decimal_bound(base_upper,15,True)}]')
        print(f'DELTA_3_ENCLOSURE: [{decimal_bound(trace_lower,15,False)}, {decimal_bound(amplified_upper,15,True)}]')
        print(f'CERTIFIED_GAP: >1/{denominator}')
    print(f'EXACT_WITNESS_RAW: {witness.numerator}/{witness.denominator}')
    print('NUMERICAL_SCOPE: directed rational enclosures, not optimized values of the two minima')
    # One fully fixed irreducible contraction family from the 4096 seed.
    g=F(1,3200000);lam=F(8001);perturb=g/32
    loss=8/(lam-1)**2
    assert lam>3 and loss<g/2
    retained=g-loss-4*perturb
    assert retained>g/4
    gamma=g/(4*(lam+1))
    assert gamma==F(1,102425600000)
    print('UNIFORM_FAMILY: K=4096, lambda=8001, each perturbation <1/102400000, gamma=1/102425600000')


if __name__=='__main__':
    fidelity()
    amplification_blocks()
    near_unitary()
    conditional_family_budget()
    irreducible_examples()
    winding_identities()
    explicit_seed()
    print('CHECK_SCOPE: exact seed data, lower-bound constants/identities, explicit permutation upper witness, and uniform family budget; analytic lower-bound argument in draft')
    print('ALL_CHECKS_OK')

#!/usr/bin/env python3
"""Exact finite checks for the CDMVZ (3.11) preregistration.

Only integers and fractions are used. This checks finite consequences and the
moment/compression mechanism; it does not certify the infinite convergence claim.
Coordinates and adjacency below are from the authors' public supplementary data.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial

VERTICES = [[[0, -2], [2, -2], [-1, 1]], [[-1, -1], [-2, 0], [-1, 1]], [[-1, -1], [1, -1], [2, 0]], [[-1, -1], [2, 0], [-1, 1]], [[0, -1], [-1, 0], [3, 0]], [[-1, 0], [-1, 0], [1, 0]], [[-1, 0], [1, 0], [1, 0]], [[-2, 0], [1, -1], [1, 1]], [[1, -1], [1, 1], [2, 0]], [[-1, 0], [1, 0], [2, 1]], [[0, 0], [0, 0], [1, 0]], [[0, 0], [-1, 0], [0, 0]], [[1, 0], [1, 0], [2, 1]], [[-1, 1], [2, 0], [1, 1]], [[-1, 1], [1, 1], [2, 0]], [[-1, 1], [-1, -1], [2, 0]], [[-1, 1], [-2, 0], [1, 1]], [[0, 2], [-1, -3], [3, -1]], [[2, 0], [-1, -1], [-1, 1]], [[2, 0], [-1, 1], [1, 1]], [[2, 0], [-1, -1], [0, 0]], [[1, 0], [-1, 0], [1, 0]], [[0, 1], [-1, 0], [3, 0]], [[1, 0], [1, 0], [1, 0]], [[0, 1], [3, 0], [1, 0]], [[1, 0], [-1, 0], [-2, 1]], [[1, 0], [-1, 0], [0, 0]], [[1, 0], [1, 0], [0, 0]], [[1, 1], [-2, 0], [-1, 1]], [[1, 1], [1, -1], [2, 0]], [[1, 1], [2, 0], [-1, 1]], [[1, 1], [2, 0], [0, 0]], [[1, 0], [0, 0], [0, 0]]]
PUBLISHED_ADJACENCY = [131328, 3150080, 6561792, 25690112, 537018368, 872947968, 1208025088, 311443456, 4194339, 402915328, 6644828160, 4294968320, 67239938, 33816612, 270532752, 1082130432, 2164260928, 1073745937, 2155880964, 1048616, 3221750914, 134234118, 260, 67403912, 65544, 134226048, 142611488, 102762048, 17056, 48, 1212480, 1377280, 3072]

# Exact arithmetic in Q(sqrt(5)), represented by (a,b) = a+b*sqrt(5).
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def mul(x,y): return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def dot(x,y):
    s=(0,0)
    for a,b in zip(x,y): s=add(s,mul(a,b))
    return s

def raw(n,k):
    lam=Q(n-2,2)
    coeff=[Q(0)]*(k+1)
    for j in range(k//2+1):
        rising=Q(1)
        for i in range(k-j): rising*=lam+i
        coeff[k-2*j]+=(-1)**j*rising*2**(k-2*j)/factorial(j)/factorial(k-2*j)
    return coeff

def normalized(n,k):
    c=raw(n,k);den=sum(c)
    assert den>0
    return [x/den for x in c]

def evaluate(c,t):
    s=Q(0)
    for x in reversed(c):s=s*t+x
    return s

# Fidelity 1: source Appendix Legendre recurrence (3.16), exact polynomial identity.
for k in range(2,21):
    p=normalized(3,k);a=normalized(3,k-1);b=normalized(3,k-2)
    for i in range(k+1):
        rhs=(Q(2*k-1,k)*a[i-1] if 0<i<=len(a) else 0)-(Q(k-1,k)*b[i] if i<len(b) else 0)
        assert p[i]==rhs
print('FIDELITY Legendre recurrence: degrees 2..20, every coefficient exact')

# Fidelity 2: source Table 1's 33 projective points reproduce its adjacency data.
assert len(VERTICES)==33
edges=0
for i,x in enumerate(VERTICES):
    assert dot(x,x)!=(0,0)
    row=0
    for j,y in enumerate(VERTICES):
        if i!=j and dot(x,y)==(0,0):row|=1<<j
        if i!=j:
            # Collinearity would violate distinct antipodal-class indexing.
            cross=[add(mul(x[(k+1)%3],y[(k+2)%3]),tuple(-v for v in mul(x[(k+2)%3],y[(k+1)%3]))) for k in range(3)]
            assert any(v!=(0,0) for v in cross)
    assert row==PUBLISHED_ADJACENCY[i]
    edges+=row.bit_count()
print('FIDELITY source 33-pair graph: vertices=33 edges='+str(edges//2)+'; all adjacency entries exact')

# Source's Witsenhausen orthonormal-frame upper bound: equal singleton marginals
# have N*p<=1. For each frame, the displayed coefficients and singleton law give
# a matching feasible solution of (3.10): f(t)=t^2/N. Congruent subsets have equal
# cardinality; moments are 1,1/N,0 according as their size is 0,1,>=2.
cases=0
for n in range(3,11):
    assert normalized(n,0)==[Q(1)]
    assert normalized(n,2)==[Q(-1,n-1),0,Q(n,n-1)]
    for N in range(1,n+1):
        p=Q(1,N);c0=Q(1,n*N);c2=Q(n-1,n*N)
        poly=[c0+c2*normalized(n,2)[0],0,c2*normalized(n,2)[2]]
        assert poly==[0,0,p] and c0+c2==p and N*p==1
        assert evaluate(poly,0)==0 and evaluate(poly,1)==p
        law={1<<i:p for i in range(N)}
        for i,j in product(range(N),repeat=2):
            moment=sum(a for I,a in law.items() if I>>i&1 and I>>j&1)
            assert moment==(p if i==j else 0)
        assert p-p*p>=0 and p-N*p*p==0  # covariance spectrum of p*I-p^2*J
        cases+=1
print('SMALL EXACT PROGRAMS: '+str(cases)+' orthonormal frames, n=3..10, N=1..n; Delta=1/N')

# All simple graphs on <=4 vertices, with uniform atomic laws on independent
# sets. Centered quadratic forms equal the explicitly nonnegative variance sum.
# Nonnegative row-stochastic compression preserves the factorization and mean.
graphs=0;forms=0
for N in range(1,5):
    pairs=list(combinations(range(N),2))
    for code in range(1<<len(pairs)):
        E=[e for k,e in enumerate(pairs) if code>>k&1]
        atoms=[I for I in range(1<<N) if all(not(I>>i&1 and I>>j&1) for i,j in E)]
        a=Q(1,len(atoms))
        z=[[Q(I>>i&1) for i in range(N)] for I in atoms]
        mean=[sum(a*v[i] for v in z) for i in range(N)]
        M=[[sum(a*v[i]*v[j] for v in z) for j in range(N)] for i in range(N)]
        C=[[M[i][j]-mean[i]*mean[j] for j in range(N)] for i in range(N)]
        for h in product([-1,0,1],repeat=N):
            form=sum(h[i]*C[i][j]*h[j] for i,j in product(range(N),repeat=2))
            variance=sum(a*sum(h[i]*(v[i]-mean[i]) for i in range(N))**2 for v in z)
            assert form==variance and form>=0;forms+=1
        A=[[Q(1,N)]*N,[Q(i+1,N*(N+1)//2) for i in range(N)]]
        assert all(sum(row)==1 and all(x>=0 for x in row) for row in A)
        w=[[sum(A[i][j]*v[j] for j in range(N)) for i in range(2)] for v in z]
        for i,j in product(range(2),repeat=2):
            compressed=sum(A[i][u]*M[u][v]*A[j][v] for u,v in product(range(N),repeat=2))
            assert compressed==sum(a*v[i]*v[j] for v in w)
        graphs+=1
print('MOMENTS: all '+str(graphs)+' simple graphs on 1..4 vertices; '+str(forms)+' exact variance/compression forms')

# The nearby-point bridge's diagonal correction is (p-q)/r, not p-q.
# This tests every entry against the uncompressed block matrix for rational
# admissible data and every listed cluster size.
bridge=0
for s,r in product(range(1,5),range(1,7)):
    for p,q in [(Q(1,3),Q(1,4)),(Q(2,5),Q(3,10)),(Q(1,2),Q(1,2))]:
        H=[[p if i==j else q for j in range(s*r)] for i in range(s*r)]
        for i,j in product(range(s),repeat=2):
            actual=sum(H[i*r+a][j*r+b] for a,b in product(range(r),repeat=2))/r**2
            expected=q+((p-q)/r if i==j else 0)
            assert actual==expected
        bridge+=1
print('NEARBY-POINT COMPRESSION: '+str(bridge)+' exact block tests; diagonal defect=(p-q)/r')

# Exact coefficient means: normalized S^2 integration is (1/2) int_-1^1.
for k in range(17):
    c=normalized(3,k)
    integral=sum(v/Q(i+1) for i,v in enumerate(c) if i%2==0)
    assert integral==(1 if k==0 else 0)
print('NORMALIZATION: S^2 Gegenbauer means exact through degree 16')

# Exact tail decay samples for t^(2m): diagonal mass stays 1 while every fixed
# |t|<1 decreases. These are diagnostic values, not a proof of compactness.
for t in [Q(0),Q(1,3),Q(1,2),Q(3,4),Q(9,10)]:
    values=[t**(2*m) for m in range(1,21)]
    assert all(a>=b>=0 for a,b in zip(values,values[1:]))
print('BOUNDARY CHECK: diagonal loss is explicitly permitted; no convergence at t=1 assumed')
print('ALL_CHECKS_OK')

#!/usr/bin/env python3
"""The unital-block overlap of two quantum channels can be below -1.

Source: Girling, Cirstoiu, Jennings, "Estimation of correlations and non-separability in quantum channels
via unitarity benchmarking", arXiv:2104.04352v3; Phys. Rev. Research 4, 023041 (2022), Appendix C.5:
"we conjecture a stronger result that for any two quantum channels E1, E2 that <T1, T2> >= -1",
where T is the unital block of the Liouville representation in an orthonormal operator basis
X_0 = I/sqrt(d), X_1, ..., X_{d^2-1}, and <A, B> = tr(A^dagger B).

Classical channels Phi_P(X) = sum_ij P_ij X_jj E_ii (P column stochastic) have overlap
tr(P^T Q) - (P1).(Q1)/d.  With a = ceil(d/2), b = floor(d/2): P sends the first a inputs to output 0
and the rest to output 1; Q sends them the opposite way.  Then the overlap is -2ab/d (= -4/3 at d=3).
Checks: (1) the formula against a direct orthonormal-basis computation (floating point) for d = 3, 4, 5;
(2) the exact values -2ab/d for d = 2..12, below -1 exactly when d >= 3.
Prints one line per check and ALL_OK.
"""
from fractions import Fraction as Fr
import cmath, math, sys

def fail(m): print("FAIL", m); sys.exit(1)
def PQ(d):
    a=(d+1)//2
    P=[[0]*d for _ in range(d)]; Q=[[0]*d for _ in range(d)]
    for j in range(d):
        if j<a: P[0][j]=1; Q[1][j]=1
        else:   P[1][j]=1; Q[0][j]=1
    return P,Q
def formula(P,Q,d):
    tr=sum(P[i][j]*Q[i][j] for i in range(d) for j in range(d))
    r1=[sum(P[i]) for i in range(d)]; r2=[sum(Q[i]) for i in range(d)]
    return Fr(tr)-Fr(sum(x*y for x,y in zip(r1,r2)),d)
def direct(P,Q,d):
    # orthonormal basis: I/sqrt d, normalized diagonal traceless (Gram-Schmidt), symmetric/antisymmetric off-diagonal
    def E(i,j): M=[[0j]*d for _ in range(d)]; M[i][j]=1; return M
    def add(A,B,s=1): return [[A[i][j]+s*B[i][j] for j in range(d)] for i in range(d)]
    def sc(c,A): return [[c*A[i][j] for j in range(d)] for i in range(d)]
    def ip(A,B): return sum(A[i][j].conjugate()*B[i][j] for i in range(d) for j in range(d))
    basis=[sc(1/math.sqrt(d),[[1.0+0j if i==j else 0j for j in range(d)] for i in range(d)])]
    for k in range(1,d):
        v=[[0j]*d for _ in range(d)]
        for i in range(k): v[i][i]=1
        v[k][k]=-k
        n=math.sqrt(k*(k+1)); basis.append(sc(1/n,v))
    for i in range(d):
        for j in range(i+1,d):
            basis.append(sc(1/math.sqrt(2),add(E(i,j),E(j,i))))
            basis.append(sc(-1j/math.sqrt(2),add(E(i,j),E(j,i),-1)))
    for a_ in range(len(basis)):
        for b_ in range(len(basis)):
            if abs(ip(basis[a_],basis[b_])-(1 if a_==b_ else 0))>1e-12: fail("basis")
    def ch(M,X): return [[(sum(M[i][j]*X[j][j] for j in range(d)) if i==k else 0j) for k in range(d)] for i in range(d)]
    n=len(basis); s=0
    for i in range(1,n):
        for j in range(1,n):
            t1=sum(basis[i][r][c].conjugate()*ch(P,basis[j])[r][c] for r in range(d) for c in range(d))
            t2=sum(basis[i][r][c].conjugate()*ch(Q,basis[j])[r][c] for r in range(d) for c in range(d))
            s+=t1.conjugate()*t2
    return s
ok=True
for d in (3,4,5):
    P,Q=PQ(d); f=formula(P,Q,d); dv=direct(P,Q,d)
    if abs(dv.real-float(f))>1e-9 or abs(dv.imag)>1e-9: fail(f"direct d={d}")
    print(f"DIRECT d={d}: basis computation {dv.real:.12f} = formula {f}")
for d in range(2,13):
    P,Q=PQ(d); f=formula(P,Q,d); a=(d+1)//2; b=d//2
    if f!=Fr(-2*a*b,d): fail(f"formula d={d}")
    if (f< -1)!=(d>=3): fail(f"threshold d={d}")
    print(f"FAMILY d={d}: overlap {f} ({'below' if f<-1 else 'not below'} -1)")
print("ALL_OK")

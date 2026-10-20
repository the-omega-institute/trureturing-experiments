"""Chirvasitu-Cusick (arXiv:2304.12734v1) Conjecture 6.1: ord(R(i,j)/sqrt 2) = 2 lcm(4, i-j, i+j), eigenvalues sqrt2-scaled
roots of unity. R(i,j) is built from Eq. (6-3): 2^(i-1) row blocks rho_k pi_{2^(i-1)-1} H_eps(k), eps(k) = floor(k/2^(j-1)) mod 2.
Checks:
  1. R(i,j)^i is Hadamard for all i <= 10 (a property the source states);
  2. exact integer check at (8,4): R^16 != 2^8 I, R^24 != 2^12 I (trace 0), R^48 = 2^24 I, so the order is 48, not 24;
  3. for all 1 <= j < i <= 10, the exact order of R/sqrt2 (floating powers, divisor descent from the candidate) equals
     the corrected formula 2 lcm(4, i-j, i+j, 2 gcd(i,j)); the conjectured formula fails only at (8,4).
"""
import numpy as np
from math import lcm, gcd
from sympy import primefactors
H0=np.array([[1,1],[1,-1]]); H1=np.array([[1,1],[-1,1]])
def R(i,j):
    n=2**i; blocks=[]
    for k in range(2**(i-1)):
        H = H1 if (k // 2**(j-1)) % 2 == 1 else H0
        M=np.zeros((2,n),dtype=np.int64)
        for c in range(2): M[:,c*2**(i-1)]=H[:,c]
        blocks.append(np.roll(M,k,axis=1))
    return np.vstack(blocks)
def mpow_int(M,e):
    n=M.shape[0]; P=np.identity(n,dtype=np.int64); B=M.copy()
    while e:
        if e&1: P=P@B
        B=B@B; e>>=1
    return P
ok=True
for i in range(2,11):
    for j in range(1,i):
        P=mpow_int(R(i,j),i); n=2**i
        had=set(np.unique(P)) <= {1,-1} and (P@P.T == n*np.identity(n,dtype=np.int64)).all()
        ok &= bool(had)
print("R(i,j)^i Hadamard for all 1<=j<i<=10:", ok)
M=R(8,4); I=np.identity(256,dtype=np.int64); res={}
for e in (16,24,48):
    P=mpow_int(M,e); res[e]=bool((P==(2**(e//2))*I).all()); print(f"(8,4) exact: R^{e} == 2^{e//2} I ? {res[e]}  trace={int(np.trace(P))}")
ok &= (not res[16]) and (not res[24]) and res[48]
def is_identity_power(U,e):
    n=U.shape[0]; P=np.identity(n); B=U.copy()
    while e:
        if e&1: P=P@B
        B=B@B; e>>=1
    return np.abs(P-np.identity(n)).max() < 1e-6
def exact_order(U, D):
    if not is_identity_power(U,D): return None
    for p in primefactors(D):
        while D % p == 0 and is_identity_power(U, D//p): D//=p
    return D
bad=[]
for i in range(2,11):
    for j in range(1,i):
        U=R(i,j).astype(float)/np.sqrt(2)
        corr=2*lcm(4,i-j,i+j,2*gcd(i,j)); conj=2*lcm(4,i-j,i+j)
        o=exact_order(U, corr)
        if o!=corr or o!=conj: bad.append((i,j,o,corr,conj))
        print(f"i={i} j={j}: order={o} corrected={corr} conjectured={conj}", flush=True)
mm=[b for b in bad if b[2]!=b[3]]; cf=[(b[0],b[1],b[2],b[4]) for b in bad if b[2]!=b[4]]
print("mismatches with corrected:", mm); print("conjecture failures:", cf)
ok &= (mm==[]) and (cf==[(8,4,48,24)])
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)

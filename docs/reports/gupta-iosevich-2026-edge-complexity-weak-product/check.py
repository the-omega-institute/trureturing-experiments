# Gupta–Iosevich–Iosevich–Song–Tian (arXiv:2607.15598) weak-product closure question: K3 x K3 check.
import itertools, random
import numpy as np
from fractions import Fraction as Fr
N=9
def F(n): return np.array([[np.exp(-2j*np.pi*m*x/n)/np.sqrt(n) for x in range(n)] for m in range(n)])
def FR(A):
    n=A.shape[0]; M=F(n)@A@F(n); return np.abs(M).sum()/np.sqrt((np.abs(M)**2).sum())
def energy(A): return np.abs(np.linalg.eigvalsh(A)).sum()
K3=np.ones((3,3))-np.eye(3)
V=[(i,j) for i in range(3) for j in range(3)]
P=np.array([[1.0 if (u[0]!=v[0] and u[1]!=v[1]) else 0.0 for v in V] for u in V])
s3=3; sP=int(P.sum()/2)
print("K3: FR =",FR(K3),"E/sqrt(2s) =",energy(K3)/np.sqrt(2*s3))
print("K3xK3: E =",energy(P),"s =",sP,"bound =",energy(P)/np.sqrt(2*sP))
J=np.ones((9,9)); I=np.eye(9)
assert np.allclose(P@P,2*J+2*I-P)
Sg=(6*P+3*I-2*J)/9; assert np.allclose(Sg@Sg,I); assert np.isclose(np.trace(Sg@P),16)
best=1e9; random.seed(1)
for t in range(20000):
    p=list(range(9)); random.shuffle(p); Q=np.eye(9)[p]; best=min(best,FR(Q@P@Q.T))
print("min FR over 20000 random labelings =",best)
# exhaustive over labelings up to the 72-element automorphism group is 9!/72 = 5040 orbits; do all 9! / fix vertex 0 -> 8! = 40320
best2=1e9
for p in itertools.permutations(range(1,9)):
    q=(0,)+p; Q=np.eye(9)[list(q)]; best2=min(best2,FR(Q@P@Q.T))
print("min FR over all labelings fixing label 0 (8! = 40320; circulant shift makes this exhaustive) =",best2)
# no circulant 0/1 matrix C on Z9 with C^2 = 2J+2I-C
cnt=0
for a in itertools.product([0,1],repeat=9):
    C=np.array([[a[(y-x)%9] for y in range(9)] for x in range(9)],dtype=float)
    if np.allclose(C@C,2*J+2*I-C): cnt+=1
print("circulant solutions of C^2 = 2J+2I-C on Z9:",cnt)

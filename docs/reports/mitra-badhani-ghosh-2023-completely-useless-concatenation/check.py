#!/usr/bin/env python3
"""Concatenating two completely useless qubit channels can be useful under the quantum switch.

Source: Mitra, Badhani, Ghosh, "Improvement in quantum communication using quantum switch",
arXiv:2108.14001v6; Phys. Scr. 98, 045101 (2023), Section III.B, Conjecture 1:
"Concatenation of two completely useless channels is always completely useless."
A channel is completely useless when it is entanglement breaking and both switch branches
C_{+-}(rho) = (1/4) sum_{x,y} (A_x A_y +- A_y A_x) rho (A_x A_y +- A_y A_x)^dagger
are entanglement-breaking CP maps (two-qubit Choi matrix: separable iff PPT).

Family: A_{eta,u} = eta D_u + (1 - eta) D with D_u(rho) = (rho + s_u rho s_u)/2 and D(rho) = Tr(rho) I/2;
A = A_{eta,z}, B = A_{eta,n}, n = (4/5, 0, 3/5), T = B o A.  Exact Gaussian-rational arithmetic:
the Choi matrices are polynomial in the squared Kraus weights, so they are exact for rational eta.
Checks, for eta in {19/20, 39/40, 99/100}:
  * the partial transposes of the Choi matrices of C_+(A), C_-(A), C_+(B), C_-(B) have all principal minors >= 0 (PSD);
  * for C_-(T) the vector w = (-2,-1,-1,2) gives w^dagger J^Gamma w / 10 = (9 eta^4 - 150 eta^2 + 125)/2000 < 0.
Prints one line per eta and ALL_OK.
"""
from fractions import Fraction as Fr
from itertools import combinations, product
import sys

class C:
    __slots__=("r","i")
    def __init__(s,r,i=Fr(0)): s.r=Fr(r); s.i=Fr(i)
    def __add__(s,o): o=o if isinstance(o,C) else C(o); return C(s.r+o.r,s.i+o.i)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,C) else C(o); return C(s.r-o.r,s.i-o.i)
    def __mul__(s,o): o=o if isinstance(o,C) else C(o); return C(s.r*o.r-s.i*o.i,s.r*o.i+s.i*o.r)
    __rmul__=__mul__
    def conj(s): return C(s.r,-s.i)
    def __eq__(s,o): o=o if isinstance(o,C) else C(o); return s.r==o.r and s.i==o.i
def fail(m): print("FAIL",m); sys.exit(1)
def mat(n,m,f): return [[f(i,j) for j in range(m)] for i in range(n)]
def mm(A,B): return mat(len(A),len(B[0]),lambda i,j: sum((A[i][k]*B[k][j] for k in range(len(B))),C(0)))
def H(A): return mat(len(A[0]),len(A),lambda i,j: A[j][i].conj())
def add(A,B): return mat(len(A),len(A[0]),lambda i,j: A[i][j]+B[i][j])
def scal(c,A): return mat(len(A),len(A[0]),lambda i,j: c*A[i][j])
def kron(A,B): return mat(len(A)*len(B),len(A[0])*len(B[0]),lambda i,j: A[i//len(B)][j//len(B[0])]*B[i%len(B)][j%len(B[0])])
I=[[C(1),C(0)],[C(0),C(1)]]; X=[[C(0),C(1)],[C(1),C(0)]]; Y=[[C(0),C(0,-1)],[C(0,1),C(0)]]; Z=[[C(1),C(0)],[C(0),C(-1)]]
# Kraus families as (weight, operator) with weight = squared coefficient (channel = sum w K rho K^dag)
def fam(eta,n):
    sn=add(add(scal(C(n[0]),X),scal(C(n[1]),Y)),scal(C(n[2]),Z))
    return [(eta/2,I),(eta/2,sn)]+[((1-eta)/4,O) for O in (I,X,Y,Z)]
def compose(Bf,Af): return [(wa*wb, mm(b,a)) for (wa,a) in Af for (wb,b) in Bf]
def branch(F,sgn):
    # C_s(rho) = (1/4) sum_{x,y} sqrt(w_x w_y)^2 (A_x A_y + s A_y A_x) rho (...)^dag : weights multiply
    return [(wx*wy/4, add(mm(ax,ay), scal(C(sgn),mm(ay,ax)))) for (wx,ax) in F for (wy,ay) in F]
def choi(F):
    # J = sum_w (I (x) K) |Om><Om| (I (x) K)^dag with |Om> = (|00>+|11>)/sqrt2 -> use unnormalized omega and divide by 2
    om=[[C(1)],[C(0)],[C(0)],[C(1)]]; P=mm(om,H(om)); J=mat(4,4,lambda i,j:C(0))
    for w,K in F:
        KK=kron(I,K); J=add(J,scal(C(w/2),mm(mm(KK,P),H(KK))))
    return J
def pt_out(J): return mat(4,4,lambda a,b: J[(a//2)*2+(b%2)][(b//2)*2+(a%2)])
def det(M):
    n=len(M)
    if n==1: return M[0][0]
    return sum(((C(-1) if j%2 else C(1))*M[0][j]*det([row[:j]+row[j+1:] for row in M[1:]]) for j in range(n)),C(0))
def psd(M):
    n=len(M)
    for k in range(1,n+1):
        for S in combinations(range(n),k):
            d=det([[M[i][j] for j in S] for i in S])
            if d.i!=0: fail("non-real minor")
            if d.r<0: return False
    return True
ok=True
for eta in (Fr(19,20),Fr(39,40),Fr(99,100)):
    A=fam(eta,(0,0,1)); B=fam(eta,(Fr(4,5),0,Fr(3,5))); T=compose(B,A)
    flags={}
    for name,F in (("A",A),("B",B)):
        for s,sg in (("+",1),("-",-1)):
            flags[name+s]=psd(pt_out(choi(branch(F,sg))))
    J=pt_out(choi(branch(T,-1))); w=[C(-2),C(-1),C(-1),C(2)]
    val=sum((w[i].conj()*J[i][j]*w[j] for i in range(4) for j in range(4)),C(0))
    pred=(9*eta**4-150*eta**2+125)/2000
    good=all(flags.values()) and val.i==0 and val.r/10==pred and pred<0
    ok&=good
    print(f"eta={eta}: PPT(C+-(A), C+-(B))={flags} witness/10={val.r/10} predicted={pred} negative={pred<0}")
print("ALL_OK" if ok else "FAIL"); sys.exit(0 if ok else 1)

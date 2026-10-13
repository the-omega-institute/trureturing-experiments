#!/usr/bin/env python3
"""Four mutually unbiased bases detect bound entanglement in dimension eight (Conjecture 1 refuted).

Source: Bae, Bera, Chruscinski, Hiesmayr, McNulty, "How many mutually unbiased bases are needed to detect bound
entangled states?", arXiv:2108.01109v2; J. Phys. A 55, 505303 (2022), Conjecture 1: for d = 2^r, r > 1, the minimal
number of MUBs required for bound entanglement detection equals d/2 + 1 (= 5 at d = 8).

Witness (source Eq. (4) and its partial transpose): W^Gamma(M_m, s) = ((d+m-1)/d) I - sum_l |l><l| (x) |l+s><l+s|
- sum_{alpha=1}^{m-1} sum_i P_{alpha i} (x) conj(P_{alpha i}), here d = 8, m = 4, s = 1.
MUBs: B0 = computational basis, B1 = H/sqrt8, B2 = diag(a) H/sqrt8, B3 = diag(b) H/sqrt8 (H the 8x8 Walsh matrix).
State: rho = M/10134, block diagonal over the XOR sectors h, data in certificate.json.
Exact Gaussian-rational checks:
  (1) the four bases are orthonormal and mutually unbiased (|<u,v>|^2 = 1/8);
  (2) rho is Hermitian with trace 1 and positive definite (exact LDL pivots > 0);
  (3) rho^Gamma (partial transpose on the second factor) is positive definite (exact LDL pivots > 0);
  (4) Tr(W rho) = -14/563 exactly.
Prints one line per check and ALL_OK.
"""
import json, os, sys
from fractions import Fraction as Fr
def fail(m): print("FAIL", m); sys.exit(1)
here=os.path.dirname(os.path.abspath(__file__)); cert=json.load(open(os.path.join(here,"certificate.json")))
d=8
class G:  # Gaussian rational
    __slots__=("r","i")
    def __init__(s,r,i=0): s.r=Fr(r); s.i=Fr(i)
    def __add__(s,o): o=o if isinstance(o,G) else G(o); return G(s.r+o.r,s.i+o.i)
    __radd__=__add__
    def __sub__(s,o): o=o if isinstance(o,G) else G(o); return G(s.r-o.r,s.i-o.i)
    def __mul__(s,o): o=o if isinstance(o,G) else G(o); return G(s.r*o.r-s.i*o.i,s.r*o.i+s.i*o.r)
    __rmul__=__mul__
    def conj(s): return G(s.r,-s.i)
    def __eq__(s,o): o=o if isinstance(o,G) else G(o); return s.r==o.r and s.i==o.i
def parse(t): return {"1":G(1),"-1":G(-1),"i":G(0,1),"-i":G(0,-1)}[t]
a=[parse(t) for t in cert["a"]]; b=[parse(t) for t in cert["b"]]
def dot(x,k): return bin(x&k).count("1")%2
# basis vectors scaled by sqrt8: w[alpha][k][x]; true vectors are w/sqrt8 (alpha>=1)
def vec(alpha,k):
    if alpha==0: return [G(1) if x==k else G(0) for x in range(d)]
    ph=[G(1)]*d if alpha==1 else (a if alpha==2 else b)
    return [ph[x]*((-1)**dot(x,k)) for x in range(d)]
scale={0:Fr(1),1:Fr(1,8),2:Fr(1,8),3:Fr(1,8)}   # squared normalisation
def ip(u,v): return sum((u[x].conj()*v[x] for x in range(d)),G(0))
for al in range(4):
    for be in range(4):
        for k in range(d):
            for l in range(d):
                g=ip(vec(al,k),vec(be,l)); n2=(g.r*g.r+g.i*g.i)*scale[al]*scale[be]
                if al==be:
                    if n2!=(1 if k==l else 0): fail(f"orthonormal {al}")
                elif n2!=Fr(1,8): fail(f"unbiased {al},{be}")
print("MUBS: four orthonormal, mutually unbiased bases in C^8")
order=[(int(p[0]),int(p[1])) for p in cert["upper_triangle_order"]]
def sym(up):
    M=[[0]*4 for _ in range(4)]
    for (i,j),v in zip(order,up): M[i][j]=M[j][i]=v
    return M
N=64; R=[[Fr(0)]*N for _ in range(N)]
for blk in cert["blocks"]:
    h=blk["h"]; A=sym(blk["A_upper"]); B=sym(blk["B_upper"])
    Mb=[[ (A if (i<4)==(j<4) else B)[i%4][j%4] for j in range(8)] for i in range(8)]
    for xi in range(8):
        for yi in range(8):
            R[xi*d+(xi^h)][yi*d+(yi^h)]=Fr(Mb[xi][yi],cert["denominator"])
if sum(R[i][i] for i in range(N))!=1: fail("trace")
if any(R[i][j]!=R[j][i] for i in range(N) for j in range(N)): fail("symmetric")
def posdef(M):
    M=[row[:] for row in M]; n=len(M)
    for k in range(n):
        if M[k][k]<=0: return False
        for i in range(k+1,n):
            if M[i][k]!=0:
                f=M[i][k]/M[k][k]
                for j in range(k,n): M[i][j]-=f*M[k][j]
    return True
if not posdef(R): fail("rho")
print("RHO: real symmetric, trace 1, positive definite")
RG=[[R[(i//d)*d+(j%d)][(j//d)*d+(i%d)] for j in range(N)] for i in range(N)]
if not posdef(RG): fail("rho^Gamma")
print("RHO_GAMMA: positive definite")
# Tr(W rho) with W = 11/8 I - C0 - sum_alpha sum_k P (x) conj(P)
tr=Fr(11,8)*sum(R[i][i] for i in range(N))
tr-=sum(R[x*d+(x+1)%d][x*d+(x+1)%d] for x in range(d))
tot=G(0)
for al in (1,2,3):
    for k in range(d):
        u=vec(al,k)
        for x in range(d):
            for z in range(d):
                for y in range(d):
                    for w in range(d):
                        r=R[y*d+w][x*d+z]
                        if r==0: continue
                        # (P (x) conj P)[(x,z),(y,w)] = u[x] conj u[y] conj u[z] u[w] / 64
                        tot=tot+(u[x]*u[y].conj()*u[z].conj()*u[w])*(r/64)
if tot.i!=0: fail("imaginary part")
val=tr-tot.r
if val!=Fr(-14,563): fail(f"trace value {val}")
print(f"WITNESS: Tr(W rho) = {val}")
print("ALL_OK")

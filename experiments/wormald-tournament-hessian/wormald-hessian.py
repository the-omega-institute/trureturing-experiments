import sys,random
from fractions import Fraction as Fr
from math import isqrt
def det(M):
  M=[r[:] for r in M]; n=len(M); d=Fr(1)
  for c in range(n):
    p=next((r for r in range(c,n) if M[r][c]!=0),None)
    if p is None: return Fr(0)
    if p!=c: M[c],M[p]=M[p],M[c]; d=-d
    d*=M[c][c]; inv=1/M[c][c]
    for r in range(c+1,n):
      f=M[r][c]*inv
      if f:
        for k in range(c,n): M[r][k]-=f*M[c][k]
  return d
def calH(A):
  d=len(A); z=Fr(2,d)
  C=[[Fr(1,2)*(i==j)+A[i][j] for j in range(d)] for i in range(d)]  # (1/2 I + A)
  B=[[Fr(int(i==j))-C[i][j]*z for j in range(d)] for i in range(d)]
  def repl(cols):
    M=[r[:] for r in B]
    for j in cols:
      for i in range(d): M[i][j]=-C[i][j]
    return det(M)
  K0=det(B)
  Kj=[repl([j]) for j in range(d)]
  Kjk=[[repl([j,k]) if j!=k else Fr(0) for k in range(d)] for j in range(d)]
  F=lambda j: z*Kj[j]
  Fjk=lambda j,k: z*z*Kjk[j][k]+(z*Kj[j] if j==k else 0)
  D=d-1; Fd=F(D); G=[-F(j)/Fd for j in range(D)]
  H=[[-(Fjk(j,k)+Fjk(j,D)*G[k]+Fjk(D,k)*G[j]+Fjk(D,D)*G[j]*G[k])/Fd for k in range(D)] for j in range(D)]
  return K0,(-1)**D*det(H)
def issq(q):
  if q<0: return False
  a,b=q.numerator,q.denominator; return isqrt(a)**2==a and isqrt(b)**2==b
def circ(d,Dset): return [[1 if (j-i)%d in Dset else 0 for j in range(d)] for i in range(d)]
if __name__=='__main__':
  for d in (3,5,7):
    A=circ(d,set(range(1,(d+1)//2))); K0,H=calH(A); r=Fr(d)/H
    print(d,K0,r,issq(r), (isqrt(r.numerator),isqrt(r.denominator)) if issq(r) else '')

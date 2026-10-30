#!/usr/bin/env python3
"""Exact finite checks of four Hanna EGF congruence pairs; no universal proof."""
from fractions import Fraction as Q
from math import comb, factorial
import hashlib

LIMIT = 300
SPEC = {
 'A232694': (4, 1, [1,1,1,9,177,5601,249681,14569545,1062623265,93853717761,9810385567329,1192614883442889,166310354311947345,26308546859152889697,4677436610087462937393,927353710845763536487305,203648424149429271943770945,49245501579619466882211194625,13045520297945193508654786790337]),
 'A232695': (5, 1, [1,1,1,11,266,10326,562926,40058076,3554828286,381374161166,48366170807276,7128626213386476,1204840675597360776,230986547885416953936,49777541426984300127816,11964954349177005321013976,3186498480002528225295506276,934756070179948684556233837476]),
 'A232696': (-1, 1, [1,1,1,-1,2,6,-264,5370,-93750,1315706,-3543736,-880688376,56549341380,-2612825765748,99009763750128,-2593891139126560,-31860555469490020,12585468136754891100,-1364794494618494224128,114095029934565534862680,-7984695190944325311086112,419424013080533747232201968]),
 'A232686': (1, 2, [1,1,2,14,196,4352,135864,5550560,283148560,17454554624,1269351110560,106925372050688,10284022937294400,1116367988620697600,135485394041604223360,18237986446534462042112,2704852975724947252789504,439428075979518197809971200,77809530527904746552436748800])}
CHOOSE = [[comb(n,j) for j in range(n+1)] for n in range(LIMIT+1)]

def exact_div(x, y):
 q,r=divmod(x,y)
 assert r==0, (x,y)
 return q

def sequence(k,p,N):
 # Lagrange inversion for B=A': E_n B=E_n A^(p*(k*n+1))/(k*n+1).
 # Powers are computed by A*(A^t)'=t*A'*A^t; all divisions checked.
 a=[1]
 for n in range(N):
  d=k*n+1
  if d==0:
   assert (k,p,n)==(-1,1,1)
   a.append(1)  # [x] A(x/B)=a(1), since B(0)=1.
   continue
  t=p*d; h=[1]
  for m in range(1,n+1):
   h.append(exact_div(sum(((t+1)*j-m)*CHOOSE[m][j]*a[j]*h[m-j]
                         for j in range(1,m+1)),m))
  a.append(exact_div(h[n],d))
 return a

def mul(f,g,N):
 return [sum((f[j]*g[n-j] for j in range(n+1)),Q(0)) for n in range(N+1)]
def power(f,k,N):
 if k<0:
  assert f[0]==1
  inv=[Q(1)]
  for n in range(1,N+1): inv.append(-sum((f[j]*inv[n-j] for j in range(1,n+1)),Q(0)))
  return power(inv,-k,N)
 out=[Q(1)]+[Q(0)]*N
 for _ in range(k): out=mul(out,f,N)
 return out

def compose(f,g,N):
 assert g[0]==0
 out=[Q(0)]*(N+1); gp=[Q(1)]+[Q(0)]*N
 for c in f[:N+1]:
  out=[x+c*y for x,y in zip(out,gp)]
  gp=mul(gp,g,N)
 return out

def rhs(f,k,p,N):
 b=[(j+1)*f[j+1] if j+1<len(f) else Q(0) for j in range(N+1)]
 z=[Q(0)]+power(b,k,N)[:N]
 return power(compose(f,z,N),p,N)

def independent_iteration(k,p,N):
 # Literal rational Q[[x]] fixed-point iteration of the defining differential equation.
 f=[Q(1),Q(1)]+[Q(0)]*N
 for _ in range(N+1):
  h=rhs(f,k,p,N)
  f=[Q(1)]+[h[n]/(n+1) for n in range(N+1)]
 return [exact_div((f[n]*factorial(n)).numerator,(f[n]*factorial(n)).denominator)
         for n in range(N+1)]

def residual(a,k,p,N):
 f=[Q(a[j],factorial(j)) if j<len(a) else Q(0) for j in range(N+2)]
 h=rhs(f,k,p,N)
 r=[((n+1)*f[n+1]-h[n])*factorial(n) for n in range(N+1)]
 assert all(z.denominator==1 for z in r)
 return [int(z) for z in r]

def v2(n):
 assert n
 n=abs(n);return (n & -n).bit_length()-1

def main():
 sequences={}
 for sid,(k,p,published) in SPEC.items():
  a=sequence(k,p,LIMIT); sequences[sid]=a
  assert a[:len(published)]==published, sid
  assert a[:13]==independent_iteration(k,p,12),sid
  if sid=='A232694': tests=[(3,3,0),(8,0,1)]
  elif sid=='A232695': tests=[(5,0,1),(2,4,0)]
  elif sid=='A232696': tests=[(2,4,0),(4,10,0)]
  else: tests=[(2,2,0),(4,4,0)]
  for modulus,start,value in tests:
   assert all(x%modulus==value for x in a[start:]), (sid,modulus)
  print(f'{sid}: fidelity {len(published)} published terms; independent Q[[x]] iteration n=0..12 agrees')
  print(f'{sid}: exactly recomputed n=0..{LIMIT}; requested congruences {tests} PASS')
  print(f'{sid}: exact terms follow')
  for n,x in enumerate(a): print(sid,n,x)
 # Finite residual data for the mathematical outlines.
 for sid,modulus,cut,N in [('A232694',3,2,10),('A232695',2,3,40),('A232696',2,3,40),('A232696',4,9,40),('A232686',2,1,8),('A232686',4,3,24)]:
  k,p,_=SPEC[sid]; f=sequences[sid][:cut+1]; r=residual(f,k,p,N)
  assert all(x%modulus==0 for x in r),(sid,modulus,r)
  print(f'RESIDUAL F\'-F(x*(F\')^{k})^{p} {sid} modulus={modulus} truncation={cut} n=0..{N}: {r}')
 # Weighted denominator hypotheses used by the infinite-tail arguments.
 for sid,cut,num,den in [('A232695',3,1,2),('A232696',3,1,2),('A232696',9,3,4)]:
  a=sequences[sid][:cut+1]
  for j,x in enumerate(a):
   c=Q(x,factorial(j)); assert v2(c.denominator) <= (num*j)//den
  for j in range(cut):
   c=Q(a[j+1],factorial(j)); assert v2(c.denominator) <= (num*j)//den
  print(f'WEIGHT {sid} truncation={cut}: F and F\' denominator v2 <= floor({num}*degree/{den}) PASS')
 # Exact exponential approximate-solution identities on a finite range.
 for k,modulus in [(4,8),(5,5)]:
  vals=[sum(comb(n,j)*(k*j)**(n-j) for j in range(n+1)) for n in range(LIMIT+1)]
  assert all(h%modulus==1 for h in vals)
  print(f'EXP_RESIDUAL k={k} modulus={modulus} n=0..{LIMIT} PASS')
 print('These are finite exact checks, not proofs of the unbounded claims.')
 print('ALL_CHECKS_OK')
if __name__=='__main__': main()

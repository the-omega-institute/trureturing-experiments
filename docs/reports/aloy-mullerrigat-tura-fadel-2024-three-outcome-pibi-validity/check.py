#!/usr/bin/env python3
"""Aloy–Müller-Rigat–Tura–Fadel arXiv:2406.11792v1, Table III: classical validity of the five proposed
three-outcome permutationally invariant Bell inequalities (3PIBIs).
(1) For 1 <= N <= NMAX (argv[1]) and every vector of counts c_{a,a'} >= 0 summing to N, evaluate the
    symmetrized observables through Table II (the LDS values) and check B_k >= 0 for k = 1..5, with the
    minimum attained (equal to 0) for every N >= 2.
(2) For N <= 4, enumerate every local deterministic strategy (each party's pair of outcomes) and compute
    the observables from the definitions P_{a|x} = sum_i p(a_i|x_i), P_{ab|xy} = sum_{i != j} p(a_i|x_i) p(b_j|y_j),
    checking agreement with (1).
(3) Check the sum-of-squares style identities of the proof route symbolically (sympy)."""
import sys, itertools
from math import comb
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 14
ROWS=[(1,1,0,-2,0,0),(1,1,-2,-2,2,0),(-2,1,2,2,0,4),(-6,1,4,4,2,12),(-6,1,4,0,0,24)]
def obs_from_counts(c):
    A=c[0][0]+c[0][1]+c[0][2]; B=c[1][0]+c[1][1]+c[1][2]; C=c[0][0]+c[1][0]+c[2][0]; D=c[0][1]+c[1][1]+c[2][1]
    P0=A+B+C+D
    P00=(A*A-A)+(C*C-C)+(B*B-B)+(D*D-D)
    P01=(A*D-c[0][1])+(B*C-c[1][0])
    P10=(A*C-c[0][0])+(B*D-c[1][1])
    P11=A*B+C*D
    return (P0,P00,P01,P10,P11)
def val(row,o): return sum(r*x for r,x in zip(row[:5],o))+row[5]
def compositions(n,k):
    if k==1: yield (n,); return
    for i in range(n+1):
        for rest in compositions(n-i,k-1): yield (i,)+rest
bad=0; mins={}
for N in range(1,NMAX+1):
    m=[None]*5
    for v in compositions(N,9):
        c=[list(v[0:3]),list(v[3:6]),list(v[6:9])]; o=obs_from_counts(c)
        for k,row in enumerate(ROWS):
            x=val(row,o); m[k]=x if m[k] is None else min(m[k],x)
    mins[N]=m
    if any(x<0 for x in m): bad+=1
    if N>=2 and any(x!=0 for x in m): bad+=1
# (2) direct definitions
for N in range(1,5):
    for strat in itertools.product(range(9),repeat=N):
        pairs=[(s//3,s%3) for s in strat]   # (outcome at x=0, outcome at x=1)
        def p(i,a,x): return 1 if pairs[i][x]==a else 0
        def P1(a,x): return sum(p(i,a,x) for i in range(N))
        def P2(a,b,x,y): return sum(p(i,a,x)*p(j,b,y) for i in range(N) for j in range(N) if i!=j)
        o=(P1(0,0)+P1(0,1)+P1(1,0)+P1(1,1),
           P2(0,0,0,0)+P2(0,0,1,1)+P2(1,1,0,0)+P2(1,1,1,1),
           P2(0,1,0,1)+P2(1,0,0,1), P2(0,0,0,1)+P2(1,1,0,1), P2(0,1,0,0)+P2(0,1,1,1))
        c=[[0]*3 for _ in range(3)]
        for a,b in pairs: c[a][b]+=1
        if o!=obs_from_counts(c): bad+=1
# (3) identities
import sympy as sp
c=sp.symbols('c00 c01 c02 c10 c11 c12 c20 c21 c22'); C3=[c[0:3],c[3:6],c[6:9]]
A=C3[0][0]+C3[0][1]+C3[0][2]; B=C3[1][0]+C3[1][1]+C3[1][2]; Cc=C3[0][0]+C3[1][0]+C3[2][0]; D=C3[0][1]+C3[1][1]+C3[2][1]
o=obs_from_counts(C3); R=A+B; S=Cc+D; K=C3[0][0]+C3[0][1]+C3[1][0]+C3[1][1]
def G(k,u,v): w=u+v; return 6*(k-1)*(k-2)+w**2+(6*k-7)*w+2*u*v
ids=[sp.expand(val(ROWS[1],o)-((R-S)**2+2*K)),
     sp.expand(2*val(ROWS[2],o)-((A-B)**2+(Cc-D)**2+(R-2)**2+(S-2)**2+4*(R-2)*(S-2)+2*(R-2)+6*(S-2)+4*(R-K))),
     sp.expand(val(ROWS[3],o)-G(K,C3[0][2]+C3[1][2],C3[2][0]+C3[2][1])),
     sp.expand(val(ROWS[4],o)-(G(C3[0][1],C3[0][0]+C3[0][2],C3[1][1]+C3[2][1])+G(C3[1][0],C3[1][1]+C3[1][2],C3[0][0]+C3[2][0])))]
if any(x!=0 for x in ids): bad+=1; print('identity residues',ids)
print(f"NMAX={NMAX} minima(N=1..4)={[mins[n] for n in range(1,5)]} identities_ok={all(x==0 for x in ids)} bad={bad}")
sys.exit(1 if bad else 0)

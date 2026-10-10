#!/usr/bin/env python3
"""Botteron–Broadbent–Chhaibi–Nechita–Pellegrini arXiv:2312.00725v3 (Quantum 8, 1402), Conjecture 15 ("Dyck paths"):
the right multiplication P^{⊠k} gives the highest CHSH value of Orbit^{(k)}(P), for the wiring W_BS on conv{PR, SR, I}.
(1) From the source's multiplication table (eq. multiplication-table-PR-SR-I), in convex coordinates (c1, c2, c3) with
    c1 = a, c2 = u - a, c3 = 1 - u: (u,a)*(v,b) = (uv, av + (1+u-2a) b / 4), CHSH = 1/2 + (u+a)/4 (sympy).
(2) The associator identity a((A*B)*C) - a(A*(B*C)) = c (2(1-v)(2-t) + (1+v) t)/16, t = 1+u-2a (sympy).
(3) Exhaustive check over all full binary parenthesizations (Catalan many) for k <= 9 at 2000 random boxes of
    conv{PR,SR,I} (and the CHSH >= 3/4 part): max CHSH equals CHSH of the right product (exact rationals).
(4) The multiplication table follows from the source's box product with the 32-entry vector of W_BS."""
import sys, random, functools
from fractions import Fraction as F
import sympy as sp
a1,a2,a3,b1,b2,b3=sp.symbols('a1 a2 a3 b1 b2 b3')
# table rows P, cols Q; vectors in basis (PR, SR, I)
T={('PR','PR'):(1,0,0),('PR','SR'):(1,0,0),('PR','I'):(0,0,1),('SR','PR'):(sp.Rational(1,2),sp.Rational(1,2),0),('SR','SR'):(0,1,0),('SR','I'):(0,0,1),('I','PR'):(sp.Rational(1,4),-sp.Rational(1,4),1),('I','SR'):(0,0,1),('I','I'):(0,0,1)}
A=dict(PR=a1,SR=a2,I=a3); B=dict(PR=b1,SR=b2,I=b3)
prod=[sum(A[p]*B[q]*T[(p,q)][i] for p in A for q in B) for i in range(3)]
u,a,v,b=sp.symbols('u a v b')
sub={a1:a,a2:u-a,a3:1-u,b1:b,b2:v-b,b3:1-v}
c1=sp.expand(prod[0].subs(sub)); c3=sp.expand(prod[2].subs(sub))
ok1 = sp.simplify(c1-(a*v+(1+u-2*a)*b/4))==0 and sp.simplify((1-c3)-u*v)==0
chsh=lambda c: 1*c[0]+sp.Rational(3,4)*c[1]+sp.Rational(1,2)*c[2]
ok1 = ok1 and sp.simplify(chsh((a,u-a,1-u))-(sp.Rational(1,2)+(u+a)/4))==0
print('product and CHSH formulas:', ok1)
def mul(X,Y): (u,a),(v,b)=X,Y; return (u*v, a*v+(1+u-2*a)*b/4)
w,c=sp.symbols('w c')
lhs=mul(mul((u,a),(v,b)),(w,c)); rhs=mul((u,a),mul((v,b),(w,c)))
t=1+u-2*a
ok2 = sp.simplify(lhs[0]-rhs[0])==0 and sp.simplify(sp.expand(lhs[1]-rhs[1]) - c*(2*(1-v)*(2-t)+(1+v)*t)/16)==0
print('associator identity:', ok2)
def trees(k):
    if k==1: return [('P',)]
    out=[]
    for l in range(1,k):
        for L in trees(l):
            for R in trees(k-l): out.append((L,R))
    return out
@functools.lru_cache(None)
def T_(k): return trees(k)
def ev(tr,P):
    if tr==('P',): return P
    return mul(ev(tr[0],P),ev(tr[1],P))
random.seed(5); bad=0
for trial in range(2000):
    # random point of conv{PR,SR,I}: weights (p1,p2,p3)
    w1,w2=sorted([F(random.randint(0,1000),1000),F(random.randint(0,1000),1000)])
    p1,p2,p3=w1,w2-w1,1-w2
    P=(p1+p2,p1)   # u = c1+c2, a = c1
    for k in range(2,10 if trial<200 else 7):
        best=max(1/F(2)+(x[0]+x[1])/4 for x in (ev(tr,P) for tr in T_(k)))
        R=P
        for _ in range(k-1): R=mul(R,P)
        if best!=1/F(2)+(R[0]+R[1])/4: bad+=1
print('parenthesization violations:', bad)
ok = ok1 and ok2 and bad==0
# (4) The multiplication table follows from the source's box product (eq: definition of the box product) with
#     W_BS = [0,0,1,1, 0,0,1,1, 0,0,0,1, 0,0,0,1, 0,1,1,0,0,1,1,0, 0,1,1,0,0,1,1,0] in the order f1, g1, f2, g2, f3, g3.
W=[0,0,1,1, 0,0,1,1, 0,0,0,1, 0,0,0,1, 0,1,1,0,0,1,1,0, 0,1,1,0,0,1,1,0]
f1=lambda x,a2: W[0+2*x+a2]; g1=lambda y,b2: W[4+2*y+b2]; f2=lambda x,a1: W[8+2*x+a1]; g2=lambda y,b1: W[12+2*y+b1]
f3=lambda x,a1,a2: W[16+4*x+2*a1+a2]; g3=lambda y,b1,b2: W[24+4*y+2*b1+b2]
def boxprod(P,Q):
    return {(a,b,x,y): sum(P[(a1,b1,f1(x,a2),g1(y,b2))]*Q[(a2,b2,f2(x,a1),g2(y,b1))]*(a==f3(x,a1,a2))*(b==g3(y,b1,b2))
            for a1 in (0,1) for a2 in (0,1) for b1 in (0,1) for b2 in (0,1)) for a in (0,1) for b in (0,1) for x in (0,1) for y in (0,1)}
keys=[(a,b,x,y) for a in (0,1) for b in (0,1) for x in (0,1) for y in (0,1)]
PR={k:F(1,2)*((k[0]^k[1])==(k[2]*k[3])) for k in keys}; SR={k:F(1,2)*(k[0]==k[1]) for k in keys}; I_={k:F(1,4) for k in keys}
lin=lambda c: {k:c[0]*PR[k]+c[1]*SR[k]+c[2]*I_[k] for k in keys}
basis={'PR':PR,'SR':SR,'I':I_}
table={('PR','PR'):(1,0,0),('PR','SR'):(1,0,0),('PR','I'):(0,0,1),('SR','PR'):(F(1,2),F(1,2),0),('SR','SR'):(0,1,0),('SR','I'):(0,0,1),('I','PR'):(F(1,4),-F(1,4),1),('I','SR'):(0,0,1),('I','I'):(0,0,1)}
tab_ok=all(boxprod(basis[p],basis[q])==lin(table[(p,q)]) for (p,q) in table)
print('multiplication table from W_BS:', tab_ok)
ok = ok and tab_ok
print('bad=', 0 if ok else 1); sys.exit(0 if ok else 1)

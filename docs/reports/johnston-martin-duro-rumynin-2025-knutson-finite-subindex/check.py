#!/usr/bin/env python3
"""Johnston-Martin Duro-Rumynin (arXiv:2403.18498v6; Comm. Algebra 53 (2025)), Question 2.5: is there a finite S with
K_r(x) = K_r^S(x) for all x?  Counterexample ring R = {(c, z) in Z x B : c = phi(z) mod 4}, B = Z[u,v]/(u^2, v^2),
phi(u) = phi(v) = 2; alpha(c, z) = c; r = (4, 0). Elements are (c, d, a, b, f) for (c, d + a u + b v + f uv).
Checks: membership and closure of the basis; r regular (s r = alpha(s) r); x_n y_n = r for even n (so K_r(x_n) = 1);
for random finite S (including elements eligible for other x_m), at n = 2(1 + sum |b(s)|) every eligible index is even,
so K_r^S(x_n) != 1; and K_r^S(x_m) = 1 at the m whose partner y_m lies in S. Exit 0 and RESULT PASS when all hold.
"""
import random, math, sys
def mem(x): c,d,a,b,f=x; return (c-d-2*a-2*b)%4==0
def mul(x,y):
    c1,d1,a1,b1,f1=x; c2,d2,a2,b2,f2=y
    return (c1*c2, d1*d2, d1*a2+a1*d2, d1*b2+b1*d2, d1*f2+f1*d2+a1*b2+b1*a2)
r=(4,0,0,0,0)
def is_mult_of_r(z):   # z = m r  <=> z = (4m,0,0,0,0)
    return z[1:]==(0,0,0,0) and z[0]%4==0, (z[0]//4 if z[0]%4==0 else None)
def basis_elems(B):
    # random elements of R via basis e=(1,1,0,0,0), p=(2,0,1,0,0), q=(2,0,0,1,0), w=(0,0,0,0,1), r
    e=(1,1,0,0,0); p=(2,0,1,0,0); q=(2,0,0,1,0); w=(0,0,0,0,1)
    out=[]
    for _ in range(B):
        co=[random.randint(-6,6) for _ in range(5)]
        out.append(tuple(co[0]*e[i]+co[1]*p[i]+co[2]*q[i]+co[3]*w[i]+co[4]*r[i] for i in range(5)))
    return out
random.seed(2403); bad=0
for n in range(2,60,2):
    xn=(2,0,1,n,0); yn=(2,0,1,-n,0)
    assert mem(xn) and mem(yn) and mul(xn,yn)==r, n
# regularity: s r = alpha(s) r for random s
for s in basis_elems(500):
    assert mem(s) and mul(s,r)==(4*s[0],0,0,0,0)
# for random finite S, check that for n > 2(1+sum|b|) every eligible index is even
for trial in range(300):
    S=basis_elems(random.randint(1,12))
    n=2*(1+sum(abs(s[3]) for s in S))
    xn=(2,0,1,n,0); idx=[]
    for s in S:
        ok,m=is_mult_of_r(mul(s,xn))
        if ok: idx.append(m)
    g=0
    for m in idx: g=math.gcd(g,m)
    if idx and g==1: bad+=1; print('fail',S,n,idx)
print('membership/product/regularity ok; trials with odd S-index:',bad)
# stronger: S contains eligible-type elements y_{n'} = (2,0,1,-n',0), r, and (2k, 0, a, -n'a, f) variants
bad=0; witnessed=0
for trial in range(300):
    ns=random.sample(range(2,40,2),random.randint(1,6))
    S=[(2,0,1,-m,0) for m in ns]+[r]+[(2*random.randint(-3,3)+2*a % 4,0,a,-m*a,random.randint(-3,3)) for m in ns for a in (1,2,3)]
    S=[s for s in S if mem(s)]
    n=2*(1+sum(abs(s[3]) for s in S))
    xn=(2,0,1,n,0); idx=[]
    for s in S:
        ok,m=is_mult_of_r(mul(s,xn))
        if ok: idx.append(m)
    g=0
    for m in idx: g=math.gcd(g,m)
    if idx: witnessed+=1
    if idx and g==1: bad+=1
    # sanity: at n' in ns the S-index is 1 because y_{n'} in S
    m0=ns[0]; x0=(2,0,1,m0,0); idx0=[is_mult_of_r(mul(s,x0))[1] for s in S if is_mult_of_r(mul(s,x0))[0]]
    g0=0
    for m in idx0: g0=math.gcd(g0,m)
    assert g0==1,(ns,idx0)
print('trials with eligible elements:',witnessed,' trials with S-index 1 at the large n:',bad)
print('RESULT PASS' if bad==0 else 'RESULT FAIL'); sys.exit(0 if bad==0 else 1)

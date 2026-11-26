"""Can the available Higgs and chiral fields themselves induce the canonical product-group trace?"""
from fractions import Fraction as F
from itertools import product
from math import prod

# Rows are twice the Dynkin index in (X,4,L,R,S); F is orthogonal to hypercharge.
rows={
    "left":(0,6,12,0,0),"right":(0,6,0,12,0),"PS scalar":(0,6,0,12,0),
    "X-S scalar":(5,0,0,0,60),"R-S scalar":(0,0,0,16,30),
    "bidoublet":(0,0,2,2,0),"family adjoints":(0,0,0,0,0)}
def missing(v):
    x,c,l,r,s=v
    return 2*c-l-r-F(16,3)*x-F(4,9)*s
assert all(missing(v)<=0 for v in rows.values())
assert all(missing(rows[k])<0 for k in ("X-S scalar","R-S scalar","bidoublet"))
base=tuple(map(sum,zip(*rows.values())))
assert base==(5,18,14,42,90)
print("unit-trace levels, degree of missing linear combination",base,missing(base))
# Only add single-factor representations of X, 4, L, or S. Minimize total added color index
# and, at that color index, the dimension of pure-factor trace supplements.
choices=[]
for dx,ds in product(range(10),range(10)):
    c=F(16,3)*(base[0]+dx)+base[3]+F(4,9)*(base[4]+ds)
    if c.denominator==1:
        choices.append((int(c),dx,ds))
target,dx,ds=min(choices)
assert (target,dx,ds)==(110,0,3)
def minimal(total,blocks):
    best={0:(0,())}
    for k in range(total):
        if k not in best:
            continue
        cost,items=best[k]
        for index,dim,label in blocks:
            n=k+index
            if index and n<=total and (n not in best or cost+dim<best[n][0]):
                best[n]=cost+dim,items+(label,)
    return best[total]
colors=[]
for u,v,w in product(range(12),repeat=3):
    if u<w:
        continue
    labels=(u,v,w)
    d=prod(sum(labels[k:l])+l-k for k in range(4) for l in range(k+1,4))//12
    c2=F(3*u*u+4*v*v+3*w*w+4*u*v+2*u*w+4*v*w+12*u+16*v+12*w,8)
    index=F(2*d,15)*c2
    assert index.denominator==1
    if 1<=index<=target-base[1]:
        colors.append((int(index),d,labels))
weak=[(n*(n+1)*(n+2)//6,n+1,n) for n in range(1,14)]
C=minimal(target-base[1],colors)
L=minimal(target-base[2],weak)
S=minimal(ds,weak)
print("least color target in the pure-factor repair",target)
print("minimum dimensions and Dynkin labels, color / L / S",C,L,S)

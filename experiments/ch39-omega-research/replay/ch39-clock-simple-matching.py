"""Canonical hypercharge, integral embedding levels, and the lowest color level in SU(N)."""
from itertools import product
from fractions import Fraction as F
from math import prod

choices=[]
for ix in range(1,9):
    for ir in range(1,42):
        for iss in range(1,100):
            i4=F(16,3)*ix+ir+F(4,9)*iss
            if i4.denominator==1 and i4<=41:
                choices.append((int(i4),ix,ir,iss))
assert min(choices)==(9,1,1,6)
assert [v for v in choices if v[0]==9]==[(9,1,1,6)]
# The two index-one factors X and R each require a separate fundamental block.
# Remaining active indices: color 9, L 9, S 6. Family needs at least one
# nontrivial center charge; family conjugation leaves dimension and index costs unchanged.
colors=((1,0),(4,1),(6,2),(10,6),(15,8))
weak=((1,0),(2,1),(3,4))
targets=(9,9,6,1)
blocks=[]
for a,b,c,f in product(colors,weak,weak,(1,3)):
    dims=(a[0],b[0],c[0],f)
    n=prod(dims)
    cost=tuple(r[1]*(n//r[0]) for r in (a,b,c))+(int(f==3),)
    if any(cost) and all(x<=y for x,y in zip(cost,targets)):
        blocks.append((dims,cost,n))
states={(0,0,0,0):(0,())}
# Increasing total charge grade makes this dynamic program acyclic even when family is saturated.
for i,j,k,f in product(range(10),range(10),range(7),range(2)):
    state=(i,j,k,f)
    if state not in states:
        continue
    n,words=states[state]
    for dims,d,extra in blocks:
        dest=(i+d[0],j+d[1],k+d[2],min(1,f+d[3]))
        if dest==state or any(x>y for x,y in zip(dest,targets)):
            continue
        if dest not in states or n+extra<states[dest][0]:
            states[dest]=(n+extra,words+(dims,))
N,witness=states[targets]
print("lowest color embedding level",min(choices),"minimum SU(N) at this level",N+5)
print("completing blocks (color, L, S, family dimensions)",witness)

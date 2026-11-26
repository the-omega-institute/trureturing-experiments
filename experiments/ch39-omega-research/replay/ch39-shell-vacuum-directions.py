"""Zero-energy realizations of a positive alignment action built from orbit centers and chord directions."""
from itertools import permutations,combinations,product
from collections import Counter
from pathlib import Path
import runpy

g=runpy.run_path(str(Path(__file__).with_name('ch39-shell-family-alignment.py')))
X,T=g["sites"],g["shell"](867)
dot,cross=g["dots"],g["cross"]
swaps=((0,1),(1,2),(0,2))
def swap(x,j,k):
    y=list(x)
    y[j],y[k]=y[k],y[j]
    return tuple(y)
indices=[[X.index(swap(x,*pair)) for x in X] for pair in swaps]
def sub(x,y):
    return [a-b for a,b in zip(x,y)]
def scale(a,x):
    return [a*b for b in x]
def parallel(x,y):
    return cross(x,y)==[0,0,0] and dot(x,y)>0
def valid_edge(i,j,u,seq,offset):
    k=indices[j][i]
    if k>=i:
        return True
    a,b=swaps[j]
    dx=X[i][a]-X[i][b]
    return parallel(scale(dx,u),sub(seq[i-offset],seq[k-offset]))
def aligned(seq,offset,axes,z):
    total=[sum(x[d] for x in seq) for d in range(3)]
    return parallel(z,total) and all(valid_edge(i,j,u,seq,offset)
        for i in range(offset,offset+len(seq)) for j,u in enumerate(axes))
def reflected(x,u):
    d=dot(u,u)
    nums=sub(scale(d,x),scale(2*dot(x,u),u))
    if any(n%d for n in nums):
        return None
    return tuple(n//d for n in nums)
def traverse_outer(first,axes):
    seq={0:first}
    stack=[0]
    while stack:
        i=stack.pop()
        for j,u in enumerate(axes[:2]):
            nxt=indices[j][i]
            pt=reflected(seq[i],u)
            if pt not in T or (nxt in seq and seq[nxt]!=pt):
                return ()
            if nxt not in seq:
                seq[nxt]=pt
                stack.append(nxt)
    return tuple(seq[i] for i in range(6))
solutions=set()
for x,y,w in permutations(T,3):
    u=scale(-1,sub(x,y))
    v=scale(-1,sub(x,w))
    c=dot(u,v)
    if c>=0 or 4*c*c!=dot(u,u)*dot(v,v):
        continue
    axes=(u,v,sub(scale(dot(u,u),v),scale(2*c,u)))
    z=cross(u,v)
    outer=traverse_outer(x,axes)
    if len(set(outer))!=6 or not aligned(outer,0,axes,z):
        continue
    for inner in permutations(sorted(set(T)-set(outer)),3):
        if aligned(inner,6,axes,z):
            solutions.add(outer+inner)
shapes=Counter((tuple(sorted(seq[0])),tuple(sorted(seq[6]))) for seq in solutions)
print("oriented zero-energy injections",len(solutions),"shapes and counts",sorted(shapes.items()))
assert len(solutions)==12 and len(shapes)==4 and set(shapes.values())=={3}
assert tuple(g["mapped"]) in solutions

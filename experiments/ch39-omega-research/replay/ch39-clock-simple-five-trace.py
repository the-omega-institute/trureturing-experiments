"""Fifth-trace obstruction for all minimum-rank, level-nine SU(N) defining representations."""
from itertools import product,combinations,combinations_with_replacement
from functools import lru_cache
from pathlib import Path
import runpy
from flint import fmpz_mat

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-simple-matching.py")))
blocks=sorted(m["blocks"])
solutions=[]
def search(start,state,N,used):
    if state==(9,9,6,1):
        if N==33:
            solutions.append(used)
        return
    for i in range(start,len(blocks)):
        dims,d,size=blocks[i]
        dest=tuple(x+y for x,y in zip(state[:3],d[:3]))+(min(1,state[3]+d[3]),)
        if dest==state or N+size>33 or any(x>y for x,y in zip(dest,m["targets"])):
            continue
        # A prefix of a dimension-minimal solution is minimal for its state.
        if m["states"][dest][0]!=N+size:
            continue
        search(i,dest,N+size,used+(dims,))
search(0,(0,0,0,0),0,())
def spectrum(t,choice,signs):
    a,b,w,h=t
    color=(a+h,b+h,-a-b+h,-3*h)
    eigen=[4*h*signs[0],4*h*signs[0],-8*h*signs[0],3*h,-3*h]
    conjugations=iter(signs[1:])
    for dim,L,S,F in choice:
        if dim in (4,10):
            z=next(conjugations)
        if dim==1:
            e=[0]
        elif dim==4:
            e=[z*v for v in color]
        elif dim==6:
            e=[x+y for x,y in combinations(color,2)]
        elif dim==10:
            e=[z*(x+y) for x,y in combinations_with_replacement(color,2)]
        else:
            assert dim==15
            e=[0]*3+[color[i]-color[j] for i,j in product(range(4),repeat=2) if i!=j]
        eigen.extend([e0+wk*w+sk*2*h for e0 in e for wk in range(1-L,L,2)
                      for sk in range(1-S,S,2)]*F)
    assert len(eigen)==38 and sum(eigen)==0
    return eigen
def light(t):
    a,b,w,h=t
    c=[a,b,-a-b]
    eigen=[z+e*w+h for z in c for e in (-1,1)]+[-z+y*h for z in c for y in (-4,2)]
    eigen.extend([-w-3*h,w-3*h,6*h,0])
    return 3*sum(z**5 for z in eigen)
samples=((1,2,0,0),(1,2,1,0),(1,2,0,1),(1,2,1,1),(1,0,2,1),(0,0,1,1),(0,0,0,1))
orientations=0
for choices in solutions:
    for signs in product((-1,1),repeat=1+sum(d[0] in (4,10) for d in choices)):
        orientations+=1
        trace=[]
        for t in samples:
            e=spectrum(t,choices,signs)
            trace.append([sum(x**5 for x in e),sum(x**3 for x in e)*sum(x*x for x in e),light(t)])
        assert fmpz_mat(trace).rank()==3
print("minimum embedding decompositions and nonredundancy-unrestricted color/X orientations",len(solutions),orientations)
assert solutions and orientations>0

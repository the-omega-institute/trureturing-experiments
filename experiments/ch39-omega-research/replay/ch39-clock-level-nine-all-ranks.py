"""Fifth-trace test of every unitary rank at the minimum color embedding level."""
from pathlib import Path
from functools import lru_cache
from itertools import product
from flint import fmpz_mat
import runpy

prev=runpy.run_path(str(Path(__file__).with_name("ch39-clock-simple-five-trace.py")))
base=prev["m"]
samples=prev["samples"]
y=tuple(prev["light"](t) for t in samples)
blocks=[]
# On the SM Lie algebra any family spectator just makes copies of this block.
for dims,charges,_ in base["blocks"]:
    if dims[-1]!=1:
        continue
    for orientation in ((1,-1) if dims[0] in (4,10) else (1,)):
        blocks.append((dims,orientation,charges[:3]))
def get_rows(dims,sign):
    result=[]
    for t in samples:
        a,b,w,h=t
        ca=(a+h,b+h,-a-b+h,-3*h)
        d,L,S,_=dims
        if d==1:
            spectrum=[0]
        elif d==4:
            spectrum=[sign*z for z in ca]
        elif d==6:
            spectrum=[x+y for x,y in prev["combinations"](ca,2)]
        elif d==10:
            spectrum=[sign*(x+y) for x,y in prev["combinations_with_replacement"](ca,2)]
        else:
            spectrum=[0]*3+[x-y for i,x in enumerate(ca) for j,y in enumerate(ca) if i!=j]
        e=[x+j*w+k*2*h for x in spectrum for j in range(1-L,L,2) for k in range(1-S,S,2)]
        result.append((sum(x**3 for x in e),sum(x**5 for x in e)))
    return result
col=[get_rows(d,s) for d,s,_ in blocks]
# Quadratic form from the fixed levels is common to all embeddings.
def second(t):
    a,b,w,h=t
    return 9*(a*a+b*b+(a+b)**2+2*w*w+30*h*h)
# Isolate the pure-X and pure-R fundamental blocks. One may conjugate all to fix X orientation.
# X exponents are (4,4,-8), so traces are -384 h^3, -30720 h^5.
fixed=[((2*4**3-8**3)*t[3]**3,(2*4**5-8**5)*t[3]**5) for t in samples]
assert all(sum(x*x for x in prev["spectrum"](t,prev["solutions"][0],[1,1]))==second(t) for t in samples)
sec=[second(t) for t in samples]
@lru_cache(None)
def prefixes(idx,remaining):
    if not any(remaining):
        return 1
    if idx==len(blocks):
        return 0
    cost=blocks[idx][2]
    maximum=min(r//v for r,v in zip(remaining,cost) if v)
    return sum(prefixes(idx+1,tuple(r-j*v for r,v in zip(remaining,cost))) for j in range(maximum+1))
tested=0
survivors=[]
def check(idx,remaining,values,seq):
    global tested
    if not any(remaining):
        tested+=1
        tr=fmpz_mat([[v5,v3*s,z] for (v3,v5),s,z in zip(values,sec,y)])
        if tr.rank()<3:
            survivors.append(seq)
        return
    cost=blocks[idx][2]
    maximum=min(r//v for r,v in zip(remaining,cost) if v)
    for j in range(maximum+1):
        nxt=tuple(r-j*v for r,v in zip(remaining,cost))
        if prefixes(idx+1,nxt):
            vals=[(x+j*a,y+j*b) for (x,y),(a,b) in zip(values,col[idx])]
            check(idx+1,nxt,vals,seq+((idx,j),) if j else seq)
start=(9,9,6)
check(0,start,fixed,())
print("all-rank color-level-nine spectra modulo global conjugation",tested,
      "passing seven-point fifth-trace test",len(survivors))
assert tested==1584 and not survivors

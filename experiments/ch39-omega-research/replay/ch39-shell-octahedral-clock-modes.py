"""All-height balance for the rank-three elementary modification encoded by the two actual shell clocks."""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import arb,acb

home=Path(__file__).resolve().parent
m=run_path(str(home/"ch39-shell-third-shell.py"))
omega=acb(-arb(1)/2,arb(3).sqrt()/2)
v=(0,1,5)
ends=((3,1,4),(1,4,3))
assert [26-m["dot"](v,p) for p in ends]==list(m["shells"][26][0][1:])
zplus,zminus=arb(2).sqrt(),-1/arb(2).sqrt()


def comps(s,end):
    X,Y,Z=[coord+s*(dest-coord) for coord,dest in zip(v,end)]
    return (acb(X+Y+Z),X+omega*Y+omega.conjugate()*Z,X+omega.conjugate()*Y+omega*Z)


def minors(s,t):
    u,vv=comps(s,ends[0]),comps(t,ends[1])
    return [zminus*u[j]*vv[(j-1)%3]-zplus*vv[j]*u[(j-1)%3] for j in range(3)]


ranges=[]
for j in range(16):
    row=[]
    for k in range(16):
        s=arb(j).union(j+1)/16
        t=arb(k).union(k+1)/16
        row+= [abs(x).lower() for x in minors(s,t)]
    ranges.append(min(row))
assert min(ranges)>3
print("uniform lower on the absolute 2-by-2 evaluation minors with raw coordinate rows",min(ranges).str(13))
for n in (25,255,719):
    h=acb.zeta_zero(n).imag
    p,q=(h*f for f in m["shells"][26][2])
    vals=minors(p.sin()**2,q.sin()**2)
    print("zero",n,"minor magnitudes",*[abs(a).str(10) for a in vals])

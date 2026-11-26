"""Facewise Yang-Mills logarithms for the smooth S(U(4) x U(3)) lift and its chiral index."""
from fractions import Fraction as F
from itertools import product,permutations
from pathlib import Path
from runpy import run_path

parent=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-conjugates.py")))
field=parent["equivariant_field"](12)
faces=parent["faces"]


def traceless_cost(n,total):
    residue=total%n
    return F(residue*(n-residue),n)


def cost(k,face,c):
    # trace_A=k=-trace_B. Each entry is measured in units 2*pi.
    shift=(1,0) if face=="cap" else (-1,-1)
    return traceless_cost(4,k-shift[0])+traceless_cost(3,-k-shift[1])+c*k*k


def minima(face,c):
    vals={k:cost(k,face,c) for k in range(-6,7)}
    return [k for k,v in vals.items() if v==min(vals.values())]


def scope(face,chosen):
    lo=F(0)
    hi=None
    for k in range(-6,7):
        if k==chosen:
            continue
        intercept=cost(k,face,0)-cost(chosen,face,0)
        slope=k*k-chosen*chosen
        if slope==0:
            assert intercept>0
        elif slope>0:
            lo=max(lo,-intercept/slope)
        else:
            v=-intercept/slope
            hi=v if hi is None else min(hi,v)
    return lo,hi


cap,side=-3,-1
regimes=[scope("cap",cap),scope("side",side)]
lo=max(x for x,y in regimes)
hi=min(y for x,y in regimes if y is not None)
assert (lo,hi)==(F(1,36),F(1,12))
print("open determinant stiffness interval",lo,hi,
      "endpoint face minima",[(minima("cap",c),minima("side",c)) for c in (lo,hi)])
for face,k in (("cap",cap),("side",side)):
    assert minima(face,F(1,18))==[k]
assert minima("cap",F(7,12))==[0] and minima("side",F(7,12))==[-1]
assert -(2*cap+6*side)==12
sides=[]
charges=[]
for f in faces:
    seq=list(zip(f,f[1:]+f[:1]))
    turn=sum(field[x,y] for x,y in seq)%12
    k=cap if turn==3 else side
    if turn==3:
        a,b=F(-3,4),[F(1)]*3
    else:
        assert turn==5
        # This side face has one triangle vertex not lying on its D=1 edge.
        short=[(x,y) for x,y in seq if 26-parent["dot"](x,y)==1]
        assert len(short)==1
        v=next(v for v in f if v not in short[0])
        omit=v.index(3)
        sides.append(omit)
        a,b=F(-1,4),[F(-1,3) if j==omit else F(2,3) for j in range(3)]
    assert 4*a==k==-sum(b)
    assert sum((x-F(-k,3))**2 for x in b)+k*k*F(1,18)==cost(k,"cap" if turn==3 else "side",F(1,18))
    charges.append([a+x for x in b])
    for x in charges[-1]:
        assert (12*x-turn)%12==0
assert sorted(sides)==[0,0,1,1,2,2]
assert [sum(v) for v in zip(*charges)]==[1,1,1]
splits=set()
for n1 in range(7):
    for n2 in range(7-n1):
        splits.add(tuple(sorted((3-n1,3-n2,3-(6-n1-n2)))))
counts={(sum(max(0,v) for v in s),sum(max(0,-v) for v in s)) for s in splits}
assert counts=={(3,0),(4,1),(5,2),(6,3)}
print("balanced mixed degrees",list(map(sum,zip(*charges))),"all Cartan minima zero-mode counts",sorted(counts))

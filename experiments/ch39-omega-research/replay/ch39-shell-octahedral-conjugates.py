"""Conjugate center charges on the octahedral core and across the full twelve-site flag complex."""
from itertools import combinations,permutations,product
from pathlib import Path
from runpy import run_path
from flint import arb

parent=run_path(str(Path(__file__).with_name("ch39-shell-octahedral-flux.py")))
core,faces,phases,_=parent["main"]()
dot,cross=parent["dot"],parent["cross"]
pi=arb.pi()


def degree(ks,m,charge):
    turns=[]
    for k,p in zip(ks,phases):
        x=(charge*k)%m
        if 2*x>m or (2*x==m and p>0):
            x-=m
        assert -pi<p+2*pi*x/m<pi
        turns.append(x)
    assert sum(turns)%m==0
    return sum(turns)//m


def pairs(m,invariant=False):
    # Face charges can be assigned arbitrarily subject to total charge zero.
    # Keep integer branch sums rather than enumerating m^8 tuples.
    blocks=[list(range(8))]
    if invariant:
        groups={}
        for i,f in enumerate(faces):
            typ=tuple(sorted(26-dot(p,q) for p,q in combinations(f,2)))
            groups.setdefault(typ,[]).append(i)
        blocks=list(groups.values())
    else:
        blocks=[[j] for j in range(8)]
    states={(0,0)}
    for group in blocks:
        nxt=set()
        for k in range(m):
            vals=[]
            for sign in (1,-1):
                n=(sign*k)%m
                vals.append(sum(n-m if 2*n>m or (2*n==m and phases[j]>0) else n for j in group))
            nxt.update((a+vals[0],b+vals[1]) for a,b in states)
        states=nxt
    return {(a//m,b//m) for a,b in states if a%m==0 and b%m==0}


def equivariant_field(m):
    assert m in (8,12)
    # Put +1 on the cyclic sides of each oriented equilateral cap, a half turn on
    # each unit-length class, and no center field on the remaining edges.
    field={}
    for f in faces:
        if all(26-dot(x,y)==7 for x,y in combinations(f,2)):
            for x,y in zip(f,f[1:]+f[:1]):
                field[x,y]=1
                field[y,x]=-1
    for f in faces:
        for x,y in zip(f,f[1:]+f[:1]):
            if (x,y) not in field:
                field[x,y]=field[y,x]=m//2 if 26-dot(x,y)==1 else 0
    for (x,y),k in field.items():
        for p in permutations(range(3)):
            xp=tuple(x[i] for i in p)
            yp=tuple(y[i] for i in p)
            assert (field[xp,yp]-k)%m==0
    return field


def pullback(m,field):
    S=parent["shells"][26]
    def ret(p):
        if p in core:
            return p
        mapping={0:1,1:3,5:4}
        return tuple(mapping[j] for j in p)
    def ok(vs):
        return all(26-dot(v,w)<=7 for v,w in combinations(vs,2))
    signed={}
    for f in combinations(sorted(S),3):
        if not ok(f):
            continue
        a,b,c=f
        den=arb(26).sqrt()*(26+dot(a,b)+dot(b,c)+dot(c,a))
        phase=(arb(dot(a,cross(b,c)))/den).atan()
        assert abs(phase)<pi/6
        image=list(map(ret,f))
        assert ok(image)
        k=sum(field.get((u,v),0) for u,v in zip(image,image[1:]+image[:1]))%m
        k=(k+m//2)%m-m//2
        assert abs(phase)+2*pi*abs(k)/m<pi
        signed[f]=k
    for t in combinations(sorted(S),4):
        if ok(t):
            assert sum((-1)**j*signed[t[:j]+t[j+1:]] for j in range(4))==0
    vals=[sum(field[u,v] for u,v in zip(f,f[1:]+f[:1]))%m for f in faces]
    assert degree(vals,m,1)==3 and degree(vals,m,-1)==-3
    print("strict permutation-invariant edge background, order",m,"face charges",sorted(vals))


def main():
    for m in range(1,9):
        allpairs=pairs(m)
        strict_shape=pairs(m,True)
        assert all(abs(a+b)<=4 for a,b in allpairs)
        print("order",m,"degrees (3,-3), arbitrary / invariant face residues",
              (3,-3) in allpairs,(3,-3) in strict_shape)
        assert ((3,-3) in allpairs)==(m in (5,7,8))
        assert ((3,-3) in strict_shape)==(m==8)
        if m==4:
            assert {b for a,b in allpairs if a==3}=={1}
    for mod in (8,12):
        pullback(mod,equivariant_field(mod))


if __name__=="__main__":
    main()

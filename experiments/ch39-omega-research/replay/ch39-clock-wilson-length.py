"""Shortest positive adjoint Wilson interaction able to reverse both adjacent-zero energy comparisons."""
from itertools import combinations,permutations
from collections import defaultdict
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-spin-clock-neighbor.py")))
q,sh = m["q"],m["sh"]
nodes = sh["sites"]
links = {(i,j):(d,z) for i,j,d,z in sh["edges"]}
I = m["m"]["I"]
ns = (718,719,720)
turns = [m["rotations"](n)[0] for n in ns]
def key(word):
    # Exact reduction in the free group; (1,2) designate U,V and minus is inverse.
    stack = []
    for v in word:
        if stack and stack[-1]==-v:
            stack.pop()
        else:
            stack.append(v)
    while len(stack)>=2 and stack[0]==-stack[-1]:
        stack=stack[1:-1]
    return tuple(stack)

from functools import lru_cache
@lru_cache(None)
def traces(word):
    vals = []
    for T in turns:
        rot = I
        for v in word:
            h = T[(724,804)[abs(v)-1]]
            rot = rot*(h if v>0 else h.conjugate().transpose())
        tr = rot[0,0]+rot[1,1]
        assert tr.imag.contains(0)
        vals.append(tr.real)
    return vals

def cycles(n):
    for group in combinations(range(19),n):
        for rest in permutations(group[1:]):
            if rest[0]>rest[-1]:
                continue
            loop = (group[0],*rest)
            Z = q(1)
            word,ds=[],[]
            for e in zip(loop,loop[1:]+loop[:1]):
                pair=tuple(sorted(e))
                if pair not in links:
                    break
                d,z=links[pair]
                ds.append(d)
                Z*=z if e==pair else z.conjugate()
                if d in m["axes"]:
                    word.append((1 if d==724 else 2)*(1 if e==pair else -1))
            else:
                yield loop,tuple(sorted(ds)),Z.real,key(word)

def gaps(z,word):
    ts = traces(word)
    return [-(ts[i]-ts[1])*(ts[i]+ts[1]+2*z) for i in (0,2)]

def check():
    sums = []
    for _,_,z,word in cycles(3):
        if word:
            diff=gaps(z,word)
            sums.append(sum(diff).upper())
    assert max(sums)<0
    print("active triangles",len(sums),"two-neighbor sum upper",max(sums).str(10))
    groups=defaultdict(list)
    for _,ds,z,word in cycles(4):
        groups[ds].append(gaps(z,word) if word else (q(0),q(0)))
    keep=[]
    for ds,pairs in groups.items():
        means=[sum(v)/len(pairs) for v in zip(*pairs)]
        if all(v>0 for v in means):
            keep.append(ds)
            print("selected length class",ds,"count",len(pairs),"mean gaps",*(v.str(12) for v in means))
            assert means[0]>q("2.01821") and means[1]>q("2.66987")
        else:
            assert any(v<0 or v.is_zero() for v in means)
    assert keep==[(244,244,724,804)]
    # The original spin and neutral determinations supply uniform lower bounds too.
    d718=3*(q("6.86775")-q("9.37987"))+q(".28869")
    d720=3*(q("5.97493")-q("9.37987"))-q(".50872")
    assert d718+q("4.1")*q("2.01821")>1 and d720+q("4.1")*q("2.66987")>q(".2")
    print("after adding 4.1e-12 in Q^4 units, positive gaps in 1e-12 Q^4",
          (d718+q("4.1")*q("2.01821")).lower().str(10),
          (d720+q("4.1")*q("2.66987")).lower().str(10))

if __name__=="__main__":
    check()

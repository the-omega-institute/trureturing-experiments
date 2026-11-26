"""Recoverability of the successful order-nine reference from the four squared-chord differences."""
from itertools import product
from math import comb
chords=((6,12,14),(244,724,804))
goal=(2,1,4,6)
targets={tuple(s[i//2]*z%9 for i,z in enumerate(goal)) for s in product((-1,1),repeat=2)}
def column(k,ds):
    return tuple((comb(d,k)-comb(ref,k))%9 for ref,*pair in ds for d in pair)
def solve(ds):
    state={(0,0,0,0):(0,())}
    for k in range(1,29):
        col=column(k,ds)
        nxt={}
        for v,(cost,coeff) in state.items():
            for u in range(-4,5):
                w=tuple((x+u*y)%9 for x,y in zip(v,col))
                cand=(cost+abs(u),coeff+(u,))
                if w not in nxt or cand<nxt[w]:
                    nxt[w]=cand
        state=nxt
        if targets & state.keys():
            matches=[(state[v],v) for v in targets if v in state]
            return k,len(state),min(matches)
for ds in (chords,tuple(tuple(n//2 for n in row) for row in chords)):
    assert all(row[0]%3==row[1]%3 for row in ds)
    result=solve(ds)
    print("squared chord convention",ds,"first integer-valued polynomial degree, lattice size, witness",result)

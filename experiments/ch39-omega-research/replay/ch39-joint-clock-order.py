"""Coherent bounds on the neutral and signed spin-tower energies at the fitted zero."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
neutral = runpy.run_path(str(root/"ch39-neutral-clock-consecutive.py"))
q,r,raw = (neutral[k] for k in ("q","r","raw"))
lookup = {(i,j):d for i,j,d,_,_ in raw}
from flint import arb_mat
I = arb_mat([[int(a==b) for a in range(3)] for b in range(3)])
def common(n):
    losses = []
    R = neutral["rotations"](n)[0]
    S = neutral["rotations"](719)[0]
    for i,j,k in combinations(range(19),3):
        pairs = ((i,j),(j,k),(i,k))
        if not all(p in lookup for p in pairs):
            continue
        def loss(turns):
            a,b,c = [turns.get(lookup[p],I) for p in pairs]
            return 2*(3-(a*b*c.transpose()).trace())
        lower = min(loss(R).lower(),loss(S).lower())
        if lower>0:
            losses.append((pairs,lower))
    return losses

def edge_bounds(z,upper):
    diag = [z for _ in range(19)]
    weights = {}
    for i,j,d,ell,w in raw:
        if upper:
            y = (z.sqrt()*ell.lower()).lower()
            h = (y/y.tanh()+w.lower()).lower()
            b = (y/y.sinh()+w.upper()).upper()
        else:
            y = (z.sqrt()*ell.upper()).upper()
            h = (y/y.tanh()+w.upper()).upper()
            b = (y/y.sinh()+w.lower()).lower()
        diag[i]+=h
        diag[j]+=h
        weights[i,j]=b
    return {p:b/(diag[p[0]]*diag[p[1]]).sqrt() for p,b in weights.items()}

def correction(n):
    losses = common(n)
    total=q(0)
    lo=q(".25")
    ratio=q(2)**(q(1)/32)
    for _ in range(32*32):
        hi = lo*ratio
        W=edge_bounds(lo,True)
        w=edge_bounds(hi,False)
        val=sum(k*(W[a]*W[b]*W[c]-w[a]*w[b]*w[c]) for (a,b,c),k in losses)
        total += (hi*hi-lo*lo)*val/2
        lo=hi
    return (total/(32*arb.pi()**2*r**4)).lower()

def check():
    # These outward-rounded inputs were independently certified by the archive's
    # consecutive-neutral and signed-neighbor scripts on the same full source box.
    for n,N,F in ((718,"-.66664e-12","2.49872e-12"),(720,"-1.69535e-12","1.32837e-12")):
        delta=correction(n)
        upper = 3*(q(N)-delta)+q(F)
        print(n,"common triangle correction",delta.str(12),"joint upper bound",upper.upper().str(12))
        assert upper<0

if __name__=="__main__":
    check()

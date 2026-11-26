"""Independent original arcs and exact density of an arithmetic SU(2) torus sampler."""
from functools import lru_cache
from pathlib import Path
import runpy

from flint import arb, acb, fmpz

root = Path(__file__).resolve().parent
geo = runpy.run_path(str(root/"ch39-arc-relative-determinant.py"))
theta = geo["theta"]
pi = arb.pi()
special = (12,200,244,724,804)
edge_set = {d for _,_,d,_ in geo["edges"]}
assert all(d in edge_set for d in special)
characters = []
for d in special:
    odd = {int(p) for p,k in fmpz(d*(3468-d)).factor() if k%2}
    characters.append(odd|{-1})
for mask in range(1,1<<len(special)):
    odd = set()
    for i,z in enumerate(characters):
        if mask & (1<<i):
            odd ^= z
    assert odd
    # No squared root-of-unity eigenvalue is needed for this test.
assert all(2*(1734-d)%1734 != 0 for d in special)

def partition(time, cutoff=32):
    t = arb(time)
    terms = [(-t*n*(n+2)).exp() for n in range(cutoff+2)]
    odd = sum((-1)**j*s for j,s in enumerate(terms[:-1]))
    z = odd+arb(0,terms[-1].upper())
    base = sum((j+1)**2*s for j,s in enumerate(terms[:-1]))
    ratio = (arb(cutoff+3)/(cutoff+2))**2*(-t*(2*cutoff+5)).exp()
    assert ratio<1
    tail = ((cutoff+2)**2*terms[-1]/(1-ratio)).upper()
    peak = base+arb(0).union(tail)
    return z,peak,z/(2*pi*peak)

def first_hits(n=128):
    hits=[]
    for i in range(1,n+1):
        h = acb.zeta_zero(i).imag
        p,q = geo["a"]*h,geo["b"]*h
        third = theta(12)*h
        fourth = theta(200)*h/(2*pi)
        phase = fourth-fourth.floor()
        assert 0<phase and phase<1
        gap = p.sin()**2*q.sin()**2*abs(third.sin())-phase
        assert gap>0 or gap<0
        if gap>0:
            hits.append(i)
    return hits

if __name__=="__main__":
    print("Haar maximum positive-zero density",(1/(2*pi)).str(15))
    for t in (".25","1","2"):
        z,peak,rho=partition(t)
        print("odd-sector heat time, partition, maximal density",t,z.str(13),rho.str(13))
    hits=first_hits()
    print("Haar hits through positive zero 128",len(hits),hits)

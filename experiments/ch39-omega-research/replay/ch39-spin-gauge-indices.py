"""Weak-color projection for mass-invariant heavy representations of the two spin towers."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q,base = short["q"],short["base"]
model = base["model"]
G = [1/x**2 for x in model["gauge_inputs"]]
beta = model["gauge"]["beta"]
time = base["T0"].union(base["T1"])
target = 8*arb.pi()**2*(G[1]-G[2])-(beta[1]-beta[2])*time
wl = base["wl"]
# Uniform bounds for the eighteen actual right threshold logarithms.
mu = short["shiftbox"][0]
_,norm = base["arc_budget"]()
emax = (mu+norm).upper()
eta = q("1.018")**2
k = arb.pi()*(emax*eta/base["data"]["d"]).sqrt()
slope = 18*k*(k/2).tan()/emax
lower = (1-10/mu.lower())/((1+slope)*(1+base["data"]["beta"]*eta*(emax+10)))
upper = 1+norm/mu.lower()
shift = (-upper.log()/2).lower().union((-lower.log()/2).upper())
wr = 18*(base["s0"].union(base["s1"])+shift)

def indices(dim):
    # Real or realified SU(2) irreps, including the singlet for padding.
    blocks = [(1,0,0),(6,0,1),(8,0,3)]
    for j2 in range(1,dim):
        copies = 1 if j2%2==0 else 2
        size = copies*(j2+1)
        if size<=dim:
            index = F(size*j2*(j2+2),12)
            assert index.denominator==1
            blocks.append((size,int(index),0))
    sets = [{(0,0)}]+[set() for _ in range(dim)]
    for size,w,c in blocks:
        for n in range(size,dim+1):
            sets[n].update((a+w,b+c) for a,b in sets[n-size])
    return sorted(sets[dim])

def possible():
    rows = []
    for il,cl in indices(8):
        for ir,cr in indices(10):
            l,r = q(F(2*(il-cl),3)),q(F(ir-cr,6))
            if (l*wl+r*wr-target).contains(0):
                rows.append((il,cl,ir,cr))
    return rows

if __name__=="__main__":
    rows = possible()
    print("weak-color target",target.lower().str(10),target.upper().str(10))
    print("log exposures",wl.lower().str(10),wl.upper().str(10),wr.lower().str(10),wr.upper().str(10))
    print("necessary index patterns",rows)

"""Impose the Z6 hypercharge lattice on both physical spin spectra, with neutral Proca fields."""
from fractions import Fraction as F
from pathlib import Path
import runpy

m = runpy.run_path(str(Path(__file__).with_name("ch39-spin-gauge-indices.py")))
q,wl,wr,T,G,beta = (m[k] for k in ("q","wl","wr","time","G","beta"))
from flint import arb

def blocks(dim,charged_color):
    out = set()
    for cdim,index,triality in ((1,0,0),(3,F(1,2),1),(3,F(1,2),2),(8,3,0)):
        if not charged_color and cdim!=1:
            continue
        for twice in range(dim):
            weak_dim = twice+1
            for h in range(51):
                if (h+3*twice+2*triality)%6:
                    continue
                copies = 2 if h or triality or twice%2 else 1
                size = weak_dim*cdim*copies
                if size<=dim and size*h*h<=2500:
                    out.add((size,int(F(size*twice*(twice+2),12)),int(weak_dim*copies*index),size*h*h))
    return sorted(out)

def options(dim,max_charge,color):
    states = [{(0,0,0)}]+[set() for _ in range(dim)]
    for size,t,c,h in blocks(dim,color):
        for k in range(size,dim+1):
            states[k].update((w+t,u+c,z+h) for w,u,z in states[k-size] if z+h<=max_charge)
    return states[-1]

def crossings():
    states = []
    patterns = m["possible"]()
    R12 = 8*arb.pi()**2*(G[0]-G[1])-(beta[0]-beta[1])*T
    max_total = max(((q(F(2*il,3))*wl+q(F(ir,6))*wr+R12).upper()) for il,cl,ir,cr in patterns)
    assert max_total*90/wl.lower()<2500 and max_total*360/wr.lower()<2500
    left,right = options(8,2500,True),options(10,2500,False)
    for il,cl,ir,cr in patterns:
        lw = q(F(2*(il-cl),3))
        rw = q(F(ir-cr,6))
        fixed_wr = ((m["target"]-lw*wl)/rw).intersection(wr)
        for _,_,hl in sorted(v for v in left if v[:2]==(il,cl)):
            for _,_,hr in sorted(v for v in right if v[:2]==(ir,cr)):
                u = q(F(hl,90)-F(2*il,3))*wl + q(F(hr,360)-F(ir,6))*fixed_wr
                if (u-R12).contains(0):
                    states.append((il,cl,hl,ir,cr,hr))
    return states

if __name__=="__main__":
    assert not crossings()
    print("no neutral-vector one-loop assignment")

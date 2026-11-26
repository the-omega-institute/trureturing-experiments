"""Charged gauge-vector threshold, hypercharge lattice, and perturbativity at the first right mass."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-spin-gauge-hypercharge.py")))
b = m["m"]
q,wl,wr,T,G,beta = (m[k] for k in ("q","wl","wr","T","G","beta"))
pi = arb.pi()
scale = 8*pi**2
# At the crossing the strong coupling bounds every perturbative hypercharge.
max_inv_alpha = 4*pi*(G[2]-(beta[2]*T+q(8)/3)/scale)
assert max_inv_alpha<41 and q(50)**2>60*max_inv_alpha
# Just before the earliest right threshold the charged vectors are absent.
peak = b["base"]["s1"]+b["shift"].upper()
max_l = 90*(scale*G[0]-beta[0]*(T-peak))/(wl-18*peak)
assert max_l < 520

def assignment(dim,max_index,allowed):
    # (weak index, color index, squared-integer-charge sum) -> smallest maximal |6Y|.
    states = [dict() for _ in range(dim+1)]
    states[0][0,0,0] = 0
    for col,tc,trial in ((1,0,0),(3,F(1,2),1),(3,F(1,2),2),(8,3,0)):
        for j in range(dim):
            for h in range(50):
                if (h+3*j+2*trial)%6:
                    continue
                copies = 2 if h or trial or j%2 else 1
                size = copies*col*(j+1)
                tw = F(size*j*(j+2),12)
                color = tc*copies*(j+1)
                if size>dim or tw.denominator!=1 or F(color).denominator!=1:
                    continue
                wi,ci = int(tw),int(color)
                for k in range(size,dim+1):
                    for (sw,sc,sh),hmax in states[k-size].items():
                        key = (sw+wi,sc+ci,sh+size*h*h)
                        if key[2]<=max_index and any(key[0]<=a and key[1]<=c for a,c in allowed):
                            val = max(h,hmax)
                            states[k][key] = min(val,states[k].get(key,val))
    data = defaultdict(list)
    for (i,c,h),hmax in states[dim].items():
        if (i,c) in allowed:
            data[i,c].append((h,hmax))
    return data

def search():
    patterns = b["possible"]()
    L = assignment(8,520,set((x,y) for x,y,_,_ in patterns))
    R = assignment(10,10*49**2,set((x,y) for _,_,x,y in patterns))
    dt,ct = beta[1]-beta[2],beta[0]-beta[1]
    weak = scale*(G[1]-G[2])
    hyper = scale*(G[0]-G[1])
    survivors = []
    for il,cl,ir,cr in patterns:
        lw,rw = q(F(2*(il-cl),3)),q(F(ir-cr,6))
        need_w = ((b["target"]-lw*wl)/rw).intersection(wr)
        max_peak = need_w.upper()/18+b["shift"].upper()-b["shift"].lower()
        t_on = T.lower()-max_peak
        for lh,hlmax in L[il,cl]:
            l1 = q(F(lh,90))
            G_pre = G[0]-(beta[0]*t_on+l1*(wl.lower()-18*max_peak))/scale
            h2_bound = 240*pi*G_pre.upper()
            if h2_bound<hlmax**2 or G_pre<1/(4*pi):
                continue
            L12 = l1-q(F(2*il,3))
            v = weak-lw*wl
            u = hyper-L12*wl
            for rh,hrmax in R[ir,cr]:
                if h2_bound<hrmax**2:
                    continue
                for Y in range(1,9):
                    if h2_bound<36*Y**2:
                        continue
                    # Complex vector with gyromagnetic ratio 2, including its longitudinal mode.
                    r1 = q(F(rh-1512*Y*Y,360))
                    R12 = r1-q(F(ir,6))
                    determinant = dt*R12-ct*rw
                    if determinant.contains(0):
                        assert not (u*rw-v*R12).contains(0)
                        continue
                    t,z = (v*R12-u*rw)/determinant,(u*dt-v*ct)/determinant
                    if t.overlaps(T) and z.overlaps(need_w):
                        survivors.append((Y,il,cl,lh,ir,cr,rh,hlmax,hrmax))
    return survivors

if __name__=="__main__":
    rows = sorted(search())
    assert rows == [(2,0,0,144,3,0,4716,6,30),(2,0,0,216,3,0,4068,6,27),
                    (2,0,0,360,3,0,2916,12,24)]
    print("perturbative charged-vector index patterns",rows)

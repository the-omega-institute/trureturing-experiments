"""A weighted Kirchhoff SU(2) graph: positive KK gap, one-loop threshold count, and harmonic one-forms."""
from pathlib import Path
import runpy
from flint import acb_mat, arb

old = runpy.run_path(str(Path(__file__).with_name("ch39-all-energy-weak-tower.py")))
short,q,pi,k = (old[x] for x in ("short","q","pi","k"))
inertia = runpy.run_path(str(Path(__file__).with_name("ch39-common-source-tower.py")))["inertia"]
raw = [(i,j,short["r"]*(old["sh"]["rho"]*old["sh"]["theta"](d)).sqrt()
        *short["fsbox"][short["band"](d)]) for i,j,d,_ in old["sh"]["edges"]]
speed = q(2)
def count_below(x,large_lengths):
    # Below the first edge pole, eliminating continuous edge interiors is exact.
    K=acb_mat(19,19)
    for i,j,L in raw:
        y=x*(L.upper() if large_lengths else L.lower())/speed
        assert y<pi
        sine,cosine = y.sin_cos()
        diag=y*cosine/sine
        K[i,i]+=diag
        K[j,j]+=diag
        K[i,j]=K[j,i]=-y/sine
    return inertia(K,q(0))

def gap_bound():
    lo,hi=q(1),q("1.5")
    assert count_below(lo,True)==1 and count_below(hi,False)>1
    ends=[]
    for large_lengths in (True,False):
        a,b=lo,hi
        for _ in range(28):
            mid=(a+b)/2
            if count_below(mid,large_lengths)==1:
                a=mid
            else:
                b=mid
        ends.append(a if large_lengths else b)
    return ends[0].lower().union(ends[1].upper())

def running():
    # Positive graph spectrum interlaces the Dirichlet spectrum with the full
    # graph's unique constant mode subtracted. Higher energy counts need no inertia.
    intervals=[n*pi/L for _,_,L in raw for n in range(1,10)]
    assert min((10*pi/L).lower() for _,_,L in raw)>8
    loss=q(0)
    for j in range(7*64):
        lo=1+q(j)/64
        hi=lo+q(1)/64
        dr=sum(not (m>hi) for m in intervals)
        kk=max(0,sum(speed*m<lo for m in intervals)-1)
        beta=q("18.5")+dr-4*kk
        loss+=max(q(0),beta.upper())*(hi/lo).log()/(2*pi)
    assert q("562.5")-k*8/2<0
    floor=q("40.6142")-loss
    assert floor>2
    return floor

if __name__=="__main__":
    gap=gap_bound()
    floor=running()
    betti=len(raw)-19+1
    assert betti==117
    shift=q(betti)/3*q("26.8220")/(2*pi)
    print("first positive mass / Q",gap.lower().str(14),gap.upper().str(14))
    print("positive-mode model: minimum inverse alpha2",floor.lower().str(14))
    print("free cycle scalars: loss of inverse alpha2 by Q",shift.lower().str(14))

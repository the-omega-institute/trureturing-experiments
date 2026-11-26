"""All perturbative scales, identical left/right spectra at r=400, two-loop gauge running."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
src = runpy.run_path(str(root/"ch39-common-height-gauge.py"))
g,old,q = src["inputs"],src["old"],src["q"]
D = old["D"]
cap = 4*arb.pi()
proj = src["data"]["g"]
dot = proj["dot"]
nonneg = proj["nonneg"]
P = proj["P"]
B = g["B"]
heavy = [[x+y for x,y in zip(left,right)] for left,right in zip(g["BD"],old["BR"])]
charged = [x+y for x,y in zip(g["d"],map(q,("8/15",0,"1/3")))]
assert all(v.overlaps(q("2/3")) for v in charged)
base = [dot(p,g["G0"])*D/2 for p in P]
projected = [[proj["proj"](p,m) for m in (B,g["D"],heavy,g["BF"])] for p in P]
tr = [nonneg(3*cap)]*3
# Contains both the bare graph and the graph with the arithmetic-mean source.
total_log,width = q(95),q(6)
assert src["left_sum"]<total_log and max(src["exposures"])<width

def crossing_caps():
    two,four = q(2),q(4)
    trial = [cap,two,cap]
    assert g["beta"][1]+dot(B[1],trial)/D < 0
    loss = (charged[1]+dot(heavy[1],trial)/D)*total_log+(8+dot(g["BF"][1],trial)/D)/3
    assert g["G0"][1]-2*loss/D > 1/two
    assert g["beta"][0]-dot(g["D"][0],tr).upper()/D > 0
    assert g["G0"][2] > 1/q("1.4")
    for n in range(19):
        for nf in (q(0),q("1/3")):
            row = [a+n*b+nf*c for a,b,c in zip(B[2],heavy[2],g["BF"][2])]
            assert row[0]>0 and row[1]>0
            leading = g["beta"][2]+charged[2]*n+8*nf
            upper = leading+dot(row,[two,two,q("1.4")])/D
            lower = leading+(4*row[2]-dot(g["D"][2],tr).upper())/D
            assert upper<0 or lower>0
    return [nonneg(v) for v in (two,two,four)]

def correction(row,gs,traces,el,ef,t):
    sm,y,h,f = row
    return ((dot(sm,gs)-dot(y,traces))*t + dot(h,gs)*el + dot(f,gs)*ef)/D

def time_window():
    floor = g["beta"][0]-dot(g["D"][0],tr).upper()/D
    assert floor>0
    t = nonneg((g["G0"][0]*D/(2*floor)).upper())
    for _ in range(5):
        corr = correction(projected[0],crossing_caps(),tr,nonneg(total_log),nonneg(q("1/3")),t)
        t = proj["cut"](t,base[0]-corr)
    return t

def bound(t,steps=64):
    # At fraction f of any trajectory, EL <= both the full exposure and
    # 18 times the maximum elapsed time above its earliest active threshold.
    gs = crossing_caps()
    for _ in range(4):
        hull = []
        for i in range(steps):
            frac = q(i).union(q(i+1))/steps
            elapsed, rel = frac*t,(frac-1)*t
            el = nonneg(min(total_log,(18*g["positive"](width+rel)).upper()))
            ef = nonneg(g["positive"](1+rel).upper()/3)
            box = []
            for j in range(3):
                one = g["beta"][j]*elapsed+charged[j]*el+8*ef
                tw = correction((B[j],g["D"][j],heavy[j],g["BF"][j]),gs,tr,el,ef,elapsed)
                inv = g["G0"][j]-2*(one+tw)/D
                low = max(1/gs[j].upper(),inv.lower())
                box.append(1/low.union(inv.upper()))
            hull = [a.union(b) for a,b in zip(hull,box)] if hull else box
        gs = [proj["cut"](a,b) for a,b in zip(gs,hull)]
    residual = base[1]-correction(projected[1],gs,tr,nonneg(total_log),nonneg(q("1/3")),t)
    return residual.lower(),gs

if __name__=="__main__":
    T = time_window()
    cuts = [T.lower()+(T.upper()-T.lower())*j/16 for j in range(17)]
    checked = [bound(a.union(b))[0] for a,b in zip(cuts,cuts[1:])]
    print("candidate time window",T.lower().str(10),T.upper().str(10))
    print("second projection lower",min(checked).str(12))
    assert min(checked)>0

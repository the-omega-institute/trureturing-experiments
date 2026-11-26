"""Upper energy limit of the perturbative product theory above the common shell crossing."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
joint = runpy.run_path(str(root/"ch39-joint-higgs-threshold.py"))
m,b,q = joint["fixed"],joint["base"],joint["q"]
pi = arb.pi()
D = 16*pi**2
g = b["model"]["inputs"]
rmin = runpy.run_path(str(root/"ch39-shell-scale-bounds.py"))["minimum"]
alpha = joint["spectrum"]()[-1].upper()
shift = (9*pi**2/b["data"]["d"]).sqrt()
assert shift<q(".021")
right_high = q("1.502") + m["norm"]/rmin**2
right_low = (-2*m["v0"]).exp()*(400/rmin)**2 + m["norm"]/rmin**2
left_static = [x.upper()/rmin**2 for x in m["left_static"]]
right_static = [right_low]*18+[right_high]*20
cap = 4*pi

def exposure(endpoint,right):
    eta = m["eta"] if right else q(1)
    dilate = 1/(1-shift*eta.sqrt())
    def interval_log(mass):
        cut = max(q(1),mass.upper())
        return max(q(0),(endpoint/cut).log().lower())
    total = sum(interval_log(dilate*x.sqrt()) for x in (right_static if right else left_static))
    for n,d in b["data"]["factors"]:
        for k in range(1,6):
            mass = k*(d/eta).sqrt() / (rmin*(1-k*k*b["h"]**2))
            total += n*interval_log(dilate*mass)
    return total

def upper(endpoint):
    t = endpoint.log()
    left,right = exposure(endpoint,False),exposure(endpoint,True)
    # All represented heavy modes have nonnegative two-loop entries.
    # The color row of the SM contains the only negative gauge coefficient.
    bounds=[]
    for i,ri in enumerate(map(q,("8/15",0,"1/3"))):
        correction = sum(min(q(0),v) for v in g["B"][i])*cap-3*cap*sum(g["D"][i])
        one = (g["beta"][i]+correction/D)*t+g["d"][i]*left+ri*right+8*t/3
        bounds.append(alpha-one/(2*pi))
    return left,right,bounds

if __name__=="__main__":
    r = q("3.85")
    L,R,A = upper(r)
    print("left/right exposure lower",L.str(12),R.str(12))
    print("inverse gauge upper",*(v.upper().str(12) for v in A))
    assert min(v.upper() for v in A)<1

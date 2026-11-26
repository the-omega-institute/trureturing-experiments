"""A crossing using the flavor clock square and the 1 TeV MS-bar input intervals."""
from pathlib import Path
import runpy
from types import FunctionType

from flint import arb

root = Path(__file__).resolve().parent
spectrum = runpy.run_path(str(root / "ch39-common-height-scalar-levels.py"))
# The archived integrator initializes its reference example on import.
# Recreate that auxiliary table if this command is the first replayed.
if not (root / "ch39-quaternion-low-masses.txt").is_file():
    runpy.run_path(str(root / "ch39-quaternion-scalar-levels.py"), run_name="__main__")
old = runpy.run_path(str(root / "ch39-right-scalar-rge.py"))
data, gauge, q, r = (old[k] for k in ("data", "g", "q", "r"))
low = spectrum["window_levels"]()
beta, ret = data["beta"], data["ret"]
low = [(x.lower() / (1 + beta*(x.lower()+10))).lower().union((x.upper()/ret).upper())
       for x in low]
exposures = [r.log() - x.log()/2 for x in low]
left_sum = sum(exposures)
assert min(exposures)>0 and data["edge19"] > r*r
# Table 2, MS-bar at 1 TeV. All 1-sigma marginal input boxes are enclosed.
yukawa_errors = (".00000014", ".00005", ".0043", ".0000002", ".000003", ".00009",
                 ".0000000006", ".00000013", ".000002")
gauge_inputs = [arb(v,e) for v,e in zip((".467453",".63811","1.0583"),
                                         (".000027",".00003",".0030"))]
yukawa_inputs = [arb(v,e) for v,e in zip(gauge["yc"],yukawa_errors)]
inputs = dict(gauge, G0=[1/x**2 for x in gauge_inputs], y0=[2*x.log() for x in yukawa_inputs])
def rates(gs, y, n):
    assert all(0 < z < 4*arb.pi() for z in gs)
    assert all(z.exp() < 4*arb.pi() for z in y)
    return gauge["rates"](gs,y,n)

inputs["rates"] = rates
fn = old["run"]
integrate = FunctionType(fn.__code__, dict(fn.__globals__, Ls=exposures, WL=left_sum, g=inputs),
                         argdefs=fn.__defaults__)

def moments(s):
    y19 = data["edge19"]/r**2
    count = sum(n for n,factor in data["factors"] if factor/(1-data["h"]**2)**2/r**2 < 2)
    assert 1 < y19 < 2 and count > 18
    y = (-2*s).exp()
    light = 18*(-y).exp()
    rest = sum((-(x/r**2)).exp() for x in low) + (count-18)*(-q(2)).exp() - light
    excess = -18 + (count-18)*(y19-1)*(-y19).exp() + (1-y)*light
    assert rest > 0 and excess > 1/q(1).exp()
    return rest.lower(), excess.lower()

def faces(t0,t1,s0,s1):
    t0,t1,s0,s1 = map(q,(t0,t1,s0,s1))
    T, S = t0.union(t1), s0.union(s1)
    assert min(t0-v.upper() for v in exposures)>0 and t0>s1 and s0>0
    results = [integrate(t,s) for t,s in ((t0,S),(t1,S),(T,s0),(T,s1),(T,S))]
    signs = [results[0][0],results[1][0],results[2][1],results[3][1]]
    for label,x in zip(("T low","T high","s low","s high","inverse alpha"),signs+[results[-1][2]]):
        print(label,x.lower().str(12),x.upper().str(12))
    print("left exposure",left_sum.lower().str(10),left_sum.upper().str(10))
    print("GeV", (1000*T.exp()).lower().str(10), (1000*T.exp()).upper().str(10))
    print("right mass / Q",(-S).exp().lower().str(10),(-S).exp().upper().str(10))
    print("heat lower bounds",*(x.str(10) for x in moments(S)))
    assert signs[0]>0 and signs[1]<0 and signs[2]>0 and signs[3]<0
    assert results[-1][2]>1

if __name__ == "__main__":
    faces("26.88","27.008","2.50","2.63")

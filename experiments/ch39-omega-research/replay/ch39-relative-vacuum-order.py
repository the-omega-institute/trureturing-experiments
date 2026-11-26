"""First surviving UV coefficient with the full, non-diagonal endpoint stiffness."""
from itertools import product
from pathlib import Path
import runpy

m = runpy.run_path(str(Path(__file__).with_name("ch39-relative-vacuum.py")))
q,w,r = m["q"],m["w"],m["r"]
ref = m["potential"][2]
# For an added edge with stiffness a, A_ij A_ji = a^2 I even with holonomy.
square = [q(0) for _ in m["V"]]
for (i,j), idx in m["source_ids"].items():
    weight = m["base"]["weights"][idx]
    square[i] += weight**2
    square[j] += weight**2
centered = [v+a-ref for v,a in zip(m["V"],m["added"])]
mean_square = [sum(centered[i]*centered[i]+square[i] for i in g)/len(g)
               for g in m["groups"]]

def coefficient(u,v,s):
    p = m["branch"](u,v,s)
    total = q(0)
    for i,n in enumerate(m["count"]):
        d = m["length"][i]
        e = m["du"][i]*(u-1)+m["dv"][i]*(v-1)
        f = d+e
        diff3 = e*(f*f+f*d+d*d)
        diff5 = e*(f**4+f**3*d+f*f*d*d+f*d**3+d**4)
        a = m["mean_A"][i]-ref
        # C3=0 replaces Delta tr(D A) by Delta tr(D^3)/3.
        total += n*(e*mean_square[i]+f*p[i]*(2*a+p[i])-diff3*(a+ref/3)-f**3*p[i]+diff5/5)
    return 2*total/r**5

def pieces(center,radius,n):
    return [center+radius*(2*q(k).union(q(k+1))/n-1) for k in range(n)]

def check():
    (u,v),(du,dv) = w["center"],w["radii"]
    s0,s1 = w["s0"],w["s1"]
    box = [pieces(u,du,8), pieces(v,dv,8), pieces((s0+s1)/2,(s1-s0)/2,8)]
    values = [coefficient(*point) for point in product(*box)]
    lo = min(x.lower() for x in values)
    hi = max(x.upper() for x in values)
    print("C5/r^5",lo.str(12),hi.str(12))
    assert lo > 0

if __name__ == "__main__":
    check()

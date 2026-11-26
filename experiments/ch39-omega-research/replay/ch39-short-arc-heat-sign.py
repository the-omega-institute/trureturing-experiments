"""The zero relative volume involves opposite signs at intermediate proper time."""
from pathlib import Path
import runpy

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
order = runpy.run_path(str(Path(__file__).with_name("ch39-relative-vacuum-order.py")))
m,w,q,r = (short[k] for k in ("m","w","q","r"))
ref = order["ref"]

def c5(box,s):
    x = short["coords"](box,s)
    inc = short["increments"](x)
    p = short["branch"](inc,s)
    val = q(0)
    for n,d,e,pj,A,B in zip(m["count"],m["length"],inc,p,m["mean_A"],order["mean_square"]):
        f = d+e
        centered = A-ref
        diff3 = e*(f*f+f*d+d*d)
        diff5 = e*(f**4+f**3*d+f*f*d*d+f*d**3+d**4)
        val += n*(e*B+f*pj*(2*centered+pj)-diff3*(centered+ref/3)-f**3*pj+diff5/5)
    return 2*val/r**5

def check():
    init = [short["full"]]*3+[w["s0"].union(w["s1"])]
    stack = [(init,0)]
    vals = []
    while stack:
        parts,depth = stack.pop()
        v = c5(parts[:3],parts[3])
        if v>q(".00002") and v<q(".0005"):
            vals.append(v)
            continue
        assert depth<20
        i = depth%4
        ax = parts[i]
        mid = ax.mid()
        for seg in (ax.lower().union(mid),mid.union(ax.upper())):
            child = parts[:]
            child[i] = seg
            stack.append((child,depth+1))
    print("C5/r^5 bounds",min(v.lower() for v in vals).str(10),max(v.upper() for v in vals).str(10))
    # The right-to-left comparison at large proper time only needs a spectral gap.
    base = m["base"]
    assert min(x.upper() for x in base["left"]) < 200
    mu = short["shiftbox"][0]
    # Eliminate the finite-dimensional endpoints at E=200, strictly below the first arc pole.
    from flint import arb
    k = arb.pi()*(200*q("1.018")**2/base["data"]["d"]).sqrt()
    assert k<arb.pi() and mu-10-200-18*k*(k/2).tan()>0
    print("left spectral bottom < 200 M^2 < right spectral bottom")

if __name__=="__main__":
    check()

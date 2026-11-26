"""Uniform error bounds for a radical approximation to the nonsolvable physical mass."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
q = arb


def calculate(index):
    inv = [4*q.pi()/s for s in flow["start"]]
    b = flow["beta"]
    slope = arb_mat([[b[0]-b[1],-q(1)/15],[b[1]-b[2],q(1)/6]])
    t,E = slope.solve(2*q.pi()*arb_mat([[inv[0]-inv[1]],[inv[1]-inv[2]]])).entries()
    logC = 106*((q(1000)/125).log()+t)-2*E
    h = q(125)**2
    model = flow["full"]["model"]
    if index==43:
        sing, family, _, _ = model["m"]["transported_spectra"]()
    else:
        sing, family, _ = model["scalar_models"]()
    errors = q("12e-8" if index==43 else "4e-8")
    gaps = [v+arb(0,errors.upper()) for v in flow["full"]["all_levels"](sing,family)]
    a = [1/d for d in gaps]
    mean = sum(a)/len(a)
    variance = sum((v-mean)*(v-mean) for v in a)/len(a)
    radius = max(abs(v-mean).upper() for v in a)
    R = ((logC-sum(d.log() for d in gaps))/53).exp()
    assert variance > 0 and R > 2*radius
    # w=y+mean is between R and R+radius. Taylor's theorem applied to
    # log(w+b_j) and the inequality (w-R)/w <= log(w/R) <= (w-R)/R.
    low = R*variance/(2*(R+2*radius)**2)
    high = (R+radius)*variance/(2*(R-radius)**2)
    # Avoid subtracting overlapping enclosures of the common variance/(2R).
    # Use rational forms to avoid subtracting almost equal quantities.
    r = radius/R
    relative_low = 4*r*(1+r)/(1+2*r)**2
    relative_high = r*(3-r)/(1-r)**2
    last = h*variance/(2*R)*max(relative_low.upper(),relative_high.upper())
    coarse = h*(mean-low.union(high))
    print("window",index,"geometric mean mass² excess GeV²",
          float(coarse.lower()),float(coarse.upper()),
          "mean-corrected mass² shortfall GeV²",float((h*low).lower()),float((h*high).upper()),
          "variance correction absolute error below",last.upper().str(8))
    assert coarse > q("1e4") and h*high < q("1e-12") and last < q("1.75e-29")


if __name__=="__main__":
    for n in (43,195):
        calculate(n)

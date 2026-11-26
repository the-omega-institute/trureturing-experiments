"""Holonomy-preserving one-loop vacuum cancellation and its gauge-threshold costs."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

folder = Path(__file__).resolve().parent
flow = run_path(str(folder/"ch39-shell-gw-first-heavy-flow.py"))
q = arb


def charged_singlets():
    Q = q("4.352e13").union(q("4.473e13"))
    inv = q("38.82").union(q("38.98"))/(4*arb.pi())
    ky = 5*inv/3
    t = (q("2.435e18")/Q).log()/(8*arb.pi()**2)
    for n in (18,20,22):
        y = ky - t*(q(47)/3+n)
        if n == 22:
            assert y < 0
            continue
        critical = (267*t+(16*t)**2/y)/inv
        print("flat-paired singlet scalar count", n, "critical kX between",
              float(critical.lower()), float(critical.upper()))
        assert y > 0 and 267*t+(16*t)**2/y < 100


def paired_shell(n):
    beta = flow["beta"]
    inverse = [4*arb.pi()/s for s in flow["start"]]
    slope = arb_mat([[beta[0]-beta[1], -q(1)/15], [beta[1]-beta[2], q(1)/6]])
    t, required = slope.solve(2*arb.pi()*arb_mat([[inverse[0]-inverse[1]],
                                                  [inverse[1]-inverse[2]]])).entries()
    Q = 1000*t.exp()
    full = flow["full"]
    model = full["model"]
    if n == 43:
        sing, family, _, _ = model["m"]["transported_spectra"]()
    else:
        sing, family, _ = model["scalar_models"]()
    gaps = full["all_levels"](sing,family)
    error = q("12e-8" if n == 43 else "4e-8")
    # Include the perturbation of every eigenvalue; the singlet ground is fixed.
    gaps = [d+arb(0,error.upper()) for d in gaps]
    assert all(d > 0 for d in gaps)
    distinct = all(a < b for a,b in zip(gaps,gaps[1:]))
    print("window", n, "all 53 heavy thresholds separated", distinct)
    # The explicit minimal duplication test requires one original doublet
    # per mass. Positive exposure can also be evaluated for a chosen doubled spectrum.
    assert distinct
    def exposure(M, upper):
        return sum(arb.max(q(0), (Q/(125**2+M**2*(v.upper() if upper else v.lower())).sqrt()).log())
                   for v in gaps)
    # Bisect the two comparison functions independently.
    bounds = []
    for high_gaps in (True,False):
        left,right = q("1e12"),q("3e13")
        assert 6*exposure(left, high_gaps)>required and 6*exposure(right,high_gaps)<required
        for _ in range(30):
            mid = (left+right)/2
            value = 6*exposure(mid,high_gaps)-required
            if value > 0:
                left = mid
            else:
                assert value < 0
                right = mid
        bounds.append(left if high_gaps else right)
    scale = bounds[0].union(bounds[1])
    first, last = [(125**2+scale**2*d).sqrt() for d in (gaps[0],gaps[-1])]
    inv_u = inverse[2]-beta[2]*t/(2*arb.pi())
    pole = Q*(2*arb.pi()*inv_u/(beta[1]+53)).exp()
    print("Q, M, first threshold GeV, top threshold/Q and SU(2) pole GeV",
          [(float(v.lower()),float(v.upper())) for v in (Q,scale,first,last/Q,pole)])
    assert first > q("7.9e11") and last < Q and pole < q("1e16")


if __name__ == "__main__":
    charged_singlets()
    for index in (43,195):
        paired_shell(index)

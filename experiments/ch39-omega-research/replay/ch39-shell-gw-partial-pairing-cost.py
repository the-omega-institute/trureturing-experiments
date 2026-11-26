"""Optimal finite vacuum matching cost over shell levels paired at their original masses."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

old = run_path(str(Path(__file__).with_name("ch39-shell-gw-flat-determinant-matching.py")))
flow, q = old["flow"], old["q"]


def bounds(index):
    beta = flow["beta"]
    inverse = [4*arb.pi()/s for s in flow["start"]]
    mat = arb_mat([[beta[0]-beta[1], -q(1)/15], [beta[1]-beta[2], q(1)/6]])
    t, E = mat.solve(2*arb.pi()*arb_mat([[inverse[0]-inverse[1]],
                                          [inverse[1]-inverse[2]]])).entries()
    Q = 1000*t.exp()
    inv = (inverse[2]-beta[2]*t/(2*arb.pi()))/(4*arb.pi())
    time = (q("2.435e18")/Q).log()/(8*arb.pi()**2)
    cap = (inv/time-q(17)/3)*6/5
    assert 20 < cap < 21
    for charged, expected in ((18,3),(20,1)):
        cap_y = (5*inv/(3*time)-256*time/(100-267*time)-q(47)/3-charged)*6/5
        assert expected < cap_y < expected+1
    full = flow["full"]
    model = full["model"]
    if index == 43:
        sing, family, _, _ = model["m"]["transported_spectra"]()
    else:
        sing, family, _ = model["scalar_models"]()
    w = "12e-8" if index == 43 else "4e-8"
    ds = [d+arb(0,w) for d in full["all_levels"](sing,family)]
    assert all(a < b for a,b in zip(ds,ds[1:]))
    correction = q(125)**2/(q("1e12")**2*ds[0].lower())
    out = []
    for count in list(range(21))+[53]:
        split = 53-count
        weights = [1]*split+[6]*count
        total = sum(weights)
        # M0 is the all-active root at zero additive weak mass.
        M0 = Q*((-E-sum(d.log()*k/2 for d,k in zip(ds,weights)))/total).exp()
        lower = M0/(1+correction).sqrt()
        assert lower > q("1e12")
        scale = lower.lower().union(M0.upper())
        masses2 = [q(125)**2+scale**2*d for d in ds]
        assert masses2[-1] < Q**2
        finite = sum(m*m for m in masses2[:split])
        out.append((finite,masses2[0].sqrt()))
    baseline = out[0][0]
    ratios = [v/baseline for v,_ in out[1:21]]
    assert all(a < b for a,b in zip(ratios,ratios[1:]))
    print("window", index, "smallest cost ratio for 1..20 paired and for exactly 20",
          ratios[0].lower().str(8), ratios[-1].lower().str(8))
    print("smallest first opening at 1,6,7 paired GeV",
          *(out[j][1].lower().str(8) for j in (1,6,7)))
    assert ratios[0]>q("2.53") and ratios[-1]>600
    if index == 43:
        assert out[7][1]>q("1.48e11")
    else:
        assert out[1][1]>q("1.48e11")


if __name__=="__main__":
    for idx in (43,195):
        bounds(idx)

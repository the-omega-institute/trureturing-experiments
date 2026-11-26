"""Necessary tensor deformation and a sufficient weak heavy matching branch."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

flow = run_path(str(Path(__file__).with_name("ch39-shell-gw-first-heavy-flow.py")))
fit = flow["fit"]


def spectral_weight(F,D,w):
    P = flow["acb_mat"](fit["rows"](F)+fit["rows"](D))
    H = P*fit["adj"](P)
    rows = fit["rows"](H)
    values = np.linalg.eigvalsh(fit["common"]["flavor"]["array"](H))[::-1]
    lower, upper = [], []
    for i,t in enumerate(values[:18]):
        b = arb(float(t))-arb(".00001")
        assert b > 0 and fit["count"](rows,b) <= 17-i
        lower.append(b)
        hi = arb(float(t))+arb(".00001")
        assert fit["count"](rows,hi) >= 18-i
        upper.append(hi)
    return sum(v*p for v,p in zip(lower,w)), upper


def check(index):
    heavy,up,down,lepton,F,D,lower,*_ = flow["at_threshold"](index)
    Q,_ = flow["masses"](index)
    loop = flow["loop"]
    weight = [(9*(m/heavy[-1]).log()/loop).exp() for m in [heavy[0],*heavy]]
    weighted, eigen_upper = spectral_weight(F,D,weight)
    total = fit["hs"](F,F).real+fit["hs"](D,D).real
    dt = (Q/heavy[-1]).log()/loop
    pole_cap = 9/(1-(-9*dt).exp())
    partial_cap = 9/(1-(1-9/(48*arb.pi()))*(-9*dt).exp())
    ratios = [weighted.sqrt()/total.sqrt()-(cap/total).sqrt()/lower for cap in (pole_cap,partial_cap)]
    assert ratios[0]>arb(".063") and ratios[1]>arb(".18")
    # Choose an existing equal-site light fit at g=7.9. Give the orthogonal
    # heavy subspace a common extra factor kappa at its own matching thresholds.
    # Ky Fan bounds every cumulative mass-ordered heavy injection by the
    # corresponding prefix sum of the descending eigenvalues (light omitted).
    # Advancing a positive jump earlier increases the solution of T'=10*T².
    kappa = arb(".10")
    cap = sum(fit["hs"](m,m).real for m in (up,down,lepton))
    times = [*heavy[:18],Q]
    for j,eigenvalue in enumerate(eigen_upper):
        cap += kappa**2*arb("7.9")**2*eigenvalue
        dt_upper = arb.max(arb(0),(times[j+1]/times[j]).log().upper())/flow["loop"]
        assert 10*dt_upper*cap < 1
        cap = cap/(1-10*dt_upper*cap)
    assert cap < arb("5.40") < 8*arb.pi()/3
    print("zero",index,"relative threshold changes required to avoid pole or quark block violation",
          *[r.lower().str(8) for r in ratios],"upper total Yukawa norm squared at kappa=.10",cap.upper().str(9))


if __name__=="__main__":
    for n in (43,195):
        check(n)

"""Orientation dependence of the strict two-loop logarithm over all heavy thresholds."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

inert = run_path(str(Path(__file__).with_name("ch39-shell-gw-two-inert-heavy.py")))
flow,fit = inert["flow"],inert["fit"]
adj,rows,count = [fit[k] for k in ("adj","rows","count")]


def frob_upper(M):
    return sum(abs(v).abs_upper()**2 for v in M.entries()).sqrt()


def spectral_log_block(n):
    model = flow["full"]["model"]
    if n == 43:
        sing,family,_,_ = model["m"]["transported_spectra"]()
    else:
        sing,family,_ = model["scalar_models"]()
    light = inert["simple_mode"](sing,0)
    l0 = inert["hs"](light,sing*light).real
    eye = fit["common"]["old"]["m"]["eye"](48)
    mass = arb("1e12").union(arb("1.2e12"))
    M = family+(125**2/mass**2-l0)*eye
    minimum = arb(".002")
    assert count(rows(M),minimum)==0 and count(rows(M),arb(8))==48
    vs, vectors = np.linalg.eigh(inert["array"](M))
    V = inert["exact_array"](vectors)
    lambdas = [arb(float(v)) for v in vs]
    assert all(minimum < v < arb(8) for v in lambdas)
    L = acb_mat(48,48)
    LV = acb_mat(48,48)
    for i,t in enumerate(lambdas):
        L[i,i] = t
        LV[i,i] = t.log()
    # log(M)V - V log(L) has norm <= ||MV-VL||/minimum by the resolvent formula.
    defect = frob_upper(eye-V*adj(V))
    residual = frob_upper(M*V-V*L)
    bound = (residual*(1+defect).sqrt()/minimum+7*defect)/2
    assert bound < arb("1e-7")
    approx = -(V*LV*adj(V))/2
    block = acb_mat([row[16:24] for row in rows(approx)[:8]])
    thresholds = [arb(s) for s in ((".65",".65",".18",".18",".15",".15",".15")
                                    if n==43 else (".50",".50",".18",".18",".15",".15",".15"))]
    gram = block*adj(block)
    for j,t in enumerate(thresholds):
        assert count(rows(gram),(t+bound)**2) <= 7-j
    return thresholds


def check(n):
    # Start from an actual fitted balanced light line; all six orthogonal
    # complements can then be rotated independently.
    inert["certify"](n)
    _,D = fit["readouts"](n)
    threshold = arb(".0374")
    for j in (0,2):
        block = acb_mat([row[9*j:9*j+9] for row in rows(D)])
        assert count(rows(adj(block)*block),threshold**2)==0
    s = spectral_log_block(n)
    spread = 6/(32*arb.pi()**3)*4*inert["g"]**2*threshold**2*sum(s)
    assert spread > arb(".0044" if n==43 else ".0038")
    delta_log = 2*arb.pi()*spread/44
    print("zero",n,"range of integrated two-loop mismatch exceeds",spread.lower().str(9),
          "formal crossing log-scale range exceeds",delta_log.lower().str(9))


if __name__ == "__main__":
    for n in (43,195):
        check(n)

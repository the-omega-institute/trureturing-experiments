"""Five phase-independent singlet thresholds force degree at least six for the physical mass."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb, arb_mat

fixed = run_path(str(Path(__file__).with_name("ch39-shell-gw-fixed-scale-density.py")))
q, flow, full, shell = (fixed[k] for k in ("q","flow","full","shell"))


def certify():
    single, _, _, _ = shell["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    intervals = []
    for i, root in enumerate(np.linalg.eigvalsh(full["model"]["array"](single))):
        lo, hi = q(float(root))-q("1e-9"), q(float(root))+q("1e-9")
        assert count(rows(single),lo)==i and count(rows(single),hi)==i+1
        intervals.append(lo.union(hi))
    gaps = [v-intervals[0] for v in intervals[1:]]
    assert gaps[0] > 0 and all(x < y for x,y in zip(gaps,gaps[1:]))
    excluded = [-q(125)**2/v for v in gaps]+[q(0)]
    assert all(x < y for x,y in zip(excluded,excluded[1:]))
    inverse = [4*q.pi()/s for s in flow["start"]]
    b = flow["beta"]
    mat = arb_mat([[b[0]-b[1],-q(1)/15],[b[1]-b[2],q(1)/6]])
    t,E = mat.solve(2*q.pi()*arb_mat([[inverse[0]-inverse[1]],
                                      [inverse[1]-inverse[2]]])).entries()
    log_ratio = 106*((q(1000)/125).log()+t)-2*E
    assert log_ratio > 2400
    for n in (43,195):
        assert fixed["check"](n)[1] > 12
    print("six real omitted squared-mass values in GeV²",
          [float(v.mid()) for v in excluded])
    print("log of required threshold product relative to its M=0 value",log_ratio.str(10))


if __name__ == "__main__":
    certify()

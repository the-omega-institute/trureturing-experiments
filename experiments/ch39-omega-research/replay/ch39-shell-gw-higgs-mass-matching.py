"""Finite one-loop heavy-doublet mass thresholds of the positive 43rd-zero branch."""
from pathlib import Path
from runpy import run_path
from flint import arb

base=run_path(str(Path(__file__).with_name("ch39-shell-gw-full-inert-potential.py")))
q=base["q"]


def check():
    thresholds=[]
    base["check"](q(81,977)**5,threshold_steps=thresholds)
    assert len(thresholds)==53 and all(p>0 for _,p in thresholds)
    jumps=[p*m**2/(16*arb.pi()**2) for m,p in thresholds]
    electroweak=q(125)**2/2
    print("first and last mass-counterterm magnitudes in GeV²",
          [(float(v.lower()),float(v.upper())) for v in (jumps[0],jumps[-1])])
    print("sum of local threshold magnitudes in GeV²",
          float(sum(jumps).lower()),float(sum(jumps).upper()),"and relative to mH²/2",
          (sum(jumps)/electroweak).lower().str(8))
    print("heavy portal traces, first, third, last",
          [(float(p.lower()),float(p.upper())) for _,p in [thresholds[i] for i in (0,2,-1)]])
    assert jumps[0]>q("4e19") and jumps[-1]>q("2e21") and sum(jumps)>q("4e22")


if __name__=="__main__":
    check()

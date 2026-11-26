"""Compare the fitted common strength with a high-energy color-singlet Yukawa partial wave."""
from pathlib import Path
from runpy import run_path
from flint import arb

folder=Path(__file__).resolve().parent
density=run_path(str(folder/"ch39-shell-gw-positive-density-cp.py"))
common=density["cp"]["common"]
adj,rows,count=(common[k] for k in ("adj","rows","count"))


def check(index,width):
    F=common["old"]["F"] if index==43 else common["operators"]()[1]
    t=common["native"]["relative_phase"](arb(676))
    D=common["conjugate"](F,t)
    G=adj(F)*F
    limit=arb(".01080")
    assert count(rows(G-60*adj(D)*D),limit)==54
    assert count(rows(F*adj(F)),arb(".5134")**2)<9
    e=density["window"](index,width)
    # If a fit used g<=8.20 at a perturbed arithmetic phase, the central
    # readouts on its physical field would obey the following stronger lower bound.
    bound=arb("8.20")
    up,down=(common["old"][name].sqrt() for name in ("u2","d2"))
    assert up>e*bound
    assert (up-e*bound)**2-60*(down+e*bound)**2>limit*bound**2
    gain=arb(".5134")-e
    perturbative=(8*arb.pi()/3).sqrt()/gain
    wave=3*(bound*gain)**2/(16*arb.pi())
    assert perturbative<arb("5.64") and wave>arb("1.05")
    print("zero",index,"fit requires common g > 8.20; gaugeless partial wave >",
          wave.lower().str(8),"tree bound on g <",perturbative.upper().str(8))


if __name__=="__main__":
    for n,w in ((43,"3e-8"),(195,"1e-8")):
        check(n,arb(w))

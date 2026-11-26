"""Retune the two onsets with the first neutral-current gauge logarithm included."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

root=Path(__file__).parent
empirical=run_path(str(root/"ch39-empirical-single-word.py"))
current=run_path(str(root/"ch39-neutral-full-light-current.py"))
logs=("24.31620433001558","26.71271598175277","26.81726493459405")


def verify():
    q=empirical["q"]
    h3,h8,t=map(q,logs)
    word=empirical["endpoint"]
    weights=empirical["events"]["orbit_weights"](46,word)
    reference=weights[empirical["reference"]]
    values=[q(abs(Fraction(a,word[0]))) for a in word[1:]]
    values.extend(q(Fraction(weights[k],reference)).sqrt() for k in empirical["chosen"])
    assert all(abs(v-c)<d for v,c,d in zip(values,empirical["data"]["centers"],
                                           empirical["data"]["widths"],strict=True))
    energy,offsets,_,_=current["light_current"](logs)
    start=t+energy.log()/2
    flow=empirical["adjoint"]["flow"]
    segments=((arb(0),h3,160,flow(False,False)), (h3,h8,28,flow(True,False)),
              (h8,t,12,flow(True,True)), (t,arb(33),64,flow(True,True)))
    empirical["triplet"]["verify_ground"](values,segments,2,t,q("4.8e-8"),q(".136"),
                                         logarithm=(start,h3,offsets))
    changes=[(a-q(b)).exp()-1 for a,b in zip((h3,h8,t),current["neutral"]["primary_logs"])]
    print("relative changes in triplet, octet and trial scale",*(v.str(12) for v in changes))
    print("retuned GeV thresholds, scale, light neutral mass",*((arb(1000)*z.exp()).str(12)
                                                               for z in (h3,h8,t,start)))


if __name__=="__main__":
    verify()

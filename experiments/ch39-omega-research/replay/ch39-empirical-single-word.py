"""Find a single short preparation word for all twelve empirical inputs."""
from fractions import Fraction
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
adjoint = run_path(str(root/"ch39-published-split-adjoint.py"))
triplet = adjoint["triplet"]
data, comp, q = triplet["data"],triplet["comp"],triplet["q"]
single = comp["single"]
events = single["events"]
words = single["source"]["shell"]["words"]
endpoint = (2636825,1232550,1682658,2787494)
reference = (20,10,6)
chosen = ((17,15,14),(15,14,12),(24,10,1),(28,18,0),(38,6,0),
          (35,1,1),(4,0,0),(30,14,0),(35,4,3))


def verify():
    widths = list(map(Fraction,data["width_strings"]))
    bands = [(c-d,c+d) for c,d in zip(words["CENTERS"],widths[:3],strict=True)]
    assert all(not words["hits"](depth,bands) for depth in range(19))
    hits = words["hits"](19,bands)
    assert len(hits) == 6
    low = {k:(events["centers"][k]-d)**2 for k,d in zip(events["names"],widths[3:],strict=True)}
    high = {k:(events["centers"][k]+d)**2 for k,d in zip(events["names"],widths[3:],strict=True)}
    successful = []
    for word, _ in hits:
        for copies in range(1,47):
            weights = events["orbit_weights"](copies,word)
            refs = single["reference_candidates"](weights,low,high)
            if refs:
                successful.append((word,copies,refs))
    assert successful == [(endpoint,46,{reference})]
    weights = events["orbit_weights"](46,endpoint)
    denominator = weights[reference]
    assert len(set((*chosen,reference))) == 10
    ys = [Fraction(weights[k],denominator) for k in chosen]
    assert all(low[k] < y < high[k] for k,y in zip(events["names"],ys,strict=True))
    multiplicities = [sum(low[k]*denominator <= x <= high[k]*denominator for x in weights.values())
                      for k in events["names"]]
    assert all(high[key] < 1 for key in low)
    assert all(high[k] < low[j] or high[j] < low[k] for i,k in enumerate(low) for j in list(low)[i+1:])
    assert multiplicities == [1,5,3,1,4,4,1,1,1] and prod(multiplicities) == 240
    values = [q(abs(Fraction(v,endpoint[0]))) for v in endpoint[1:]]+[q(y).sqrt() for y in ys]
    assert all(abs(v-c) < d for v,c,d in zip(values,data["centers"],data["widths"],strict=True))
    t = q("26.81725683017028")
    h3, h8 = q("24.3145558253697"),q("26.71304007023748")
    flow = adjoint["flow"]
    segments = ((arb(0),h3,160,flow(False,False)),(h3,h8,28,flow(True,False)),
                (h8,t,12,flow(True,True)),(t,arb(33),64,flow(True,True)))
    triplet["verify_ground"](values,segments,2,t,q("9.08e-8"),q(".136"))
    print("single-word empirical source: minimum depth 19, minimum 46 copies at that depth, 240 flavor tables")
    print("twelve signed offsets in printed half-widths",*((v-c)/d for v,c,d in zip(values,data["centers"],data["widths"],strict=True)))
    print("adjoint masses and trial scale in GeV",*((1000*x.exp()).str(11) for x in (h3,h8,t)))


if __name__ == "__main__":
    verify()

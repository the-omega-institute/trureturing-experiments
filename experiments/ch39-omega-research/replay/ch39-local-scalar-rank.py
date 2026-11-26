"""Rank, sign and vertex support required of an endpoint perturbation between the heavy towers."""
from pathlib import Path
import runpy
from flint import arb

root = Path(__file__).resolve().parent
overlap = runpy.run_path(str(root/"ch39-scalar-overlap-bound.py"))
source = overlap["source"]
q = overlap["q"]

def check(low,weights):
    T,W = overlap["ranges"](low)
    H = source["static"](weights)
    # On an unchanged vertex trial space all added endpoint terms vanish.
    # Unchanged edge interiors admit the same affine interpolation in a parallel frame.
    diagonal = [H[2*i,2*i].real for i in range(9)]
    assert max(diagonal)*18<400**2
    assert min(low)>1 and T>q(400).log()
    logs = [(q(400).log()-v.log()/2).lower() for v in diagonal]
    five = 2*sum(sorted(logs)[:5])
    outer = 2*sum(sorted(logs[:6])[:5])
    # Every left/right perturbation with at most eight positive directions
    # retains at least the log exposure of the ten heaviest low modes by min-max.
    ten = sum((q(400).log()-v.log()/2).lower() for v in low[-10:])
    assert min(five,outer,ten)>W
    print("five untouched sites; five outer; positive rank <=8",*(z.str(10) for z in (five,outer,ten)))

if __name__=="__main__":
    check(overlap["profile"]["low"],[arb(0)]*13)
    check(source["mass_bounds"](),[v.upper() for v in source["source_bounds"]()])

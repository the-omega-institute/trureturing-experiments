"""All strictly convex spectral actions at large scalar norm by Ky Fan majorization."""
from pathlib import Path
from runpy import run_path
from flint import acb_mat, arb

m=run_path(str(Path(__file__).with_name("ch39-shell-gw-positive-spectral-moments.py")))
s=m["s"]
rows,count,hs,adj=m["rows"],m["count"],m["hs"],m["adj"]
a,b=s["a"],s["b"]
unit=s["eye"](18)


if __name__=="__main__":
    K=m["kin"]["scalar_kinetic"]()
    assert count(rows(K+600*m["constraint"]),arb("1.119"))==0
    for hop in (a,b):
        dev=hop*hop-4*unit
        assert count(rows(dev),-arb("4.577"))==0 and count(rows(dev),arb("4.577"))==18
    x=acb_mat([[m["p"][i//3]/arb(18).sqrt() if i==j else 0
                for j in range(18)] for i in range(18)])
    difference=a*x-x*b
    assert count(rows(adj(difference)*difference),arb(".05875")**2)==18
    coefficient=(2*arb("1.119")).sqrt()/18-arb(".05875")
    T=(arb("4.577")/coefficient)**2
    assert T<35400
    print("all nontrivial Ky Fan sums strictly exceed the singlet comparator for T above",
          T.upper().str(9))

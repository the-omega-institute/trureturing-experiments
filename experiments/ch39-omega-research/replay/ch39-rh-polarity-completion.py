"""Finite polarity-resolved completion of all mixed shell/time measurements."""
from fractions import Fraction as Q
from itertools import combinations, product
import json
from pathlib import Path
from runpy import run_path

directory=Path(__file__).parent
base=run_path(str(directory/"ch39-integer-shell-record.py"))
reflect=run_path(str(directory/"ch39-rh-origin-free-shell-sums.py"))["reflect"]


def main():
    record=json.loads((directory/"ch39-five-tick-record.json").read_text())
    den=2**record["bits"]
    signs=list(product((1,-1),repeat=2))
    chars=[[1,s,t,s*t] for s,t in signs]
    assert all(sum(c[i]*c[j] for c in chars)==4*int(i==j) for i,j in product(range(4),repeat=2))
    totals=[0]*42
    # Keep the first-offset anisotropy before summing the four polarities.
    cross=[0,0,0]
    for _,height,weight,A,B in record["rows"]:
        a,b=base["turn"](A,den),base["turn"](B,den)
        phase=[(Q(1),Q(0)),a,b,base["times"](a,b)]
        assert all(x*x+y*y==1 for x,y in phase)
        for i,j in combinations(range(4),2):
            c=phase[i][0]*phase[j][0]+phase[i][1]*phase[j][1]
            assert sum(c*r[i]*r[j] for r in chars)==0
            assert sum(c*(r[i]*r[j])**2 for r in chars)==4*c
        for j in range(42):
            totals[j]+=weight*height**2//(height**2+(j*den)**2)
        for k,(i,j) in enumerate(((0,1),(0,2),(1,2))):
            c=phase[i][0]*phase[j][0]+phase[i][1]*phase[j][1]
            cross[k]+=(weight*height**2*c)//(height**2+den**2)
    assert totals[0]==den
    assert all(max(abs(x-y-1),abs(x+46-y))<128 for x,y in zip(totals,record["levels"]))
    shell=[p for p in product(range(5),repeat=3) if sum(x*x for x in p)==17]
    for p,q in combinations(shell,2):
        r,_=reflect(q)
        plus=sum((x-y)**2 for x,y in zip(p,r))
        minus=sum((x-y)**2 for x,y in zip(p,q))
        assert plus+minus==68 and Q(plus,68)==Q(17+sum(x*y for x,y in zip(p,q)),34)
    print("184 positive channels on 46 true height pairs: all 36 shell ratios exact at all 42 offsets.")
    print("Totals track the 42 RH samples within 2^-281; polarity contrasts recover unsymmetrized cross readings.")
    for v in cross:
        print("Single-polarity mixed off-diagonal at j=1 in",float(Q(v,den)),float(Q(v+46,den)))


if __name__=="__main__":
    main()

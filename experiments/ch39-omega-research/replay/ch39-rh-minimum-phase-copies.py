"""Four height-independent phase copies are necessary on the existing finite clock record."""
from fractions import Fraction as Q
from itertools import combinations
import json
from pathlib import Path
from runpy import run_path

directory=Path(__file__).parent
base=run_path(str(directory/"ch39-integer-shell-record.py"))
mul=base["times"]
tol=2**125


def main():
    d=json.loads((directory/"ch39-five-tick-record.json").read_text())
    den=2**d["bits"]
    pairs=list(combinations(range(4),2))
    sums={key:[[0,0],[0,0]] for key in pairs}
    for _,g,weight,A,B in d["rows"]:
        a,b=base["turn"](A,den),base["turn"](B,den)
        phases=[(Q(1),Q(0)),a,b,mul(a,b)]
        for pair in pairs:
            i,k=pair
            rel=mul((phases[i][0],-phases[i][1]),phases[k])
            for j in (1,2):
                t=Q(weight,den)*Q(g*g,g*g+(j*den)**2)
                for c in range(2):
                    sums[pair][j-1][c]+=(tol*t*rel[c])//1
    worst=Q(0)
    for pair in pairs:
        (x,y),(z,w)=[[Q(k,tol) for k in row] for row in sums[pair]]
        error=Q(len(d["rows"]),tol)
        det=abs(x*w-y*z)-error*(abs(x)+abs(y)+abs(z)+abs(w)+2*error)
        assert det>0
        # Real and imaginary parts of a universal complex copy overlap have
        # two real measured equations; this bounds the sum of their magnitudes.
        gain=(abs(x)+abs(y)+abs(z)+abs(w)+4*error)/det
        worst=max(worst,gain)
        print(pair,"inverse gain ceiling",float(gain))
    assert worst<3945
    print("An alphabet of at most three phase copies forces at least one off-diagonal error >",float(1/(3*worst)))


if __name__=="__main__":
    main()

"""Integer-only replay of the entire five-tick joint observation."""
import json
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

directory=Path(__file__).parent
base=run_path(str(directory/"ch39-integer-shell-record.py"))


def main():
    d=json.loads((directory/"ch39-five-tick-record.json").read_text())
    den=2**d["bits"]
    rows=d["rows"]
    assert len(rows)==46 and max(r[0] for r in rows)==1279
    assert min(r[2] for r in rows)>0 and sum(r[2] for r in rows)==den
    phases=[Q(0)]*4
    # All accumulation rounds down onto a fixed integer grid; at most
    # one unit per summand is lost. This bounds numerator sizes uniformly.
    totals=[0]*42
    for _,g,count,two,three in rows:
        a,b=base["turn"](two,den),base["turn"](three,den)
        assert all(x*x+y*y==1 for x,y in (a,b))
        cs=[a[0],b[0],base["times"]((a[0],-a[1]),b)[0],base["times"](a,b)[0]]
        phases=[old+Q(count,den)*v for old,v in zip(phases,cs)]
        for j in range(42):
            totals[j]+=count*g*g//(g*g+(j*den)**2)
    assert totals[0]==den and len(d["levels"])==42
    assert all(max(abs(total-target-1),abs(total+46-target))<128
               for total,target in zip(totals,d["levels"]))
    t=max(map(abs,phases))
    assert t<Q(1,2**281) and 8*t/(1-4*t)<Q(1,2**277)
    # The arithmetic discrepancy is below 2^-281; for the normalized Gram
    # use direct positive interval products and avoid square roots.
    f=[(Q(v,den),Q(v+46,den)) for v in totals]
    for j in range(1,4):
        for i in range(j):
            lo,hi=f[i]; left,right=f[j]
            raw_lo=(i*lo+j*left)/(i+j)
            raw_hi=(i*hi+j*right)/(i+j)
            k=d["gram_bins"][j][i]
            assert Q(k,2**96)**2<raw_lo**2/(hi*right)
            assert raw_hi**2/(lo*left)<Q(k+1,2**96)**2
    print("All 42 squared overlaps agree within 2^-281; the normalized nine-point shell agrees within 2^-277.")
    print("Integer replay retains all six calibrated 96-bit Gram bins.")


if __name__=="__main__":
    main()

"""Classify all 240 empirical flavour tables at the same retuned masses."""
from fractions import Fraction
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb

root=Path(__file__).parent
retuned=run_path(str(root/"ch39-neutral-retuned-ground.py"))
empirical=retuned["empirical"]
q=empirical["q"]


def table_bands():
    weights=empirical["events"]["orbit_weights"](46,empirical["endpoint"])
    reference=weights[empirical["reference"]]
    bands=[]
    for name,width in zip(empirical["events"]["names"],empirical["data"]["width_strings"][3:],strict=True):
        centre=empirical["events"]["centers"][name]
        d=Fraction(width)
        entries=sorted((Fraction(num,reference),orbit) for orbit,num in weights.items()
                       if (centre-d)**2<=Fraction(num,reference)<=(centre+d)**2)
        bands.append(entries)
    return bands


def verify():
    logs=retuned["logs"]
    h3,h8,t=map(q,logs)
    energy,offsets,_,_=retuned["current"]["light_current"](logs)
    start=t+energy.log()/2
    bands=table_bands()
    counts=[len(row) for row in bands]
    assert counts==[1,5,3,1,4,4,1,1,1] and prod(counts)==240
    word=empirical["endpoint"]
    gauge=[q(abs(Fraction(a,word[0]))) for a in word[1:]]
    adj=empirical["adjoint"]
    comp=empirical["comp"]
    for index,(square,orbit) in enumerate(bands[2]):
        corners=[]
        for upper in (False,True):
            corners.append(gauge+[q(square if j==2 else row[-1 if upper else 0][0]).sqrt()
                                  for j,row in enumerate(bands)])
        segments=((arb(0),h3,160,adj["flow"](False,False)),(h3,h8,28,adj["flow"](True,False)),
                  (h8,t,12,adj["flow"](True,True)),(t,arb(33),64,adj["flow"](True,True)))
        if orbit==(24,10,1):
            empirical["triplet"]["verify_ground"](corners[0],segments,2,t,q("4.9e-8"),q(".136"),
                                                  upper_values=corners[1],logarithm=(start,h3,offsets))
            continue
        pivot=t+q((".00002123",".00001725")[index])
        hull=[]
        for upper,values in zip((False,True),corners,strict=True):
            state=[v*v for v in values]
            for part,(a,b,n,flow) in enumerate(segments[:3]):
                end=pivot if part==2 else b
                step=(end-a)/n
                for _ in range(n):
                    state=comp["advance"](state,step,upper,flow=flow)
            hull.append([1/(r*s)-d for r,s,d in zip(adj["ratios"],state[:3],offsets,strict=True)])
        # Between and outside thresholds the previous all-box slope bounds
        # retain their signs, including the ramp: each ramp slope is <.001.
        assert all(0<d/(h3-start)<q(".001") for d in offsets)
        radius=min(hull[0][1]-hull[1][0],hull[0][1]-hull[1][2])/2
        assert radius>q(".0000027")
        print("excluded top orbit",orbit,"all-target all-scale residual above",radius.lower().str(10))
        moved_logs=(("24.31764225325352","26.71338595599353","26.81726783749831"),
                    ("24.31737256915445","26.71326031426304","26.81726729304166"))[index]
        moved3,moved8,moved_q=map(q,moved_logs)
        new_energy,new_offsets,_,_=retuned["current"]["light_current"](moved_logs)
        new_segments=((arb(0),moved3,160,adj["flow"](False,False)),
                      (moved3,moved8,28,adj["flow"](True,False)),
                      (moved8,moved_q,12,adj["flow"](True,True)),
                      (moved_q,arb(33),64,adj["flow"](True,True)))
        empirical["triplet"]["verify_ground"](corners[0],new_segments,2,moved_q,q("5.0e-8"),q(".136"),
                 upper_values=corners[1],logarithm=(moved_q+new_energy.log()/2,moved3,new_offsets))
        assert 1<(moved3-h3).exp()<q("1.001440")
        assert 1<(moved8-h8).exp()<q("1.000671")
        print("recovered 80 tables with relative threshold shifts",
              ((moved3-h3).exp()-1).str(9),((moved8-h8).exp()-1).str(9))
    assert [orbit for _,orbit in bands[2]]==[(18,5,0),(15,14,3),(24,10,1)]
    print("80 tables match and 160 fail; unique top orbit (24,10,1)")


if __name__=="__main__":
    verify()

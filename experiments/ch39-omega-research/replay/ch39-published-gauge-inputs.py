"""Test the model using the 1 TeV marginal intervals of AHS, Table 2."""
from pathlib import Path
from runpy import run_path

from flint import arb

comp = run_path(str(Path(__file__).with_name("ch39-unrestricted-top-match.py")))
model, q, ratios = comp["model"], comp["q"], comp["ratios"]
centers = [*model["centers"],*model["yukawa_centers"]]
# Antusch, Hinze & Saad, arXiv:2510.01312v2, Table 2 (PDG 2024).
width_strings = (".000027",".00003",".0030",".00000014",".00005",".0043",
                     ".0000002",".000003",".00009",".0000000006",".00000013",".000002")
widths = list(map(q, width_strings))
anchor = comp["anchor"]


def corner(radius, upper):
    sign = 1 if upper else -1
    return ([(v-sign*radius*d)**2 for v,d in zip(centers[:3],widths[:3],strict=True)] +
            [(v+sign*radius*d)**2 for v,d in zip(centers[3:],widths[3:],strict=True)])


def endpoints(states, time):
    result = []
    for state,upper in zip(states, (False,True), strict=True):
        for _ in range(192):
            state = comp["advance"](state, time/192, upper)
        result.append([1/(r*g) for r,g in zip(ratios,state)])
    return result


def validate_slopes():
    running, loop = comp["running"],comp["loop"]
    # Up to ten printed half-widths, the same scalar bootstrap keeps X_i > 1/2.
    for i in range(3):
        lower = -2*(running["BETA"][i]+2*sum(max(v,0) for v in running["B"][i])/loop)/loop
        start = 1/(centers[i]+10*widths[i])**2
        assert min(start,start+33*lower) > arb(1)/2
        assert 1+lower > q(".01")*ratios[i]
    state = corner(2,True)
    step = arb(1)/16
    for i in range(33*16):
        state,tube = comp["advance"](state,step,True,need_tube=True)
        if i*step < 23:
            continue
        g = [arb(0).union(arb(2))]*3
        y = [arb(0).union(sum(tube[k:k+3]).upper()) for k in (3,6,9)]
        speed = [-2/(loop*r)*(b+sum(a*u-d*z for a,u,d,z in zip(br,g,dr,y,strict=True))/loop)
                 for r,b,br,dr in zip(ratios,running["BETA"],running["B"],running["D"],strict=True)]
        assert speed[0] < speed[1] < speed[2]


def verify():
    validate_slopes()
    for radius,scale,floor in ((1,"26.878761165276433",".0010454"),(1.5,"26.893842628902227",".0001322")):
        lo,hi = endpoints([corner(radius,u) for u in (False,True)],q(scale))
        assert min(lo[0]-hi[1],lo[2]-hi[1])/2 > q(floor)
    for radius,scale,floor in ((1,"26.81660044100034",".0038942"),(9,"26.816798660682565",".0003339")):
        tm = q(scale)
        lo,hi = endpoints([corner(radius,u) for u in (False,True)],tm)
        target = model["target"](11431,tm)
        assert min(lo[0]-target,target-hi[1]) > q(floor)+q("1.3e-16")

    for signs,scale,free,limit in (((1,-1,1),"26.8787430550987580",True,".0011101"),
                                   ((1,-1,0),"26.8166006451804967",False,".0039547")):
        gs = [c+k*w for c,k,w in zip(centers[:3],signs,widths[:3],strict=True)]
        point = [x*x for x in (*gs,*centers[3:])]
        tm = q(scale)
        lo,hi = endpoints([point,point],tm)
        target = q("3.297877826342") if free else model["target"](11431,tm)
        assert max(abs(x-target) for x in (*lo,*hi))+q("1.3e-16") < q(limit)
    gs = [(arb(v)/comp["mass"][0]).sqrt() for v in comp["mass"][1:]]
    deviations = [(v-c)/w for v,c,w in zip(gs,centers[:3],widths[:3],strict=True)]
    assert q("8.873") < deviations[0] < q("8.874") and q("-17.437") < deviations[1] < q("-17.436")
    print("four-word source gauge offsets in printed half-widths",*(d.str(10) for d in deviations))
    print("minimum maximum residual in the published box: free target (.0010454,.0011101), 11431 (.0038942,.0039547)")


if __name__ == "__main__":
    verify()

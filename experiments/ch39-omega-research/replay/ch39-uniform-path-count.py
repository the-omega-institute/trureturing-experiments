"""在规范项系数一、端口归一的模型中确定最小等权路径数。"""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

direction = run_path(str(Path(__file__).with_name("ch39-pure-port-direction-window.py")))
port = direction["port"]
study = port["study"]
geo = study["geo"]


def exclude_small_counts():
    possible = [(a,n-a) for n in range(2,26) for a in range(1,n)
                if Fraction(2014,1000)**2 < Fraction(n-a,a) < Fraction(2072,1000)**2]
    assert possible == [(4,17)]
    h = acb.zeta_zero(11431).imag
    phase = geo["phases"](h)
    ref,ref_der = direction["shifted_response"](phase[724],phase[804])
    zmax = arb(2)**-20
    shift = arb(0).union(zmax)
    tail = 4*zmax*zmax/(arb("1.3")-4*arb.pi()/32)**3
    points = [(5861,-5327),(2262,-401),(1667,415),(1161,1108),(775,1637),(-113,2856)]
    rows = []
    for i,j in points:
        vals,ders = direction["shifted_response"](phase[724]*acb(0,arb(i)/2**20).exp(),
                                                 phase[804]*acb(0,arb(j)/2**20).exp())
        rows.append([x-y+shift*(dx-dy)+arb(0,tail.upper())
                     for x,y,dx,dy in zip(vals,ref,ders,ref_der)])
    x = -arb(13)/21
    r = 4*arb(17).sqrt()/21
    stack = [(arb(".892"),arb(".920"))]
    boxes = 0
    while stack:
        lo,hi = stack.pop()
        angle = lo.union(hi)
        y,z = r*angle.cos(),-r*angle.sin()
        boxes += 1
        if any(t+x*d+y*u+z*v+1/study["center"] < 0 for t,d,u,v in rows):
            continue
        assert hi-lo>arb(2)**-28
        mid = (lo+hi)/2
        stack.extend(((lo,mid),(mid,hi)))
    print("All totals through 25 leave 4:17; excluded using",boxes,"phase intervals")


def construct(a,b,loops, *, unclocked):
    edges = {(i,j):d for i,j,d,_ in geo["edges"]}
    phase = 1
    for tri,weight in loops:
        i,j,k = tri
        lengths = [edges[e] for e in combinations(tri,2)]
        if unclocked:
            assert all(d not in (724,804,200) for d in lengths)
        z = geo["unit"](i,j)*geo["unit"](j,k)*geo["unit"](i,k).conjugate()
        if z.imag < 0:
            z = z.conjugate()
        assert z.imag > 0
        phase *= z**weight
    vector = acb_mat([[(arb(a)/(a+b)).sqrt()],[phase*(arb(b)/(a+b)).sqrt()]])
    gap,full,energy = port["verify"](vector,allowance=arb("4e-9"))
    assert gap > arb("4.13e-9") and full > arb("1e-15")
    assert energy > -arb(2)**-20
    print("Path groups",a,b,"clock-free",unclocked,"phase",phase.arg().str(14))


if __name__ == "__main__":
    ctx.prec = 192
    exclude_small_counts()
    construct(5,21,[((5,13,15),2),((5,7,13),2),((5,6,15),2)],unclocked=False)
    construct(5,21,[((5,13,14),1),((5,13,15),4),((5,13,17),1),((5,15,16),1)],unclocked=True)

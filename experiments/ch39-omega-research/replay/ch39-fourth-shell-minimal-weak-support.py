"""七条平凡源线和十二条非平凡目的线中，给出弱差的全部最小支撑。"""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    _,_,incoming,readout,energies,weights,_=lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    triple=lib("fourth-shell-three-line-weak")
    adj,norm=triple["adj"],triple["norm"]
    shell=lib("three-ray-chiral-boundary")
    n=lib("rank-two-generator-up-quark")["N"]
    mid=lib("adjacent-pvm-transport")["spectrum"](shell,6,n+3)[1]
    target=adj(readout)
    sites=shell["geo"]["sites"]
    rotation=[sites.index((y,z,x)) for x,y,z in sites]
    groups=lib("charge-six-projector-algebra-transport")["groups"]
    src=groups(mid,rotation)[0]
    dest=groups(target,rotation)
    other=sorted(dest[1]+dest[2])
    assert 15 in dest[0]
    lo,h,k=energies
    high=incoming*acb_mat([[(h-lo).sqrt(),0],[0,(k-lo).sqrt()]])
    old=acb_mat([[(readout*high)[15,j] for j in range(2)]])
    floor=lo+norm(old)**2
    limit=16*arb("91.1876")/arb("1.22089e19")
    assert floor>limit
    print("纯非平凡支撑保持的第16层振幅",floor.str(14),flush=True)
    for ntriv in (1,2):
        floors=[]
        for left in combinations(src,ntriv):
            for right in combinations(other,2-ntriv):
                e=acb_mat([[*[mid[i,j] for j in left],*[target[i,j] for j in right]] for i in range(19)])
                a=adj(e)*high
                co=readout*e
                u=acb_mat([[co[15,j] for j in range(2)]])
                b=old-u*a
                adjugate=acb_mat([[a[1,1],-a[0,1]],[-a[1,0],a[0,0]]])
                distance=abs(norm(b*adjugate)-norm(u)*abs(a.det()))/norm(a)
                floors.append((lo+distance**2).lower())
        assert len(floors)==comb(7,ntriv)*comb(12,2-ntriv) and min(floors)>limit
        print("双线平凡数、支撑数、第16层振幅下界",ntriv,len(floors),min(floors).str(14),flush=True)
    cases=triple["verify_data"](root,mid,incoming,readout,energies,weights,trivial_counts=(1,2,3))
    counts=[(comb(12,3),0)]
    for ntriv in (1,2,3):
        off,yes=cases[ntriv]
        assert len(off)+len(yes)==comb(7,ntriv)*comb(12,3-ntriv)
        counts.append((len(off),len(yes)))
    total_off,total_yes=map(sum,zip(*counts,strict=True))
    assert total_off+total_yes==comb(19,3) and total_yes>0
    print("三线分类（按平凡线数0至3，列为排除、可达）",counts,
          "合计",total_off,total_yes,flush=True)
    return counts,cases


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)

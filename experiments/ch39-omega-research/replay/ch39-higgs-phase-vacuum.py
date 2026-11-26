"""相位 Hessian 排除任意密度势真空给出实测的顶底层级。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb,acb_mat,arb,ctx


def gram(rot,case):
    w = acb_mat(19,19)
    entries = iter(case["gram"])
    for i in range(19):
        for j in range(i,19):
            orbit = {(i,j),tuple(sorted((rot[i],rot[j]))),tuple(sorted((rot[rot[i]],rot[rot[j]])))}
            if (i,j) == min(orbit):
                value = arb(next(entries))/10**6
                for a,b in orbit:
                    w[a,b] = w[b,a] = value
    assert next(entries,None) is None
    return w


def verify():
    root = Path(__file__).parent
    model = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    ctx.prec = 256
    up,down = flavor["flavor_maps"](model,11431)
    gu,gd = flavor["adjoint"](up)*up,flavor["adjoint"](down)*down
    lap = model["laplacian"](3,acb.zeta_zero(11431).imag)
    sites = model["geo"]["sites"]
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    identity = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    data = json.loads((root/"ch39-phase-vacuum-witnesses.json").read_text())
    for name,case in data.items():
        w = gram(rot,case)
        positive(w)
        v = [arb(z)/10**6 for z in case["first"]]
        diff = acb_mat([[(w[i,i]+w[j,j]-2*w[i,j]+acb(0,v[i]-v[j]))*lap[i,j]
                         for j in range(19)] for i in range(19)])
        coefficient = arb(403) if name == "hessian" else arb("414.8")
        positive(diff+1000*gd-coefficient*gu)
        if name == "hessian":
            assert all(z.is_zero() for z in v)
            trace = sum(w[i,i].real for i in range(19))
            ratio = arb(".8573")/arb(".01410")
            positive(diff+1000*(gd-3*gu/ratio**2)-arb(".006")*trace*identity)
        print("证书、相位稳定剖面的顶底比上界",name,(3000/coefficient).sqrt().str(12))


if __name__ == "__main__":
    verify()

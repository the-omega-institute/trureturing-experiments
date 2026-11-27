"""两个已配平的右侧局部荷六见证的全部酉立方根均排除轻上夸克。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def root_data(center, radius):
    values, vectors = center.eig(right=True, algorithm="rump")
    projectors = []
    for j in range(len(values)):
        v = acb_mat([[vectors[i,j]] for i in range(len(values))])
        projectors.append(v*v.conjugate().transpose()*(1/sum(abs(x)**2 for x in v.entries())))
    gaps = [min(abs(z-w).lower() for k,w in enumerate(values) if k!=j) for j,z in enumerate(values)]
    assert all(g > 2*radius for g in gaps) and all(abs(z-1)>radius for z in values)
    # 真实块的升起本征线与中心线的距离及三次根的相位扰动。
    error = arb(4)/3*(radius/2).asin()+sum(radius/(gap-radius) for gap in gaps)
    roots = [[acb(0, -arb(2)/3*(z.arg()+2*k*arb.pi())).exp() for k in range(3)] for z in values]
    return projectors, roots, error


def verify(root):
    grav = run_path(str(root/"ch39-minimal-support-right-gravity.py"))
    blocks = grav["verify"](root)
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    local = run_path(str(root/"ch39-up-eight-site-cover.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    positive = run_path(str(root/"ch39-up-quark-arbitrary-site-phase-floor.py"))["positive"]
    n = 2_306_888_100_086
    _, labels, _ = train["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    _, old = spec(shell,-4,n)
    _, target = spec(shell,-4,n+1)
    w = labels[n,-4]
    assert all(0 < x < arb(".947555") for x in w)
    h = weighted(old,w)
    four = local["rows"](target,(14,8,7,6))
    for (sites,center,radius), goal in zip(blocks,(".000013",".000040"),strict=True):
        fixed = sorted(sites)
        free = [i for i in range(19) if i not in fixed]
        order = [sites.index(i) for i in fixed]
        sub = acb_mat([[center[i,j] for j in order] for i in order])
        projectors, roots, error = root_data(sub,radius)
        threshold = (arb(goal).sqrt()+arb(".947555").sqrt()*error)**2
        dim = len(fixed)
        count = 0
        for flags in product(range(3),repeat=dim-1):
            # 统一乘三次单位根不改读数，正交补已容许任意酉。
            phases = [row[k] for row,k in zip(roots,(0,*flags),strict=True)]
            up = sum((p*z for p,z in zip(projectors,phases,strict=True)),acb_mat(dim,dim))
            local["certify_set"](free,h,target,four,positive,threshold,up.conjugate().transpose())
            count += 1
        print("站位、荷六块半径、上型块误差上界、分支数、轻上型严格下界",
              sites,radius.str(10),error.str(12),count,goal,flush=True)
    return blocks


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

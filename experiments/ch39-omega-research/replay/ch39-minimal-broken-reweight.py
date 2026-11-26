"""在固定五位置上计算把两项引力和电磁—色匹配移到电弱背景的最小重赋权。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    joint = run_path(str(root / "ch39-all-field-signed-crossing.py"))
    moments = run_path(str(root / "ch39-signed-cutoff-gravity.py"))["moments"]
    vector_log = run_path(str(root / "ch39-em-colour.py"))["vector_log_moment"]
    center = arb("1.45137757e14")
    energy = (center-100000).union(center+100000)
    data = joint["joint_system"](root, energy, order=17, last=49)
    old = data["matrix"].solve(data["target"]).entries()
    pi = arb.pi()
    r = (arb("246.22")/arb("1.22089e19"))**2/2
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    kw, kz = g2**2/2, (g2**2+3*g1**2/5)/2
    rows = data["rows"]
    heavy = [z for _, _, z in rows[49:65]]
    shift = (-r).union(r)
    moved = [(-59,arb(2)/3,arb(0)), (-12,1,yt**2*r), (6,-1,kw*r), (3,-arb(1)/2,kz*r)]
    moved.extend(rows[1:49]+rows[65:])
    for g, a, z in rows[49:65]:
        moved.extend(((g-4,a-arb(2)/3,z), (4,arb(2)/3,z+shift)))
    scalar_error = r*sum(1/(z-r) for z in heavy)/(12*pi)
    low_moment = vector_log(kw*r/data["p0"])
    ratio = kw*r/data["p1"]
    high_moment = vector_log(ratio.lower()).union(vector_log(ratio.upper()))
    cols, changes = [], []
    for x, origin in zip(data["sites"], data["cols"], strict=True):
        k, w = moments(moved, data["heat"], x)
        delta_b = 7/(4*pi)*(high_moment-low_moment+arb(0).union(kw*r/x))
        delta_b += arb(0,scalar_error.upper())
        b0 = 5*origin[3]/3+origin[4]-8*origin[5]/3
        cols.append([1,k,w,b0+delta_b])
        changes.append((k-origin[1],w-origin[2],delta_b))
    target = arb_mat([[0], *[[-sum(c*d[i] for c,d in zip(old,changes,strict=True))]
                             for i in range(3)]])
    choices = []
    for fixed in range(5):
        active = [i for i in range(5) if i!=fixed]
        mat = arb_mat([cols[i] for i in active]).transpose()
        assert not mat.det().contains(0)
        values = mat.solve(target).entries()
        delta = [arb(0)]*5
        for i, val in zip(active, values, strict=True):
            delta[i] = val
        cost = sum(abs(v) for v in delta)/2
        choices.append((cost,delta))
        print("保持第",fixed+1,"项权的重赋权代价",cost.str(12))
    best = min(range(5), key=lambda j:float(choices[j][0].mid()))
    cost, delta = choices[best]
    print("最佳固定项、半一范数及五项权增量", best+1,cost.str(14),
          *[v.str(12) for v in delta])
    assert best == 1 and arb(".062252") < cost < arb(".062254")
    assert all(cost < v for i,(v,_) in enumerate(choices) if i!=best)


if __name__ == "__main__":
    verify()

"""真实三源上，两部门有符号谱压缩给出正定平方质量所需的负迹。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def squared_norm(v):
    return sum(abs(z)**2 for z in v)


def cone_constants(rays, selected):
    # h_j 是 (1-vv*)P_j v，s_j 是该投影的横向迹。
    v = rays[selected].entries()
    zero = [acb(0) for _ in v]
    columns = []
    slopes = [zero]
    for index, ray in enumerate(rays):
        if index == selected:
            columns.append((arb(0),zero))
            continue
        x = ray.entries()
        phase = sum(z.conjugate()*y for z,y in zip(v,x,strict=True))
        perp = [y-z*phase for y,z in zip(x,v,strict=True)]
        s = squared_norm(perp)
        assert s > 0
        h = [y*phase.conjugate() for y in perp]
        columns.append((s,h))
        slopes.append([y/s for y in h])
    c = max(squared_norm(a).upper() for a in slopes).sqrt().upper()
    b = max(squared_norm([s*y-z for y,z in zip(a,h,strict=True)]).upper()
            for a in slopes for s,h in columns).sqrt().upper()
    return c,b


def tangent_bound(coefficients, ratio, negative):
    c,b = coefficients
    gap = 1-ratio
    h = c*ratio+b*negative
    discriminant = gap**2-4*(c+h)*h
    assert discriminant > 0
    root = discriminant.sqrt()
    # 尾迹给 tan(theta)≤sqrt(r+b)，因此根间的禁区只留小根分支。
    assert (ratio+negative).sqrt() < (gap+root)/(2*(c+h))
    small = 2*h/(gap+root)
    assert small < arb(".1")
    return small.upper()


def verify(root):
    model = run_path(str(root/"ch39-three-source-positive-ckm.py"))
    up,down = model["source_rays"](root)
    distance = model["distance"]
    ctx.prec = 320
    near,far,exceptional = [],[],{}
    for i,u in enumerate(up):
        for j,d in enumerate(down):
            angle = distance(u,d)
            if (i,j) in ((0,0),(6,6)):
                exceptional[i] = angle
            elif angle < arb("1e-15"):
                near.append(angle)
            else:
                far.append(angle)
    assert len(near)==72 and len(far)==287 and set(exceptional)=={0,6}
    near_max,far_min = max(z.upper() for z in near),min(z.lower() for z in far)
    assert far_min > arb(".1360203")
    constants = {i:[cone_constants(r,i) for r in (up,down)] for i in exceptional}
    for i,cd in constants.items():
        print("所选谱层、两部门的横向斜率与负量系数",i+1,
              [[v.str(12) for v in row] for row in cd])

    # AHS 表 2 的 1 TeV 六个 Yukawa 边际 1σ HPD 窗，取对排除最宽松的比。
    ratios = [(arb(".00316")**2+arb(".00000629")**2)/arb(".8573")**2,
              (arb(".000271")**2+arb(".0000137")**2)/arb(".01392")**2]
    s23,s13 = arb(".04183"),arb(".003732")
    half = arb(".000055")
    targets = []
    for sign in (-1,1):
        a,b = s23+sign*half,s13+sign*half
        targets.append((b*b+(1-b*b)*a*a).sqrt().asin())
    ahs = [(arb(a).cos()*arb(b).cos()).acos() for a,b in (
        (".04233",".00369"),(".04317",".00385"))]

    def coarse(budgets,target):
        total = sum((r+b).sqrt().atan() for r,b in zip(ratios,budgets,strict=True))
        assert near_max+total < target[0] and far_min-total > target[1]

    def exceptions(budgets,target):
        for i,angle in exceptional.items():
            bound = angle+sum(tangent_bound(cs,r,b) for cs,r,b in
                              zip(constants[i],ratios,budgets,strict=True))
            assert bound < target[0],(i,budgets,bound,target)
            if i == 0:
                print("第一层配第一层的重线夹角上界",bound.upper().str(13))

    for label, target, equal, upper, lower in (
        ("书中 CKM",targets,".0000689",".0002146",".0001012"),
        ("同尺度表 2",ahs,".0001106",".0003444",".0001619")):
        equal,upper,lower = map(arb,(equal,upper,lower))
        print(label,"等预算、仅上型、仅下型的排除预算",
              equal.str(8),upper.str(8),lower.str(8))
        coarse((equal,equal),target)
        exceptions((equal,equal),target)
        # 在带权预算三角形中，小根是负量的递增凸函数。两坐标轴端点控制最大和。
        coarse((upper,lower),target)
        exceptions((upper,arb(0)),target)
        exceptions((arb(0),lower),target)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

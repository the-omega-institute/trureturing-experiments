"""最小割三窗内非负绕数三角作用取得全环面极小的最小阶。"""
from collections import defaultdict
from itertools import combinations
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb,arb,ctx

def terms(shell):
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    edges = {}
    for i, j, d, transport in shell["geo"]["edges"]:
        spoke = pole in (i, j)
        softened = int(spoke and (i in tuned or j in tuned))
        x, y = int(not spoke and d == 724), int(not spoke and d == 804)
        z = acb(1) if spoke else transport
        edges[i,j] = z, x, y, softened
        edges[j,i] = z.conjugate(), -x, -y, softened
    triangles = defaultdict(list)
    for i, j, k in combinations(range(19), 3):
        cycle = [(i,j), (j,k), (k,i)]
        if not all(e in edges for e in cycle):
            continue
        phases, xs, ys, weights = zip(*(edges[e] for e in cycle), strict=True)
        z, x, y = prod(phases), sum(xs), sum(ys)
        if x < 0 or (x == 0 and y < 0):
            x, y, z = -x, -y, z.conjugate()
        triangles[x,y].append((z, sum(weights)))
    return triangles


def verify_physical(root, x, y, radius):
    physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
    model, initial = physical["model_inputs"]()
    ctx.prec = 224
    joint = physical["joint"]
    chiral = joint["chiral"]
    phases = {724: acb(0,x).exp(), 804: acb(0,y).exp()}
    # 供谱例程对所有高度大于 10000 的软化量使用共同的外包围。
    rows = chiral["source"]["spectral_supplier"](model)(arb(10)**13, clocks=phases)
    bands = {}
    for q in (0, *model["charges"]):
        values = chiral["ordered"]([v for v,charge in rows if q == charge])
        bands[q] = [v + arb(0,320*q*radius) for v in values]
    assert all(bands[q][2] < arb("9.99999999") and bands[q][3] > 10 for q in model["charges"])
    _, _, alpha = joint["gauge"]["crossing_enclosure"](
        model, initial, bands, False, 724, "29.04", ".1184", recenter=True)
    heavy = 32*arb("2.435e18")/724**2*bands[2][15]
    ratio = (arb("246.22")**2/(heavy*arb("1.22089e19")))**2
    vacuum = joint["vacuum"]
    paired = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](bands), min_gap="1e-18")
    fields = paired + vacuum["low_sector"](4,ratio)
    assert vacuum["newton"](fields,arb(".03")) < 1 < vacuum["newton"](fields,arb(".033"))
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    curve = ratio**2*coefficient/(arb.pi()*arb("1e-122"))*(1-arb(0).union(ratio/(arb(".03")*coefficient)))
    assert arb(".11844168") < alpha < arb(".11845151")
    assert arb(".99949733532") < curve < arb(".99949733536")
    print("极小点方框的强耦合、曲率比", alpha.str(12),curve.str(14))


def coeff(tri,q,delta):
    result = []
    for key in ((1,0),(0,1),(1,-1)):
        group = tri[key]
        bare = sum(v**q for v,w in group)
        response = sum(v**q for v,w in group if w==1)
        assert all(w in (0,1) for v,w in group)
        result.append(bare-delta*response)
    return result


def energy(q, coefficients, x, y):
    rotations = [acb(0,q*t).exp() for t in (x,y,x-y)]
    return sum((a*(1-z)).real for a,z in zip(coefficients, rotations, strict=True))


def lower_orders(model,tri):
    delta = arb(0).union(arb(23)/25*arb(2)**-20/14**2)
    floor = arb(".023")
    radius = arb(".10124")
    assert radius/arb(2).sqrt() > arb(".07") and 25*arb.pi()/16 < 2*arb.pi()-arb(".07")
    curvature_upper = arb("18.944146")
    assert model["sign"]["density"]["eigen_band"]()[1] < curvature_upper
    count = model["cell_count"](2,curvature_upper)
    pending = [(1,0,0)]
    while pending:
        n,i,j = pending.pop()
        a,b = [max(0, abs(2*k+1-n)-1) for k in (i,j)]
        if (a*a+b*b)*(arb.pi()/(16*n))**2 > radius**2 or count(n,i,j)[1] <= 15:
            continue
        assert n < 2**19
        pending.extend((2*n,2*i+ii,2*j+jj) for ii in (0,1) for jj in (0,1))
    for q in range(1,26):
        a,b,d = coeff(tri,q,delta)
        p,r,s = a.real+d.real-floor,b.real+d.real-floor,d.real
        assert p>0 and p*r-s*s>0
        error = sum(abs(c.imag).upper() for c in (a,b,d))
        assert floor*2*(arb(".035").sin())**2-error > arb(".000055")
    # 对第 26 圈的能量亏额，零相位为同一个竞争点。在不能比较的格子才检查荷二曲率。
    c = coeff(tri,26,delta)
    cap = [arb("18.94088212"),arb("18.94414547")]
    band = model["sign"]["density"]["eigen_band"]()
    assert cap[0]<band[0]<band[1]<cap[1]
    below,above = [model["cell_count"](2,t) for t in cap]
    pending=[(1,0,0)]
    leaves=[0,0]
    while pending:
        n,i,j=pending.pop()
        x,y=[arb.pi()/16*(arb(2*k+1,1)/n-1) for k in (i,j)]
        if energy(26,c,x,y)>arb(".10"):
            leaves[0]+=1
        elif below(n,i,j)[0]>=16 or above(n,i,j)[1]<16:
            leaves[1]+=1
        else:
            assert n<2**18, (n,i,j)
            pending.extend((2*n,2*i+ii,2*j+jj) for ii in (0,1) for jj in (0,1))
    print("首25阶、26阶相对零相位的统一严格下界，26阶的两类叶格", ".000055", ".10", leaves)


def derivatives(modes,x,y):
    gx=gy=hxx=hxy=hyy=arb(0)
    for q,cs in modes:
        u,v,w=[c*acb(0,q*t).exp() for c,t in zip(cs,(x,y,x-y),strict=True)]
        gx+=q*(u+w).imag
        gy+=q*(v-w).imag
        hxx+=q*q*(u+w).real
        hxy-=q*q*w.real
        hyy+=q*q*(v+w).real
    return gx,gy,hxx,hxy,hyy


def global_cover(tri):
    # 高度至少 10^12。用有界辐边亏额（含零极限）同时包住整个族。
    delta=arb(0).union(arb(23)/25*arb(2)**-20/10**24)
    modes=[(1,coeff(tri,1,delta)),(27,[arb(13)/100*c for c in coeff(tri,27,delta)])]
    center = [arb(".0707075643615173"),arb(".07570197643661009")]
    radius=arb(".0001")
    def value(x,y):
        return sum(energy(q,c,x,y) for q,c in modes)
    trial=value(*center)
    assert trial < 0
    pending=[(1,0,0)]
    leaves=[0,0,0]
    while pending:
        n,i,j=pending.pop()
        x,y=[arb.pi()*(arb(2*k+1,1)/n-1) for k in (i,j)]
        if value(x,y)>trial:
            leaves[0]+=1
            continue
        if any(abs(x-s*center[0])<radius and abs(y-s*center[1])<radius for s in (-1,1)):
            leaves[1]+=1
            continue
        gx,gy,*_=derivatives(modes,x,y)
        if not gx.contains(0) or not gy.contains(0):
            leaves[2]+=1
            continue
        assert n<2**22, (n,i,j)
        pending.extend((2*n,2*i+ii,2*j+jj) for ii in (0,1) for jj in (0,1))
    # 最小点必在这两个方框内。强凸性及中心残差将每个可能的最小点缩至更小方框。
    for sign in (1,-1):
        x,y=[sign*c for c in center]
        _,_,a,d,b=derivatives(modes,x+arb(0,radius),y+arb(0,radius))
        h=arb(10)
        assert a>h and (a-h)*(b-h)>d*d
        gx,gy,*_=derivatives(modes,x,y)
        dist=(gx.abs_upper()**2+gy.abs_upper()**2).sqrt()/h
        assert dist < arb("2.1e-17")
        print("全局极小点候选距方框中心的上界",dist.str(12))
    print("全环面能量、候选方框、梯度三种叶格数",leaves,"比较点能量差",trial.str(16))
    return center


def verify(root):
    ctx.prec=192
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    model=run_path(str(root/"ch39-mincut-three-window-angle.py"))
    ctx.prec=192
    tri=terms(shell)
    assert sum(map(len,tri.values()))==529
    lower_orders(model,tri)
    center=global_cover(tri)
    for s in (-1,1):
        verify_physical(root,*(s*c for c in center),arb("3e-14"))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

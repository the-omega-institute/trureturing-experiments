"""在同一实际零点定位荷六第三层与荷三第四层的共同谱隙闭合。"""
from fractions import Fraction as F
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


LEFT = F("-0.02908853997597")
RIGHT = F("-0.02908853997584")
POLE = (17,17,17)
SEED = (7,17,23)


def ball(x):
    return arb(x.numerator)/x.denominator


def geometry_interval(sites,t,*,reference=False):
    v = 1698*t*t
    seed = tuple(((867-v)*a+1734*t*b)/(867+v) for a,b in zip(SEED,(40,-7,-7)))
    moved = dict(zip(permutations(SEED),permutations(seed),strict=True))
    points = [moved.get(p,p) for p in sites]
    assert all(min(p)>0 for p in points)
    x0 = sites.index(POLE)
    def squared(x,y):
        return sum((a-b)*(a-b) for a,b in zip(x,y))
    distances = {}
    for i,j in combinations(range(19),2):
        before = squared(sites[i],sites[j])
        after = squared(points[i],points[j])
        distances[i,j]=distances[j,i]=after
        assert after>0,(t,i,j,points[i],points[j],after)
        if sites[i] not in moved and sites[j] not in moved:
            assert before==after
            continue
        assert (after<804 if before<804 else after>804)
        assert before not in (724,804) and (after<724 or after>724)
        if reference:
            assert before!=244 and (after<244 or after>244)
        if x0 in (i,j):
            assert 68<after<544
    assert all(any(distances[i,j]>544 for j in range(19) if j!=i) for i in range(19) if i!=x0)


def verify(root):
    base = run_path(str(root/"ch39-flex-two-loop-curvature.py"))
    shape = run_path(str(root/"ch39-metric-pair-flexibility.py"))
    physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
    model,initial = physical["model_inputs"]()
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    read = run_path(str(root/"ch39-scale-free-geometric-pair.py"))["readout"]
    graph = shell["geo"]
    sites = graph["sites"]
    geom = read(sites)
    assert geom[0]==(32,724) and geom[2]==(32,724,804)
    ctx.prec,ctx.cap = 256,3
    # 从原点到临界盒的整段核验图、钟边、中心及最近/最远轨道。
    for k in range(48):
        interval = ball(LEFT*k/48).union(ball(LEFT*(k+1)/48))
        geometry_interval(sites,interval)
    geometry_interval(sites,ball(LEFT).union(ball(RIGHT)),reference=True)
    # 再向负方向延伸时，一条原平方弦长 512 的边会穿过 724。
    a,b = F("-0.032989"),F("-0.032988")
    assert F(-1,25)<a<b<LEFT
    i,j = sites.index((7,17,23)),sites.index((23,17,7))
    assert shape["distance"](sites[i],sites[j])==512
    def moving_distance(t):
        p = shape["deform"](sites,t,t)
        return shape["distance"](p[i],p[j])
    assert moving_distance(a)>724>moving_distance(b)
    print("原512边穿过724的参数夹逼、端点及原坐标",a,b,sites[i],sites[j])

    h = acb.zeta_zero(base["NUMBER"]).imag
    phases = graph["phases"](h)
    soft = shell["soft"]
    pole,tuned = soft["pole"],soft["star"]["orbit"]
    softened = {(i,j) for i,j in geom[1] if pole in (i,j) and (i in tuned or j in tuned)}
    delta = arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
    mean = arb(15)-delta/3
    cut = 9+(1-delta)**2
    ordered = physical["joint"]["chiral"]["ordered"]
    # 立体投影速度不超过 2 sqrt(1698)。正坐标保证未归一化复振幅的模至少为 867²。
    speed = arb(83)
    assert arb(1698).sqrt()*2<speed
    derivative = 4*867*30*speed+2*30**3*speed
    assert derivative<20*867**2 and arb(867).sqrt()<30

    def layers(t,charges,radius=F(0)):
        p = shape["deform"](sites,t,t)
        assert read(p)==geom
        links = base["links"](graph,p,phases,softened,delta)
        levels = {0:base["neutral_levels"](geom[1],delta)}
        for q in charges:
            eig = base["laplacian"](links,q).eig(multiple=True)
            assert len(eig)==19 and all(z.imag.contains(0) and z.real>0 for z in eig)
            # 权和对角恒定，归一化联络速度上界 20，Weyl 用最大绝对行和。
            width = 18*20*q*ball(radius)
            levels[q] = [z.real+arb(0,width) for z in eig]
        return {q:ordered([v*abs(v-mean) for v in row]) for q,row in levels.items()}

    for t in (LEFT,RIGHT):
        row = layers(t,(3,6))
        gap = row[3][3]-row[6][2]
        if t==LEFT:
            assert arb("-7.057e-12")<gap<arb("-7.056e-12")
        else:
            assert arb("7.037e-12")<gap<arb("7.038e-12")
        print("有理端点、荷三第四减荷六第三",t,gap.str(16))
    row = layers((LEFT+RIGHT)/2,model["charges"],(RIGHT-LEFT)/2)
    crossing = row[3][3].union(row[6][2])
    assert arb("9.78116")<crossing<arb("9.78117")<cut
    assert all(row[q][2]<crossing for q in model["charges"] if q!=6)
    assert all(row[q][3]>crossing for q in model["charges"] if q!=3)
    assert row[3][4]>cut and row[6][1]<crossing
    assert all(all(x<cut or x>cut for x in row[q]) for q in model["charges"])
    counts = [sum(x<cut for x in row[q]) for q in model["charges"]]
    assert list(model["charges"])==[1,4,2,3,6] and counts==[3,3,3,4,3]
    for q in model["charges"]:
        print("临界盒全域的荷与第三、第四层",q,row[q][2].str(13),row[q][3].str(13))
    print("最小割以下的实际轻模数",counts)

    # 在谱隙的边界上，按各荷最低三层延拓既有阈值及四层场内容。
    joint = physical["joint"]
    vacuum = joint["vacuum"]
    mass = 32*arb("2.435e18")/724**2*row[2][15]
    A = (arb("246.22")**2/(mass*arb("1.22089e19")))**2
    paired = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](row),min_gap="1e-18")
    rows = paired+vacuum["low_sector"](4,A)
    assert vacuum["newton"](rows,arb(".03"))<1<vacuum["newton"](rows,arb(".033"))
    C = (9*arb(3).log()-12*arb(2).log())/2
    curve = C*A*A/(arb.pi()*arb("1e-122"))*(1-arb(0).union(A/(arb(".03")*C)))
    assert arb(".99741813")<curve<arb(".99741814")<arb(".998")
    print("荷二第十六层、四层曲率比",row[2][15].str(14),curve.str(14),flush=True)
    energy,_,strong = joint["gauge"]["crossing_enclosure"](model,initial,row,False,724,"29.04",".1184",recenter=True)
    assert arb(".11839139")<strong<arb(".11840854")
    assert arb("3.72805e14")<energy<arb("3.73039e14")
    print("秩三延拓的两圈交点能标与强耦合区间",energy.lower().str(13),energy.upper().str(13),
          strong.lower().str(13),strong.upper().str(13))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

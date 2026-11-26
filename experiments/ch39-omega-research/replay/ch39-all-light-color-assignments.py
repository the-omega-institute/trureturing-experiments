"""用全部轻夸克的最少标量分量修复规范交点的一圈有限分类。"""
from collections import Counter
from itertools import permutations, product
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


gauge = run_path(str(Path(__file__).with_name("ch39-portal-color-sextet.py")))["gauge"]
choices = ((0,0),(1,1),(2,2),(0,3),(2,5))


def shell_table(model, bands):
    unit = arb("2.435e18") / 724**2
    def integers(coefficients):
        return [int((30*v).unique_fmpz()) for v in coefficients]
    rows = []
    def add(levels, factor, coefficients):
        w = integers(coefficients)
        rows.extend(((factor*unit*h/model["mz"]).log(), w) for h in levels)
    for q, weyl in zip(model["charges"],model["weyl"]):
        add(bands[q][3:], 32 if q%2==0 else 1, [2*v for v in weyl])
    add(bands[3][3:], 1, model["higgs"])
    add(bands[0], 32, [arb(0),arb(-7),-arb(21)/2])
    return sorted(rows,key=lambda row:row[0].mid())


def segments(model, initial, rows):
    lower = (arb("1e13")/model["mz"]).log().lower()
    upper = (arb("1e17")/model["mz"]).log().upper()
    # 以精确二进中点为折点。每项正部对折点 1-Lipschitz，统一补足各分量误差。
    uncertainty = [sum(abs(w[i])*t.rad()/30 for t,w in rows) for i in range(3)]
    beta = [int((30*v).unique_fmpz()) for v in model["beta"]]
    shift = [2*arb.pi()*initial[0],2*arb.pi()*initial[1],arb(0)]
    top = (arb("172.5")/model["mz"]).log()
    beta[0] += 4
    beta[2] += 25
    shift[0] += 2*top/15
    shift[2] += 5*top/6
    previous = lower
    for t, weight in rows + [(upper,[0,0,0])]:
        center = t.mid()
        end = min(center,upper)
        if end > previous:
            intervals = [s+arb(0,r) for s,r in zip(shift,uncertainty)]
            yield previous,end,intervals,tuple(beta)
            previous = end
        if previous >= upper:
            break
        for i,w in enumerate(weight):
            beta[i] += w
            shift[i] += w*center/30


def cover(model,initial,bands):
    configurations = Counter(tuple(map(sum,zip(*p))) for p in product(choices,repeat=5))
    assert sum(configurations.values())==3125 and len(configurations)==61
    candidates = []
    target = arb(".117").union(arb(".119"))
    table = shell_table(model,bands)
    for lo,hi,s,b in segments(model,initial,table):
        span = lo.union(hi)
        for (A,B),count in configurations.items():
            slope = b[0]-b[1]+2*A
            delta = s[0]-s[1]
            if slope==0:
                assert not delta.contains(0)
                continue
            root = 30*delta/slope
            if not root.overlaps(span):
                continue
            inverse = (s[1]-s[2] + (b[2]+5*B-b[1])*root/30)/(2*arb.pi())
            if not inverse.overlaps(1/target):
                continue
            assert lo<root<hi and all(not t.overlaps(root) for t,_ in table)
            assert target.contains_interior(1/inverse)
            unified = (s[1]-b[1]*root/30)/(2*arb.pi())
            assert unified>0
            candidates.append((A,B,count))
            vectors = sum(1 for h in bands[0] if 32*arb("2.435e18")/724**2*h < model["mz"]*root.exp())
            print("有效指标、带标签数、能标、强耦合、统一逆耦合、开启中性重层",
                  A,B,count,(model["mz"]*root.exp()).str(11),(1/inverse).str(10),unified.str(10),vectors)
    assert candidates==[(3,15,20)]


def check_tangencies(assignments):
    sm = gauge["sm"]
    r = Fraction(1,1000)
    g1,_,g3 = sm["hull_inputs"](".467453 .63811 1.0583",r)
    ys = sm["hull_inputs"](".00000615 .00311 .8616 .0000135 .000268 .01401 .0000028227 .00059465 .010101",r)
    p = [y**2 for y in ys]
    hyper = 3*g1**2/5
    source = [16*g3**2+4*hyper/3]*3+[16*g3**2-2*hyper/3]*3+[6*hyper]*3
    # CKM 的每一行和列平方和均为一，故覆盖每个允许的 X_i。
    mixing = [arb(0).union(p[5])]*3 + [arb(0).union(p[2])]*3 + [arb(0)]*3
    reps = {"1":[(12,0,0)],"3":[(6,arb(4)/3,arb(1)/9),(6,0,0)],
            "8":[(8,3,0),(4,0,0)],"6":[(12,arb(10)/3,arb(1)/9)]}
    trace = sum(n*u for n,u in zip([12]*6+[4]*3,p))
    for others in assignments:
        blocks = []
        fields = [*others[:2],"6",*others[2:]]
        for i in range(9):
            for n,casimir,ysq in reps[fields[i]] if i<6 else [(4,0,0)]:
                delta = source[i] - 6*(casimir*g3**2 + ysq*hyper)
                residual = -p[i]*(p[i]+3*mixing[i]+arb("1.5")+delta)
                residual -= 3*ysq*hyper**2 + arb(".001")*(trace+2*p[i])
                blocks.append((n,p[i],residual))
        def image(A):
            assert all(A**2+8*u*w>0 for n,u,w in blocks)
            return sum(n*u*2*w/(A+(A**2+8*u*w).sqrt()) for n,u,w in blocks)
        assert image(arb(2))>2 and image(arb(5))<5
    print("可同时匹配门户与平方 Yukawa 一阶导数的带标签分配",len(assignments))


def verify():
    chiral=gauge["chiral"]
    model,_,_,initial,_=chiral["source"]["mean"]["inputs"]()
    ctx.prec=224
    assignments = set(permutations(("3","8","8","8","6")))
    assert len(assignments)==20
    check_tangencies(assignments)
    read=chiral["source"]["spectral_supplier"](model)
    for n in (22027,28325):
        src=read(acb.zeta_zero(n).imag)
        bands={q:chiral["ordered"]([h for h,z in src if z==q]) for q in (0,*model["charges"])}
        print("正零点",n)
        cover(model,initial,bands)


if __name__=="__main__":
    verify()

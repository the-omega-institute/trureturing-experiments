"""精确参考味的全局三窗真空所需的绕数权总变差，及三个见证的精确距离。"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def prefix_margin(root):
    old = run_path(str(root / "ch39-monotone-winding-action.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    exact = run_path(str(root / "ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 144
    bands = [(arb(a),arb(b)) for a,b in old["BANDS"]]
    ds = old["verify_phases"](root, shell, bands)
    distance = ds[0].union(ds[1])
    counts = Counter(v for v in exact["enumerate_faces"]().values() if v[0] != (0,0))
    assert sum(counts.values()) == 96
    rows = [(i,j,z,exact["unit"](z).arg(),n) for ((i,j),z),n in counts.items()]
    locks = old["LOCKS"]
    u,v = [-exact["unit"](locks[k]).arg() for k in ((1,0),(0,1))]

    def rival(x,y,aligned):
        result = []
        for i,j,z,a,n in rows:
            angle = a+i*x+j*y
            if ((i,j),z) in aligned:
                # 精确锁定的相位：所有固定竞争点的活动角只能等于 u 或 v。
                assert angle.contains(0)
                angle=arb(0)
            result.append((angle,n))
        return result

    # 三个竞争点的混合系数为 2/5、1/5、2/5。
    rivals = [rival(u,v,set(locks.items())), rival(u,arb(0),{((1,0),locks[1,0])}),
              rival(v,v,{((0,1),locks[0,1])})]
    mixture = [(2,rivals[0]),(1,rivals[1]),(2,rivals[2])]
    geom = old["geometric_tail"]
    stats = [geom(r,identically_zero=True) for r in rivals]
    assert [n for n,_ in stats] == [6,3,3]
    tail_strength = sum(a*n for a,(n,_) in zip((2,1,2),stats))/arb(5)
    tail_error = sum(a*e for a,(_,e) in zip((2,1,2),stats))/5
    correction = 12*arb("4.48e-9")
    advertised = arb(13)/50
    margin = advertised + correction
    end = 15000
    groups = [(angle,n*factor) for factor,r in mixture for angle,n in r]
    prefix = [arb(0)] + [old["kernel_sum"](groups,n)/5 for n in range(1,end)]
    smallest = arb(100)
    for l,r in bands:
        for lift in (0,1):
            def angles(t):
                out = []
                for i,j,_,a,n in rows:
                    theta = a+(i+j)*t-j*distance
                    if lift and i+j:
                        theta += -arb.pi() if a>0 else arb.pi()
                    out.append((theta,n))
                return out

            near,err = geom(angles(l.union(r)))
            bound = tail_strength-near-(tail_error+err)/end
            assert bound > margin
            smallest = smallest.min(bound)
            pending = [(l,r,range(1,end))]
            while pending:
                left,right,orders = pending.pop()
                cells = angles(left.union(right))
                undecided = []
                for n in orders:
                    gap = prefix[n]-old["kernel_sum"](cells,n)
                    if gap > margin:
                        smallest = smallest.min(gap)
                    else:
                        undecided.append(n)
                if undecided:
                    assert right-left>arb("1e-10"),(left,right,undecided[:8])
                    mid=(left+right)/2
                    pending.extend(((left,mid,undecided),(mid,right,undecided)))
    assert smallest-correction > advertised
    print("混合竞争点归一权、端点和尾部的全前缀裕量", "2/5,1/5,2/5", (smallest-correction).lower().str(15))
    print("任意有限或无限绕数权的总变差下界", F(13,50*384))


def weight_36(lib, tri, x, y, delta):
    orders = (1,2,30,32,33,36)
    grad = {q: lib["derivatives"]([(q,lib["coeff"](tri,q,delta))],x,y)[:2] for q in orders}
    fixed = {1:arb(".03744"),2:arb(".02331"),30:arb(".08198")}
    rhs = [1-sum(fixed.values())]+[-sum(w*grad[q][k] for q,w in fixed.items()) for k in range(2)]
    matrix = arb_mat([[1,*grad[q]] for q in orders[3:]]).transpose()
    free = matrix.solve(arb_mat([[t] for t in rhs])).entries()
    return {**fixed,**dict(zip(orders[3:],free,strict=True))}


def bound_expectation(test, maximum):
    # 指标以外的项恒为零；超过最大的支撑序号后，全部前缀均值递减。
    total = F(0)
    values=[]
    for n in range(1,max(test)+1):
        total+=test.get(n,0)
        values.append(total/n)
    assert max(values)==maximum


def exact_examples(root):
    common=run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    winding=run_path(str(root / "ch39-mincut-global-winding.py"))
    local=run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    flavor=run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    shell=run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec=192
    tri=winding["terms"](shell)
    x=arb("-.01046")
    y=x-local["reference_angles"](flavor)[0]
    delta=arb(0).union(arb("4.48e-9"))
    pr=weight_36(winding,tri,x,y,delta)
    h=(1-pr[1])/35
    assert pr[1]>h>pr[2] and all(pr[q]>h for q in (30,32,33,36))
    high={q:F(1) for q in (30,32,33,36)}
    bound_expectation({**high,1:F(4,35)},F(4,35))
    p1,p2=F(".03744"),F(".02331")
    d36=1-p1-p2-4*(1-p1)/35
    assert d36==1-p1-p2+F(4,35)*p1-F(4,35)==F(2902351,3500000)

    pr=common["minimum_order_weights"](winding,tri,x,y)
    high={q:F(1) for q in pr if q>=30}
    bound_expectation(high,F(1,10))
    assert all(pr[q]<arb(1)/40 for q in pr if q<30)
    assert all(pr[q]>arb(1)/40 for q in high if q<=40)
    assert sorted(q for q in high if q<=40)==[30,35,38,40]
    low=sum(map(F,(".0198450",".0067192",".0238946",".0115864",".0107808")))
    d105=1-low-F(1,10)
    assert d105==F(413587,500000)

    pr={q:F(n,26) for q,n in zip((1,2,30,32,33,36),(1,1,2,1,15,6),strict=True)}
    high={q:F(1) for q in (30,32,33,36)}
    bound_expectation({**high,2:F(4,17)},F(2,17))
    h=(1-2*pr[1])/34
    assert pr[1]==pr[2]>h and all(pr[q]>h for q in high)
    small=1-2*pr[1]-4*h
    assert small==sum(pr[q] for q in high)+F(4,17)*pr[2]-F(2,17)==F(180,221)
    print("36圈、105圈、小整数八位码作用的精确距离",d36,d105,small)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    prefix_margin(root)
    exact_examples(root)

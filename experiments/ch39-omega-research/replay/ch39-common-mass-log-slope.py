"""共同单调重质量在第 11431 枚上一圈交汇的最小斜率和对数凸质量跨度。"""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-native-center-parity-crossing.py")))
    ctx.prec = 256
    levels,_ = model["setup"]()
    loc = sorted([(v.log(),q) for q,row in levels.items() for v in row],key=lambda pair:float(pair[0]))
    for (left,q),(right,r) in zip(loc,loc[1:]):
        assert left < right or q==r
    gaps = [b[0]-a[0] for a,b in zip(loc,loc[1:])]
    groups = {q:[F(w)*6 for w in row.split()] for q,row in zip(model["charges"],
              ("1/15 1 2/3","8/15 0 1/3","2/15 0 1/3","1/5 1/3 0","2/5 0 0"))}
    groups[0] = [0,-7,F(-21,2)]
    groups[3] = [v+w for v,w in zip(groups[3],(F(1,10),F(1,6),0))]
    direction = (-325,603,-278)
    prefix = []
    total = [F(0)]*3
    for _,charge in loc:
        total = [v+w for v,w in zip(total,groups[charge])]
        prefix.append(total)
    projections = [dot(direction,row) for row in prefix]
    positive = [i for i,z in enumerate(projections[:-1]) if z>0]
    assert positive == [3,4,5,6,7,21,22,23,41] and projections[43]==0 and projections[-1]<0
    as_arb = model["q"]
    cost = sum(as_arb(projections[i])*gaps[i] for i in positive)
    L = (arb("1e13")/model["mz"]).log()
    em,weak = arb("127.930"),arb(".23129")
    A0 = [3*(1-weak)*em/5,weak*em,1/arb(".1180")]
    rhs = [2*arb.pi()*v-b*L for v,b in zip(A0,model["beta"])]
    assert dot(direction,model["beta"]) < 0
    demand = dot(direction,rhs)
    slope = demand/cost
    assert arb("106.7546704") < slope < arb("106.7546706")
    # 正投影区间饱和，一段零投影选取内部斜率，末端曝光置零。
    remaining = [v-slope*sum(as_arb(prefix[i][j])*gaps[i] for i in positive)
                 for j,v in enumerate(rhs)]
    pivot = 43
    mat = arb_mat([[as_arb(w),2*arb.pi()] for w in prefix[pivot][:2]])
    change,common = mat.solve(arb_mat([[r] for r in remaining[:2]])).entries()
    assert 0 < change < slope*gaps[pivot] and 16 < common < 17
    assert (as_arb(prefix[pivot][2])*change+2*arb.pi()*common-remaining[2]).contains(0)
    assert change+slope*sum(gaps[i] for i in positive)<L
    # 其它两种角权相对费米核的单边允许损失。
    vector = dot(direction,groups[0])
    higgs = dot(direction,[F(1,10),F(1,6),0])
    assert vector == -1302 and higgs == 68
    corrected = demand-18*as_arb(abs(vector))*arb(5)/42
    assert corrected/cost > arb("76.23229")
    # 凸的折线斜率由各结点处的非负增量给出。末段之后无正收益。
    ratios = [sum(as_arb(projections[i])*gaps[i] for i in positive if i>=j)
              /(loc[-1][0]-loc[j][0]) for j in range(max(positive)+1)]
    assert all(ratios[3]>r for i,r in enumerate(ratios) if i!=3)
    assert demand/ratios[3] > arb("115.64515") and corrected/ratios[3]>arb("82.58089")
    print("尖阈值最小质量对数斜率、取等逆耦合",slope.str(15),common.str(12))
    print("有限减法斜率下界", (corrected/cost).str(15))
    print("对数凸质量跨度的尖阈值、有限减法下界",
          (demand/ratios[3]).str(15),(corrected/ratios[3]).str(15))


if __name__ == "__main__":
    verify()

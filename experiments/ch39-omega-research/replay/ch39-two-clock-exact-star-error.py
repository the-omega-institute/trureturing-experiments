"""用黎曼零点钟角与整数背景相位，核验保留中心星时部门长边误差的精确值。"""
from itertools import permutations, product
from math import gcd
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb


class Angle:
    """区间值与二倍角指数的证书：两个整数钟频及 Q(sqrt(-867)) 中的方向。"""

    def __init__(self, value, clock724=0, clock804=0, x=1, y=0):
        self.value = arb(value)
        self.clock724, self.clock804 = clock724, clock804
        divisor = gcd(x,y)
        self.x, self.y = x//divisor, y//divisor

    def __add__(self, other):
        return Angle(self.value+other.value,self.clock724+other.clock724,self.clock804+other.clock804,
                     self.x*other.x-867*self.y*other.y,self.x*other.y+self.y*other.x)

    def __neg__(self):
        return Angle(-self.value,-self.clock724,-self.clock804,self.x,-self.y)

    def nonnegative(self):
        if self.value>0:
            return True
        if self.value<0:
            return False
        assert abs(self.value)<1
        assert self.clock724==self.clock804==self.y==0 and self.x>0
        return True


def verify(replay):
    shell = run_path(str(replay/"ch39-arc-relative-determinant.py"))
    sites = shell["sites"]
    assert shell["R"]==867
    hu,hd = [acb.zeta_zero(n).imag for n in (2345666,964890)]
    clocks = {d:(shell["theta"](d)-shell["theta"](244))*(-4*hu-2*hd) for d in (724,804)}
    delta = acb(0,clocks[724]-clocks[804]).exp().arg()
    assert arb("-2.75")<delta<arb("-2.74")
    radius = Angle(-delta/2,-1,1)
    turn = Angle(2*arb.pi())
    pole = (17,17,17)
    assert sites[14]==pole
    edges = {}
    clock_edges = {}
    for i,j,d,z in shell["edges"]:
        if d<=372:
            continue
        if 14 in (i,j):
            t = Angle(0)
        else:
            phase = z**-6 * acb(0,clocks.get(d,0)).exp()
            # e^(2 i t) 的背景方向是原有理方向共轭的十二次幂。
            a = 867*(867+shell["dot"](pole,sites[i])+shell["dot"](sites[i],sites[j])+shell["dot"](sites[j],pole))
            b = shell["det"](pole,sites[i],sites[j])
            raw = Angle(0,x=a,y=b)
            two = raw+raw
            four = two+two
            twelve = four+four+four
            t = Angle(phase.arg(),2*int(d==724),2*int(d==804),twelve.x,twelve.y)
        edges[i,j] = t
        if d in clocks:
            clock_edges[i,j] = t

    # 在目标半径处以中心星固定外点提升，参见有理双边界证书。
    lifts = []
    for j in (6,7,8):
        keys = [e for e in clock_edges if e[1]==j]
        options = []
        for ks in product((-1,0,1),repeat=3):
            targets = {e:clock_edges[e]+Angle(k*turn.value) for e,k in zip(keys,(0,*ks),strict=True)}
            if not any(x.value-y.value>4*radius.value for x,y in permutations(targets.values(),2)):
                options.append(targets)
        lifts.append(options)

    nodes = (6,7,8,14)
    cycles = [seq for n in (2,3,4) for seq in permutations(nodes,n) if min(seq)==seq[0]]
    for patterns in product(*lifts):
        chosen = {key:v for group in patterns for key,v in group.items()}
        costs = {}
        for i in range(6):
            spokes = {j:t for (row,j),t in chosen.items() if row==i}
            spokes[14] = Angle(0)
            for j,k in permutations(spokes,2):
                cost = spokes[k]+(-spokes[j])+radius+radius
                if (j,k) not in costs or cost.value<costs[j,k].value:
                    costs[j,k] = cost
        # 任一严格较小误差会严格减小每个环的总费用。
        assert any((-sum((costs[j,k] for j,k in zip(seq,seq[1:]+seq[:1],strict=True)),Angle(0))).nonnegative()
                   for seq in cycles)

    # 有理近似只选择整数绕数，实际满足误差的顶点势由有向路径精确给出。
    units = run_path(str(replay/"ch39-two-clock-department-star-error.py"))["PHASE_TRILLIONTHS"]
    arcs = []
    for (i,j),t in edges.items():
        guess = arb(units[j]-units[i])/10**12
        k = ((guess-t.value)/turn.value+arb(1)/2).floor().unique_fmpz()
        assert k is not None
        lifted = t+Angle(int(k)*turn.value)
        arcs.extend(((i,j,lifted+radius),(j,i,(-lifted)+radius)))
    potentials = [Angle(0) for _ in range(19)]
    for _ in range(18):
        for i,j,cost in arcs:
            new = potentials[i]+cost
            if new.value < potentials[j].value:
                potentials[j] = new
    for i,j,cost in arcs:
        assert (potentials[i]+cost+(-potentials[j])).nonnegative()
    print("保留中心星的精确长边误差等于两部门有效钟角差的圆周距离的一半",radius.value.str(24))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

"""已发表最优字概率对公平自适应停时的均值、凸序与逐步位势下界。"""

from math import comb
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def target_masses(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    optimal = load("realized-vacuum-word-masses")
    target = load("height-independent-flavor-vacuum")
    winding = load("mincut-global-winding")
    geometry = load("mincut-real-triangle-tomography")
    local = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = target["minimum_order_weights"](winding,winding["terms"](shell),x,y)
    theta = geometry["unit"]((108469,109)).arg()
    inverse = arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)]).inv()
    sizes = [1]+[2]*106
    base = [inverse[j,0]-arb(n)/213 for j,n in enumerate(sizes)]
    free = [inverse[j,106] for j in range(107)]
    active = [(q,h*weights[q].sqrt()) for q,h in zip(optimal["SUPPORT"],optimal["WINNER"])]
    direction = [sum(v*inverse[j,q] for q,v in active) for j in range(107)]
    i,j = 3,14
    det = direction[i]*free[j]-direction[j]*free[i]
    assert not det.contains(0)
    scale = (base[j]*free[i]-base[i]*free[j])/det
    tail = (base[i]*direction[j]-base[j]*direction[i])/det
    groups = [inverse[j,0]+scale*direction[j]+tail*free[j] for j in range(107)]
    assert (sum(groups)-1).contains(0) and min(groups)>0
    return groups


def fixed_time_calls(n):
    # 按奇偶性求 |S_n| 的全部整数二项式计数，再用尾概率递推折线期望。
    counts = [0]*(n+1)
    k = n % 2
    numerator = comb(n,(n+k)//2)
    for s in range(k,n+1,2):
        counts[s] = (1 if s==0 else 2)*numerator
        numerator = numerator*(n-s)//(n+s+2)
    assert sum(counts)==2**n
    positive_tail = 2**n-counts[0]
    call = sum(s*c for s,c in enumerate(counts))
    calls = []
    for s in range(107):
        calls.append(arb(call)/(2**n))
        call -= positive_tail
        positive_tail -= counts[s+1]
    return calls


def verify(root):
    groups = target_masses(root)
    mean_square = sum(j*j*p for j,p in enumerate(groups))
    assert arb("3760.86449156219814840")<mean_square<arb("3760.86449156219814842")
    print("最优字概率的二阶矩及有限均值公平停时的期望",mean_square.str(21))
    calls = [sum((j-k)*groups[j] for j in range(k+1,107)) for k in range(107)]
    print("目标 E|J|",calls[0].str(18))
    for m in (4430,4431):
        gaps = [v-u for v,u in zip(fixed_time_calls(m),calls)]
        least = min(g.lower() for g in gaps)
        if m==4430:
            assert gaps[0]<arb("-.00694837")
        else:
            assert least>arb(".00032071")
        print("固定步数、k=0 处的凸函数期望差、全部 k 的最小差",m,gaps[0].str(17),least.str(17))

    # 已由鞅不等式得到的必要位势递推，在对称半轴上做逐步区间包围。
    probability = [v if j==0 else v/2 for j,v in enumerate(groups)]
    zero = arb(0)
    gap = calls
    for n in range(1,4457):
        neighbor = [gap[1]]+[(a+b)/2 for a,b in zip(gap[:-2],gap[2:])]
        gap = [(v-p).max(zero) for v,p in zip(neighbor,probability)]+[zero]
        if n==4455:
            assert all(v.is_zero() for j,v in enumerate(gap) if j!=1)
            assert arb(".0021969351")<gap[1]<arb(".0021969352")
            print("4455 步后在 ±1 各自的位势缺口下界",gap[1].str(19))
            assert 2*gap[1]>arb(".00439387")
    assert all(v.is_zero() for v in gap)
    print("位势递推全部归零的首个步数",4456)


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)

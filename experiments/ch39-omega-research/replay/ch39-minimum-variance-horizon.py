"""最小平均步数字概率的精确最坏停时，含四个零质量位置的正向检查。"""

from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def terminal(root):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    minima = run_path(str(Path(__file__).with_name("ch39-minimum-variance-fixed-weights.py")))
    support = load("realized-vacuum-word-masses")["SUPPORT"]
    model = load("mincut-global-winding")
    local = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 320
    theta = load("mincut-real-triangle-tomography")["unit"]((108469,109)).arg()
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = load("height-independent-flavor-vacuum")["minimum_order_weights"](
        model, model["terms"](shell), x, y)
    signs = dict(zip(support, minima["OPTIMAL_SIGNS"]))
    rhs = arb_mat([[int(q==0), signs[q]*weights[q].sqrt() if q in signs else 0, int(q==106)]
                   for q in range(107)])
    basis = arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)]).solve(rhs)
    det = basis[17,1]*basis[60,2]-basis[60,1]*basis[17,2]
    s = (basis[60,0]*basis[17,2]-basis[17,0]*basis[60,2])/det
    t = (basis[17,0]*basis[60,1]-basis[60,0]*basis[17,1])/det
    mass = [(basis[j,0]+s*basis[j,1]+t*basis[j,2])/(1 if j==0 else 2) for j in range(107)]
    for j in (17,60):
        assert mass[j].contains(0)
        mass[j] = arb(0)
    assert min(m for j,m in enumerate(mass) if j not in (17,60))>0
    assert (mass[0]+2*sum(mass[1:])-1).contains(0)
    return mass


def verify(root):
    target = terminal(root)
    zero = arb(0)
    gap = [2*sum(((r-j)*target[r] for r in range(j+1,107)),zero) for j in range(107)]
    active = [arb(1)]+[zero]*106
    held = [zero]*107
    crossed = {17:False, 60:False}
    for time in range(3951):
        following = [zero]*107
        released = [zero]*107
        for j,d in enumerate(gap):
            if j in crossed:
                # 目标质量精确为零。即使恰好触零，也必须把全部活动质量放行。
                assert held[j].is_zero()
                following[j] = (gap[j-1]+gap[j+1])/2
                released[j] = active[j]
                if active[j]>0:
                    crossed[j] = True
                continue
            if d.is_zero():
                held[j] += active[j]
                continue
            assert held[j].is_zero() and j<106
            unconstrained = (gap[abs(j-1)]+gap[j+1])/2-target[j]
            if unconstrained>0:
                following[j] = unconstrained
                released[j] = active[j]
            else:
                assert unconstrained<0 and active[j]>d
                released[j] = d
                held[j] += active[j]-d
        if time==3950:
            assert all(d.is_zero() for d in gap[1:]) and gap[0]>0
            assert arb(".00121105912") < gap[0] < arb(".00121105913")
            assert (released[0]-gap[0]).contains(0) and all(r.is_zero() for r in released[1:])
            print("3950 步时仅在原点为正的剩余势差、贪心策略的末步概率",gap[0].str(23))
        assert released[106].is_zero()
        active = [released[1]]+[(a+b)/2 for a,b in zip(released[:-2],released[2:])]+[released[105]/2]
        gap = following
    assert all(d.is_zero() for d in gap) and all(crossed.values())
    assert all(held[j].is_zero() and active[j].is_zero() for j in (17,60))
    errors = [h+a-p for h,a,p in zip(held,active,target)]
    assert all(e.contains(0) and abs(e)<arb("1e-35") for e in errors)
    print("末步后所有位置的势差精确为零",3951)
    print("正向终端逐字概率的区间残差绝对上界",max(abs(e).upper() for e in errors).str(10))
    print("允许经过且从未停止的指数",(-60,-17,17,60))


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)

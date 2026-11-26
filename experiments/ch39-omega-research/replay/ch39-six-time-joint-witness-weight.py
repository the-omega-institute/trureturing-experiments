"""指定联合见证的六时间点权界、尾概率及实际支撑的内界。"""
from math import prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, arb_poly, ctx


root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
bands = run_path(str(root/"ch39-six-time-height-bands.py"))
clock = run_path(str(root/"ch39-arithmetic-height-bandwidth.py"))
dictionary = run_path(str(root/"ch39-seventh-time-zero-dictionary.py"))
physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
first, number = 8359, 2306888100085
scale = bands["scale"]


def coordinate(h):
    return 1/(1+(h/scale)**2)


def moment_bounds(height, endpoint, denominator, moments):
    c = coordinate(height)
    rescaled = denominator(arb_poly([0, endpoint]))
    moment = [value/endpoint**i for i,value in enumerate(moments)]
    nodes, weights = bands["cutoff_rule"](c/endpoint, moment, rescaled)
    for upper in (False, True):
        bands["verify_indicator"](nodes, rescaled, upper)
    tail_lower = sum(w for n,w in zip(nodes,weights) if n < c/endpoint)
    tail_upper = tail_lower+weights[0]
    local = arb_mat([[endpoint*moments[i+j+1]-moments[i+j+2] for j in range(3)]
                     for i in range(3)])
    v = arb_mat([[1],[c],[c*c]])
    cap = denominator(c)/(c*(endpoint-c)*(v.transpose()*local.solve(v))[0,0])
    assert cap.overlaps(weights[0])
    print("高度、闭尾段上下界、局部 Christoffel 界", height.str(17), tail_lower.str(22),
          tail_upper.str(22),cap.str(22),flush=True)
    return tail_lower,tail_upper,cap


def witness_rule(numbers, heights, times, moments, denominator):
    coordinates = {n:coordinate(heights[n]) for n in numbers}
    coordinates[None] = arb(0)
    _, values = dictionary["certified_rule"]([*numbers,None],coordinates,times,moments,denominator)
    # 显式粘合区间：任取高度对应末节点 z≤10^-28 的实际合格零点。
    z=arb(0).union(arb("1e-28"))
    nodes=[coordinates[n] for n in numbers]+[z]
    y=arb_poly([0,1])
    weights=[]
    for i,v in enumerate(nodes):
        others=[t for j,t in enumerate(nodes) if i!=j]
        numerator=prod(y-t for t in others)
        expected=sum(a*m for a,m in zip(numerator,moments))
        weights.append(denominator(v)*expected/prod(v-t for t in others))
    assert all(w>0 for w in weights)
    return values[-2],values[-1]+values[-2]


def verify():
    ctx.prec = 300
    heights = {n:acb.zeta_zero(n).imag for n in (first,number)}
    ctx.prec = 640
    times = clock["times"]()[:7]
    times[6] += arb("1.360900010e-8")
    denominator,moments = bands["radial"]["compact_moments"](times)
    upper = coordinate(heights[first])
    low,high,_ = moment_bounds(arb("5.9e11"),upper,denominator,moments)
    assert low>arb(".39132009747930899") and high<arb(".39156811030917581")
    _,joint_high,cap = moment_bounds(heights[number],upper,denominator,moments)
    assert cap<arb(".000253654097221442") and joint_high<arb(".39156809996557095")
    assert cap/low < arb(".000648202") and low-cap>arb(".39106644338")
    # 初始支撑的四枚均已通过原分界9三窗，分别围住两处高斯内节点。
    support=[first,142268163595,142271930613,295869382715,295872524995,number]
    known=set(dictionary["lower_support"])|set(dictionary["upper_support"])|set(dictionary["spectrum_lower"])
    assert set(support[:-1])<=known
    ctx.prec = 192
    _, alpha, curve = physical["observables"](number)
    assert arb(".117")<alpha<arb(".119") and arb(".998")<curve<arb("1.002")
    for n in support:
        if n not in heights:
            heights[n]=acb.zeta_zero(n).imag
    ctx.prec=512
    atom,tail = witness_rule(support,heights,times,moments,denominator)
    assert atom>arb(".00025345503694") and tail>arb(".39156808912862")
    print("实际零点趋近构造的联合见证点权与尾权",atom.str(21),tail.str(21),flush=True)


if __name__=="__main__":
    verify()

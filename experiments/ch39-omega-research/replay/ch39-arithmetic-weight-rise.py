"""九枚物理基点同时命中前两个时间读数所需的相邻正权增长。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


joint = run_path(str(Path(__file__).with_name("ch39-joint-arithmetic-support.py")))
numbers = (8359,11661,11934,12207,22027,28325,48067,58214,95478)


def prefixes(s, nodes, ratio):
    weights, aa, bb = arb(0), arb(0), arb(0)
    pairs = []
    for i,u in enumerate(nodes):
        z = s*u
        a = z/(1+z)
        b = 3*z**2/((1+z)*(1+4*z))
        w = ratio**i
        weights += w
        aa += w*a
        bb += w*b
        pairs.append((aa/weights,bb/weights))
    return pairs


def at_first_reading(nodes, weights, target):
    lo,hi = arb(0),arb(1)
    for _ in range(100):
        mid = (lo+hi)/2
        read = sum(w*(mid*u)/(1+mid*u) for w,u in zip(weights,nodes))
        if read < target:
            lo = mid
        elif read > target:
            hi = mid
        else:
            break
    s = lo.union(hi)
    return sum(w*3*(s*u)**2/((1+s*u)*(1+4*s*u)) for w,u in zip(weights,nodes))


def exclude(nodes, first, second, ratio):
    # 选择分离直线的斜率用中点，验算对整段取外端点。
    def cover(lo, hi, depth):
        left, right = prefixes(lo,nodes,ratio), prefixes(hi,nodes,ratio)
        if max(a.upper() for a,_ in right) < first or min(a.lower() for a,_ in left) > first:
            return 1
        mid = (lo+hi)/2
        points = prefixes(mid,nodes,ratio)
        lower, upper = [], []
        for a,b in points:
            if a > first:
                lower.append(float((b-second)/(a-first)))
            elif a < first:
                upper.append(float((b-second)/(a-first)))
        ls = max(lower,default=-1.)
        us = min(upper,default=1.)
        if ls < us:
            slope = arb((ls+us)/2)
            limiting = left if slope >= 0 else right
            if all(b-second - slope*(a-first) < 0 for (_,b),(a,_) in zip(right,limiting)):
                return 1
        assert depth < 55
        return cover(lo,mid,depth+1) + cover(mid,hi,depth+1)
    odds = first/(1-first)
    return cover((odds/nodes[0]).lower(), (odds/nodes[-1]).upper(),0)


def verify():
    ctx.prec,ctx.cap = 1536,3
    times = joint["clock"]["times"]()
    first = 1-times[1]
    second = first-(1-times[2])/4
    h = [acb.zeta_zero(n).imag for n in numbers]
    nodes = [(h[0]/v)**2 for v in h]
    lo,hi = arb("1.14827"),arb("1.14828")
    for index,r in enumerate((lo,hi)):
        weights = [r**i for i in range(len(nodes))]
        read = at_first_reading(nodes,[w/sum(weights) for w in weights],first)
        assert (read < second if index == 0 else read > second)
    count = exclude(nodes,first,second,lo)
    print("九枚上最大相邻权重倍率的下确界在 (1.14827,1.14828)，尺度分段",count)
    arith = at_first_reading(nodes,[u/sum(nodes) for u in nodes],first)
    gap = 4*(second-arith)
    assert arb(".00003635") < gap < arb(".00003636")
    print("限制九枚并按逆高度平方赋权，校准首项后的第二项亏额",gap.str(12))


if __name__ == "__main__":
    verify()

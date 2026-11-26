"""有限十九枚局部调和权的成对高度矩、有限观察窗和占据对角。"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def prefix_ratio():
    # 和上游使用同一个起点，移到加权均值后只需求两次单变量矩。
    # 按递增起点取分块，避免在远高于初始高度的区间继续采用低高度求零点例程。
    boundaries = (1, 32, 128, 512, 2048, 8192, 11433)
    zeros = [z for lo, hi in zip(boundaries, boundaries[1:]) for z in acb.zeta_zeros(lo, hi-lo)]
    heights = [z.imag for z in zeros]
    assert len(heights) == 11432 and all(a < b for a,b in zip(heights,heights[1:]))
    heights = heights[:-1]
    raw = [h/(arb(1)/4+h*h) for h in heights]
    total = sum(raw)
    mean = sum(w*h for w,h in zip(raw,heights))/total
    variance = sum(w*(h-mean)**2 for w,h in zip(raw,heights))/total
    fourth = sum(w*(h-mean)**4 for w,h in zip(raw,heights))/total
    ratio = arb(3)/2 + fourth/(2*variance**2)
    print("第一枚至第 11431 枚的调和成对差矩比", ratio.str(16))
    return ratio


def verify():
    transform = run_path(str(Path(__file__).with_name("ch39-finite-fejer-observation.py")))["finite_transform"]
    ctx.prec = 128
    assert prefix_ratio() > 2
    ctx.prec = 160
    ordinates = [zero.imag for zero in acb.zeta_zeros(11421, 21)]
    assert all(a < b for a, b in zip(ordinates, ordinates[1:]))
    heights = ordinates[1:-1]
    raw = [h/(arb(1)/4 + h*h) for h in heights]
    weights = [x/sum(raw) for x in raw]
    pairs = list(combinations_with_replacement(range(len(heights)), 2))
    def moment(k):
        return sum(2*weights[i]*weights[j]*(heights[j]-heights[i])**k for i,j in pairs if i<j)
    ratio = moment(4)/moment(2)**2
    assert 2 < ratio < 3
    print("十九枚局部调和权的有序差矩比", ratio.str(16))
    masses = [weights[i]*weights[j]*(1 if i==j else 2) for i,j in pairs]
    diagonal = sum(w*w for w in masses)
    centers = [heights[i]+heights[j] for i,j in pairs]
    gaps = [(abs(centers[i]-centers[j]), 2*masses[i]*masses[j])
            for i,j in combinations(range(len(centers)), 2)]
    eta = arb(1)/1000
    assert sum(int(g < eta) for g,w in gaps) == 1
    assert all(g < eta or g > eta for g,w in gaps)
    excess = sum(w*(1-g/eta) for g,w in gaps if g < eta)
    sufficient = 64*(diagonal+excess)/(eta*excess)
    print("整线碰撞超额、对角及通用界所需半径阈值", excess.str(13), diagonal.str(13), sufficient.str(16))
    def margin(radius):
        return (diagonal*(transform(arb(0),eta,radius)-1)
                + sum(w*transform(g,eta,radius) for g,w in gaps))
    # 积分在半径上单调，两个相邻整数即可核验首次越过。
    low, high = margin(111531), margin(111532)
    assert low < 0 < high
    assert arb("2.43290996641730") < ratio < arb("2.43290996641732")
    assert 11331842 < sufficient < 11331843
    print("有限积分首次超过对角的整数半径", 111532, low.str(12), high.str(12))


if __name__ == "__main__":
    verify()

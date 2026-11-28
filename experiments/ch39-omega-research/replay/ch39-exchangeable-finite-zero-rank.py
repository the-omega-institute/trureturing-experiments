"""四个等权子壳类组成的有理正混合，在给定的两组实际编号中共同选中中央。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_poly, ctx


# 键依次给出上下三站轨道数和极点标志、交叠数、上/下边界及概率分子。
SUPPORT = (
    ((1,False,1,False,1), "D","D", 1),
    ((1,True,1,False,0), "D","I", 119393),
    ((2,False,1,False,0), "D","I", 743534),
    ((2,False,2,False,1), "D","I", 137072),
)
CENTRAL = {-4:2345666, 2:964890}


def verify(replay, tile=None, *, support=SUPPORT, denominator=10**6,
           holdouts=((-4,483437,200),(2,33560,10000)), closest=(-4,11757),
           gap_bounds=("5358.20","5358.21")):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    reference = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    ctx.prec = 160
    triple = reference["geometry"](shell,pump)
    classes = grouped["geometry_classes"]()
    weights = sum(case[-1] for case in support)
    assert weights == denominator
    polynomial = arb_poly([0])
    for key, up, down, numerator in support:
        group = classes[key]
        term = sum((reference["squared_factors"](triple(-4,a,up)+triple(2,b,down))
                    for a,b in group), arb_poly([0]))
        polynomial += term * (arb(numerator)/(weights*len(group)))

    adj = pump["adjoint"]
    frames = list(pump["rotation_bases"](shell))
    soft = shell["soft"]
    pole = soft["pole"]
    radial = acb_mat(19,19)
    for j in soft["star"]["orbit"]:
        radial[pole,pole] += 1
        radial[j,j] += 1
        radial[pole,j] -= 1
        radial[j,pole] -= 1
    blocks = {}
    for q in CENTRAL:
        base, links = pump["clock_matrices"](shell,q)
        blocks[q] = [[adj(b)*op*b for op in (base,*links,radial)] for b in frames]
    def rate(q,height):
        phases = shell["geo"]["phases"](height)
        z0,z1 = [phases[k]**q for k in (724,804)]
        loss = arb(23)/25*arb(2)**-20/(height**2+soft["point"]["study"]["offset"])
        levels = []
        for base,edge0,edge1,r in blocks[q]:
            links = edge0*z0+edge1*z1
            levels.append(reference["minimum"](base+links+adj(links)-loss*r))
        return 9**12 * sum(polynomial(energy/9-1) for energy in levels)
    central = {q:rate(q,acb.zeta_zero(n).imag) for q,n in CENTRAL.items()}
    # 两枚训练集合外的实际编号逐荷反转同一有理权重的排序。
    for charge, index, lower in holdouts:
        improvement = central[charge]-rate(charge,acb.zeta_zero(index).imag)
        assert improvement > lower
        print("集合外的荷、编号、目标减竞争作用",charge,index,improvement.str(13))
    # 先横跨相邻零点的闭区段，再遍历正零点开头的连续编号。
    intervals = [(q,n-320,n+321) for q,n in CENTRAL.items()]
    intervals += [(0,start,start+512) for start in range(1,32769,512)]
    if tile is not None:
        assert 0 <= tile < len(intervals)
        intervals = intervals[tile:tile+1]
    best = arb.pos_inf()
    nearest = None
    count = 0
    for q,start,stop in intervals:
        for n,zero in enumerate(acb.zeta_zeros(start,stop-start),start):
            for charge in (q,) if q else CENTRAL:
                if n != CENTRAL[charge]:
                    enclosure = rate(charge,zero.imag)-central[charge]
                    gap = enclosure.lower()
                    if (charge,n) == closest:
                        assert enclosure.upper() < arb(gap_bounds[1])
                    else:
                        assert gap > arb(gap_bounds[1]), (charge,n,gap)
                    if gap < best:
                        best, nearest = gap, (charge, n)
                    count += 1
    assert best > arb(gap_bounds[0])
    if tile is None:
        assert count == 66816 and nearest == closest
    print("分区、严格比较数、最小区间下界与所在编号",tile,count,best.str(14),nearest)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent,
           int(sys.argv[2]) if len(sys.argv)>2 else None)

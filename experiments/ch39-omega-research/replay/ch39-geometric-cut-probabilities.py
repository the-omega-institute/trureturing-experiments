"""按子壳的边界边数赋 Gibbs 权时，两枚实际零点排除整族正参数。"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from math import comb, prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def cut_classes(geo, masks):
    sites = geo["sites"]
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    orbits = [sorted({i,rot[i],rot[rot[i]]}) for i in range(len(sites))
              if i == min(i,rot[i],rot[rot[i]])]
    classes = defaultdict(list)
    for mask in masks:
        keep = {i for b,orbit in enumerate(orbits) if mask & (1 << b) for i in orbit}
        cut = sum((i in keep) != (j in keep) for i,j,_,_ in geo["edges"])
        assert cut % 6 == 0
        classes[cut//6].append(mask)
    return dict(sorted(classes.items()))


def bernstein_weights(classes):
    # 用 x 或 1/x 覆盖整个正半轴。在 [a/8,(a+1)/8] 上将每项幂
    # 写成 13 次 Bernstein 多项式，按每类计数加权，再将每行归一。
    degree = max(classes)
    rows = []
    for reverse, a, j in product((False,True),range(8),range(degree+1)):
        row = []
        for k, masks in classes.items():
            exponent = degree-k if reverse else k
            value = sum(Fraction(comb(exponent,h)*a**(exponent-h)*comb(j,h),comb(degree,h))
                        for h in range(min(exponent,j)+1))
            row.append(len(masks)*value/8**exponent)
        total = sum(row)
        assert total > 0
        rows.append([arb(v.numerator)/v.denominator for v in (t/total for t in row)])
    return arb_mat(rows)


def verify():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay/f"ch39-{name}.py"))

    geom = load("geometric-selector-mixture")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    ctx.prec = 320
    triple = geom["geometry"](shell,pump)
    masks = load("geometric-exchangeable-priors")["MASKS"]
    classes = cut_classes(shell["geo"],masks)
    assert len(classes) == 11 and max(classes) == 13 and sum(map(len,classes.values())) == 126
    weights = bernstein_weights(classes)
    # 公用一次抽样或上下独立抽样。双方边界的四种选择与掩码无关。
    rivals = ((-4,2345666,113,4), (2,964890,1,3))
    def moments(charge, index, modes):
        energies = light(charge,index)
        factors = [{m: [prod((e-r)**2 for r in triple(q,m,mode)) for e in energies] for m in masks}
                   for q,mode in zip((-4,2),modes,strict=True)]
        averaged = [arb_mat([[sum(f[m][j] for m in block)/len(block) for j in range(3)]
                             for block in classes.values()]) for f in factors]
        independent = averaged[0]*averaged[1].transpose()
        same = arb_mat([[sum(factors[0][m][j]*factors[1][m][j] for m in block for j in range(3))
                         /len(block)] for block in classes.values()])
        return independent, same

    lower = None
    for modes in product(("D","I"),repeat=2):
        paired = arb_mat(len(classes),len(classes))
        shared = arb_mat(len(classes),1)
        for q, center, rival, amount in rivals:
            cp,cs = moments(q,center,modes)
            rp,rs = moments(q,rival,modes)
            paired += amount*(cp-rp)
            shared += amount*(cs-rs)
        # 每个比值是其 Bernstein 分子/分母同位置系数比的凸组合。
        # 参数独立时取两维张量积；共享掩码时只需一维。
        bounds = [weights*paired*weights.transpose()/7, weights*shared/7]
        value = min(v.lower() for matrix in bounds for v in matrix.entries())
        assert value > 32398
        lower = value if lower is None else min(lower,value)
    print("截边数", [6*k for k in classes])
    print("两个实际竞争者按 4:3 加权的统一差额下界", lower.str(12))


if __name__ == "__main__":
    verify()

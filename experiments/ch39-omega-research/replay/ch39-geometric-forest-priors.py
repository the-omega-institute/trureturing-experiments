"""实际零点上型 33 与下型 70 排除普通图森林系数生成的两类子集抽样。

对每个掩码 M，在普通无磁拉普拉斯的原度主块或诱导图上取 det(I+zL)
的一项系数。对同一边界、同一次数在 126 个掩码上的非零总和归一。
上下型可各取一行而独立抽样，也可将两行逐掩码相乘、归一后共用掩码。
磁参考谱的两种边界另选，允许与所选森林规则关联，给定规则时与抽出的掩码独立。
"""
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpz_mat


def forest_rows(geo, masks):
    sites = geo["sites"]
    rotation = [sites.index((y, z, x)) for x, y, z in sites]
    orbits = [sorted({i, rotation[i], rotation[rotation[i]]}) for i in range(len(sites))
              if i == min(i, rotation[i], rotation[rotation[i]])]
    assert list(map(len, orbits)) == [3]*6 + [1]
    size = len(sites)
    full = fmpz_mat(size, size)
    for i, j, _, _ in geo["edges"]:
        full[i, i] += 1
        full[j, j] += 1
        full[i, j] -= 1
        full[j, i] -= 1

    totals = {(mode, degree): [] for mode in ("D", "I") for degree in range(size+1)}
    for mask in masks:
        keep = sorted(i for j, orbit in enumerate(orbits) if mask & (1 << j) for i in orbit)
        n = len(keep)
        principal = fmpz_mat([[full[i, j] for j in keep] for i in keep])
        induced = fmpz_mat(principal)
        for i in range(n):
            induced[i, i] = -sum(principal[i, j] for j in range(n) if j != i)
        for mode, block in (("D", principal), ("I", induced)):
            # det(I+zL) 的 k 次系数等于特征多项式的 n-k 次系数乘 (-1)^k。
            polynomial = block.charpoly()
            for k in range(size+1):
                value = int(polynomial[n-k]) * (-1)**k if k <= n else 0
                assert value >= 0
                totals[mode, k].append(value)

    rows = [weights for weights in totals.values() if sum(weights) > 0]
    assert len(rows) == 38
    return rows


def normalize(weights):
    total = sum(weights)
    assert total > 0
    return [arb(value)/total for value in weights]


def action(charge, index, modes, light, triple, masks, independent, shared):
    energies = light(charge, index)
    factors = []
    for q, mode in zip((-4, 2), modes, strict=True):
        rows = []
        for mask in masks:
            roots = triple(q, mask, mode)
            rows.append([prod((e-r)**2 for r in roots) for e in energies])
        factors.append(arb_mat(rows))

    # 每个字符先平均两边因子再相乘，给出独立分布的 38×38 个期望。
    up, down = factors
    separate = (independent*up) * (independent*down).transpose()
    # 共用掩码时先逐掩码计算完整作用，然后按两森林系数的归一化乘积平均。
    diagonal = arb_mat([[sum(up[i, k]*down[i, k] for k in range(len(energies)))]
                        for i in range(len(masks))])
    return separate, shared*diagonal


def verify():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    reference = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    ctx.prec = 320
    triple = reference["geometry"](shell, pump)
    masks = grouped["MASKS"]
    assert len(masks) == 126
    rows = forest_rows(shell["geo"], masks)
    independent = arb_mat([normalize(row) for row in rows])
    shared = arb_mat([normalize([u*d for u, d in zip(up, down, strict=True)])
                      for up, down in product(rows, repeat=2)])
    # k=0 时两种普通图边界的行相同，计数按带边界和次数的规则标签。
    assert shared.nrows() == len(rows)**2

    lower = []
    count = 0
    for modes in grouped["MODE_PAIRS"]:
        paired = arb_mat(len(rows), len(rows))
        common = arb_mat(shared.nrows(), 1)
        for q, rival, weight in ((-4, 33, 6), (2, 70, 7)):
            target = grouped["TARGET"][q]
            ap, ac = action(q, target, modes, light, triple, masks, independent, shared)
            bp, bc = action(q, rival, modes, light, triple, masks, independent, shared)
            paired += weight*(ap-bp)
            common += weight*(ac-bc)
        bounds = [min(v.lower() for v in table.entries())/13 for table in (paired, common)]
        assert all(b > 4413 for b in bounds)
        lower.append(bounds)
        count += paired.nrows()*paired.ncols() + common.nrows()
    assert count == 11552
    print("非零森林系数规则与带标号的抽样、磁边界组合数", len(rows), count)
    print("上型 33 和下型 70 按 6:7 加权的目标减竞争作用，独立及共用掩码的统一下界",
          *(min(column).str(14) for column in zip(*lower, strict=True)))


if __name__ == "__main__":
    verify()

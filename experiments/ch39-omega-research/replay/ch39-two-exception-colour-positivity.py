"""两个独立插值例外下，牛顿与真空目标强制色逆耦合为负。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

ROOT = Path(__file__).resolve().parent
PARTS = 12
UNIT = 2**16


def lower(value):
    return int((value * UNIT).lower().floor().unique_fmpz())


def upper(value):
    return int((value * UNIT).upper().ceil().unique_fmpz())


def verify():
    previous = run_path(str(ROOT / "ch39-two-copy-gram-cutoff-two-exceptions.py"))
    edges, rows, heat, _ = previous["model_inputs"]()
    cut = arb(1)/625
    pi = arb.pi()
    k0 = sum(alpha for _, alpha, _ in rows) / (2*pi)
    assert k0 < 0
    # 对所有例外选取求牛顿主项的最小一阶矩，热核误差与选取无关。
    lows = [w*a for a, _, w in edges]
    highs = [w*b for _, b, w in edges]
    total_low, total_high = sum(lows), sum(highs)
    k_error = sum((-(w*alpha)).max(0) * (cut*b) * (1-heat(2, z/(cut*b)))
                  for _, b, w in edges for _, alpha, z in rows) / (2*pi)
    k_max = arb(0)
    for i, j in combinations(range(len(edges)), 2):
        p = (total_low-lows[i]-lows[j]).min(total_high-highs[i]-highs[j])
        p += lows[i].min(highs[i])+lows[j].min(highs[j])
        k_max = k_max.max((k0*cut*p).max(0)+k_error)
    assert k_max < 1

    lib = run_path(str(ROOT / "ch39-chiral-minimal-neutral-gravity.py"))
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = run_path(str(ROOT / "ch39-signed-gravity-unification.py"))
    heavy = [(row[2], z) for row, (_, _, z) in zip(gauge["gauge_weights"](model), rows[1:], strict=True)]
    logscale = (arb(10**16)/model["mz"]).log()
    offset = initial[2]-model["beta"][2]*logscale/(2*pi)
    offset += sum(b*(arb.const_euler()+z.log()) for b, z in heavy)/(4*pi)
    log_coeff = -sum(b for b, _ in heavy)/(4*pi)
    w0 = -sum(g for g, _, _ in rows)/(8*pi)
    assert log_coeff > 0 and w0 > 0
    spin_pos = sum(max(g, 0)*z for g, _, z in rows)/(4*pi)
    spin_neg = sum(max(-g, 0)*z for g, _, z in rows)/(4*pi)

    # 每个数组元素都是一个谱线对在整段插值上的向外取整界。
    low_w = np.empty((len(edges), PARTS), dtype=np.int64)
    high_w = np.empty_like(low_w)
    high_colour = np.empty_like(low_w)
    for i, (a, b, weight) in enumerate(edges):
        pos, neg = weight.max(0), (-weight).max(0)
        for k in range(PARTS):
            left = ((PARTS-k)*a+k*b)/PARTS
            right = ((PARTS-k-1)*a+(k+1)*b)/PARTS
            w_left, w_right = w0*weight*left**2, w0*weight*right**2
            low_w[i,k] = lower(w_left.min(w_right)-(pos*spin_neg+neg*spin_pos)*right/cut)
            high_w[i,k] = upper(w_left.max(w_right)+(pos*spin_pos+neg*spin_neg)*right/cut)
            remainder = arb(0)
            for bj, z in heavy:
                y = z/(cut*left)
                ein_upper = y.min(1+y.log().max(0))
                remainder += (-(weight*bj)).max(0)*ein_upper
            log_term = (weight*left.log()).max(weight*right.log())
            high_colour[i,k] = upper(log_coeff*log_term + remainder/(4*pi))

    bound0 = upper(offset)
    # 每个整数中间结果的绝对值小于下面的保守估计，保证 int64 运算精确。
    largest = max(abs(int(x)) for table in (low_w, high_w, high_colour) for x in table.flat)
    assert (len(edges)+4)*largest+abs(bound0) < 2**62
    assert arb("1e-122")/cut**2 < arb(1)/UNIT

    def box_sums(table, total, i, j):
        shared = total-table[i]-table[j]
        return shared[:, None, None]+table[i][None, :, None]+table[j][None, None, :]

    low_sum, high_sum, colour_sum = (m.sum(axis=0) for m in (low_w, high_w, high_colour))
    largest_colour = None
    viable_boxes = 0
    for i, j in combinations(range(len(edges)), 2):
        # 目标 W/c² 严格介于 0 与一个取整单位之间。
        possible = (box_sums(low_w, low_sum, i, j) <= 0) & (box_sums(high_w, high_sum, i, j) >= 1)
        if possible.any():
            values = bound0 + box_sums(high_colour, colour_sum, i, j)[possible]
            local = int(values.max())
            largest_colour = local if largest_colour is None else max(local, largest_colour)
            viable_boxes += int(possible.sum())
    assert largest_colour is not None and largest_colour < -722*UNIT
    assert model["beta"][2] == -7
    increase = -model["beta"][2]*(arb("1.22089e19")/10**16).log()/(2*pi)
    assert 0 < increase < 8
    print("小标尺全域 K 上界及升至普朗克标尺的色支增量", k_max.str(16), increase.str(14))
    print("选择数、每组闭盒、通过真空必要条件的盒数、色支上界（精确分子/分母）",
          len(edges)*(len(edges)-1)//2, PARTS**3, viable_boxes, largest_colour, UNIT)


if __name__ == "__main__":
    verify()

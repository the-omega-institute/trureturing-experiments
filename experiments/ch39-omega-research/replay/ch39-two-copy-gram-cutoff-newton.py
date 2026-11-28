"""两份十九线 RH Gram 壳谱的乘积截止在连续公共标尺上的牛顿读数。"""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

from flint import arb


BITS = 8
MESH = 2**BITS
SCALE_BITS = 56
UNIT = 2**SCALE_BITS
LOW = -32
TABLE_LOW = -55


def quantize(value):
    """给出有理近似的分子，同时保留原区间以便检查余量。"""
    return int((value*UNIT).mid().floor().unique_fmpz())


def kernel_table(spectrum, heat):
    table = []
    pi = arb.pi()
    for octave in range(TABLE_LOW, 1):
        for m in range(MESH, 2*MESH if octave < 0 else MESH+1):
            x = arb(m)*arb(2)**(octave-BITS)
            value = x*sum(a*heat(2, z/x) for _, a, z in spectrum)/(2*pi)
            rounded = quantize(value)
            assert abs(value-arb(rounded)/UNIT) < arb(2)/UNIT
            table.append(rounded)
    return table


def stencil(entries, mantissa):
    """聚合同一标尺尾数下所有谱线对的线性插值系数，全部使用整数。"""
    coeff = {}
    for node, weight in entries:
        numerator = node*mantissa
        exponent = numerator.bit_length()-1
        shift = exponent-BITS
        left = numerator >> shift
        remainder = numerator-(left << shift)
        assert MESH <= left < 2*MESH and 0 <= shift <= SCALE_BITS
        fraction = remainder << (SCALE_BITS-shift)
        index = (exponent-BITS-SCALE_BITS)*MESH+left-MESH
        coeff[index] = coeff.get(index, 0)+weight*(UNIT-fraction)
        if fraction:
            coeff[index+1] = coeff.get(index+1, 0)+weight*fraction
    return list(coeff.items())


def bounds(first, second, coherence, table, slope_bound, curvature_bound):
    entries = []
    drift = arb(0)
    absolute_first = arb(0)
    variation = arb(0)
    for (node_a, weight_a), (node_b, weight_b) in product(first, second):
        node = node_a*node_b
        weight = weight_a.real*weight_b.real+coherence*weight_a.imag*weight_b.imag
        n, w = quantize(node), quantize(weight)
        # 数学定义中 λ_max/λ_max 恰为 1，区间运算保留了除法舍入误差。
        assert 0 < n <= UNIT
        w0 = arb(w)/UNIT
        drift += abs(weight-w0)+abs(w0)*abs(node-arb(n)/UNIT)
        absolute_first += abs(weight)*node
        variation += abs(w0)
        entries.append((n, w))
    assert drift < arb("1e-6")
    # 单个热核插值及整组读数对标尺的插值各用一次二阶导数界。
    rounded_first = sum(abs(w)*n for n, w in entries)/arb(UNIT**2)
    endpoint_error = slope_bound*drift+curvature_bound*rounded_first/(8*MESH**2)
    cell_error = curvature_bound*absolute_first/(8*MESH**2)
    table_error = 2*variation/UNIT
    upper = (slope_bound*absolute_first*arb(2)**LOW).upper()
    witness = None
    for m in range(MESH, 2*MESH):
        coefficients = stencil(entries, m)
        first_index = min(i for i, _ in coefficients)
        last = 1 if m == MESH else 0
        for octave in range(LOW, last):
            base = (octave-TABLE_LOW)*MESH
            assert base+first_index >= 0
            x = arb(m)*arb(2)**(octave-BITS)
            approx = arb(sum(v*table[base+i] for i, v in coefficients))/UNIT**3
            error = x*endpoint_error+table_error
            upper = upper.max((approx+error+x*cell_error).upper())
            if m == 35*MESH//32 and octave == -20:
                witness = (approx-error).lower()
    assert witness is not None
    return upper, witness


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    source = {}
    for charge, levels, direct, _ in native["shell_weights"](root):
        nodes = [v/levels[-1] for v in levels[:-1]]+[arb(1)]
        for pivot, row in enumerate(direct, 1):
            source[charge, pivot] = list(zip(nodes, row, strict=True))
    fields = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    _, lib, (spectrum,) = fields(root, (2_345_666,), tight_neutral=True)
    table = kernel_table(spectrum, lib["expint"])
    slope_bound = sum(abs(a) for _, a, _ in spectrum)/(2*arb.pi())
    curvature_bound = slope_bound/arb(1).exp()
    assert slope_bound < arb("31.65")
    winner = ((2, 3), (2, 3))
    results = {}
    budget_pass = 0
    for pair in combinations_with_replacement(source, 2):
        first, second = (source[key] for key in pair)
        first_tv, second_tv = (sum(abs(w.real) for _, w in group) for group in (first, second))
        cost = (first_tv*second_tv-1)/2
        assert cost > arb("18292.7063") or cost < arb("18292.7062")
        budget_pass += int(cost > arb("18292.7063"))
        for coherence in (-1, 1):
            results[pair, coherence] = bounds(first, second, coherence, table, slope_bound, curvature_bound)
    upper = max(results[winner, k][0] for k in (-1, 1))
    witness = min(results[winner, k][1] for k in (-1, 1))
    runner = max(v[0] for (pair, _), v in results.items() if pair != winner)
    print("独立实读出通过有序原子负量粗筛的配方、最佳配方", budget_pass, winner)
    print("同一标尺全混合下界、全标尺上界、其余 119 种上界",
          *[v.str(13) for v in (witness, upper, runner)])
    assert len(results) == 240 and budget_pass == 18
    assert runner < arb(".0554") < arb(".0565") < witness and upper < arb(".057")


if __name__ == "__main__":
    verify()

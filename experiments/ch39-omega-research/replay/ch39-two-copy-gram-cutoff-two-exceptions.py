"""固定交汇能标，逐对核验双副本壳谱的两个独立插值例外。"""
from itertools import combinations, product
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
CUTS = ((".004", ".9", 61_732), (".00425", ".7", 1_101))


def model_inputs():
    native = run_path(str(ROOT / "ch39-native-rh-ray-pvm.py"))
    probes = {}
    for q, levels, direct, _ in native["shell_weights"](ROOT):
        if q in (-4, 2):
            pivot = 0 if q == -4 else 2
            nodes = [value / levels[-1] for value in levels[:-1]] + [arb(1)]
            probes[q] = list(zip(nodes, (w.real for w in direct[pivot]), strict=True))
            assert all(0 < node <= 1 for node in nodes)
    edges = [(s * t, (s + t) / 2, u * v)
             for (s, u), (t, v) in product(probes[-4], probes[2])]
    assert len(edges) == 361
    assert sum(weight for _, _, weight in edges).contains(1)

    rebuild = run_path(str(ROOT / "ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    _, lib, (central,) = rebuild(ROOT, (2_345_666,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 192
    inverse_terms = run_path(str(ROOT / "ch39-upward-crossing-cutoff-inverse-moment.py"))
    constant, envelope, _, _ = inverse_terms["inverse_moment_terms"](ROOT, model, initial, central)
    assert envelope > 0
    log_scale = (arb(10**16) / model["mz"]).log()
    j_floor = (log_scale - constant) / envelope
    assert j_floor > 12_720_563
    return edges, central, lib["expint"], j_floor


def quadratic_minimum(left, middle, right):
    """从二次 Bernstein 系数求极小值，比较不确定时取系数下界。"""
    if middle >= left or middle >= right:
        return left.min(right)
    if middle < left and middle < right:
        left_gap, right_gap = left - middle, right - middle
        return middle + left_gap * right_gap / (left_gap + right_gap)
    return left.min(middle).min(right)


def section_bounds(edges, rows, heat, cut, end):
    """预计算三个区域的公共和与逐线例外修正。"""
    pi = arb.pi()
    k0 = sum(alpha for _, alpha, _ in rows) / (2 * pi)
    w0 = -sum(g for g, _, _ in rows) / (8 * pi)
    assert k0 < 0 < w0
    spin_pos = sum(g * z for g, _, z in rows if g > 0) / (4 * pi)
    spin_neg = sum(-g * z for g, _, z in rows if g < 0) / (4 * pi)

    linear = []
    quadratic = []
    vacuum_error = arb(0)
    vacuum_correction = []
    inverse = []
    k_error = arb(0)
    for a, b, weight in edges:
        v = (1 - end) * a + end * b
        positive, negative = weight.max(0), (-weight).max(0)
        x = cut * b
        k_error += sum(
            (-(weight * alpha)).max(0) * x * (1 - heat(2, z / x))
            for _, alpha, z in rows
        ) / (2 * pi)
        linear.append((weight * a, weight * b, (weight * a).min(weight * b)))
        # 公共二次式在 [0,T] 上的三个 Bernstein 系数，及例外下界。
        aa, ab, vv, bb = weight * a**2, weight * a * v, weight * v**2, weight * b**2
        quadratic.append((aa, ab, vv, aa.min(bb)))
        error = (positive * spin_neg + negative * spin_pos) / cut
        vacuum_error += error * v
        vacuum_correction.append(error * (1 - end) * (b - a))
        inverse.append((abs(weight) / v, abs(weight) / a - abs(weight) / v))

    linear0 = sum(row[0] for row in linear)
    linear1 = sum(row[1] for row in linear)
    q0, q1, q2 = (sum(row[index] for row in quadratic) for index in range(3))
    j_total = sum(value for value, _ in inverse)

    def bounds(i, j):
        left0, left1, left_min = linear[i]
        right0, right1, right_min = linear[j]
        p = (linear0 - left0 - right0).min(linear1 - left1 - right1)
        k_upper = (k0 * cut * (p + left_min + right_min)).max(0) + k_error

        left0, left1, left2, left_min = quadratic[i]
        right0, right1, right2, right_min = quadratic[j]
        q = quadratic_minimum(q0 - left0 - right0, q1 - left1 - right1, q2 - left2 - right2)
        w_lower = w0 * (q + left_min + right_min)
        w_lower -= vacuum_error + vacuum_correction[i] + vacuum_correction[j]
        j_upper = (j_total + inverse[i][1] + inverse[j][1]) / cut
        return k_upper, w_lower, j_upper

    return bounds


def verify():
    edges, rows, heat, j_floor = model_inputs()
    eps = arb("1e-122")
    pending = set(combinations(range(len(edges)), 2))
    example = (1, 19)  # 将 (i,j) 展平为 19*i+j。
    print("固定能标的 J 下界", j_floor.str(16), flush=True)
    for cut_text, end_text, expected in CUTS:
        cut, end = arb(cut_text), arb(end_text)
        assert 0 < cut < 1 and 0 < end < 1
        bounds = section_bounds(edges, rows, heat, cut, end)
        target = eps / cut**2
        excluded = set()
        for i, j in pending:
            k, w, inverse = bounds(i, j)
            if k < 1 and w > target and inverse < j_floor:
                excluded.add((i, j))
        assert len(excluded) == expected
        print(cut_text, end_text, "此轮新增排除", len(excluded), flush=True)
        if (cut_text, end_text) == CUTS[0][:2]:
            assert example in excluded
            k_upper, w_lower, j_upper = bounds(*example)
            assert k_upper < arb(".913") and w_lower > 21 and j_upper < 8_099_000
            print("例外 (0,1) 与 (1,0) 的 K 上界、W/c² 下界、J 上界",
                  *(value.str(12) for value in (k_upper, w_lower, j_upper)), flush=True)
        pending.difference_update(excluded)

    total = comb(len(edges), 2)
    print("全部选择、已排除、尚未排除", total, total - len(pending), len(pending))
    unresolved_lines = [i for i in range(len(edges))
                        if all(tuple(sorted((i, j))) in pending
                               for j in range(len(edges)) if i != j)]
    extra = sorted((divmod(i, 19), divmod(j, 19)) for i, j in pending
                   if i not in unresolved_lines and j not in unresolved_lines)
    print("尚未排除的选择包含任一下列线对",
          [divmod(i, 19) for i in unresolved_lines])
    print("其他尚未排除的例外组合", extra)


if __name__ == "__main__":
    verify()

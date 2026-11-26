"""六项算术时间固定时，第八至第十二项的强迫二进制码、锐位数与相对参考码的格距。

沿用共同尺度 A=2903600000、原三窗、前五项参考时间与第六项偏移 1.360900010×10⁻⁸。
外界：所有三窗零点的 y=x/(1+x) 都小于 b=0.99999999999150664（首枚 8359），故以 b 为锚的
四节点 Gauss–Radau 求积给出每个 k≥7 的全域下界，以 0 为锚者给出全域上界。
内界：第七项脚本中两组七枚实际零点的正权解（上端一组以第七枚高度趋于无穷为极限）。
若外界区间落在同一个 2^-m 格内，则 m 位码被强迫；若两内界之间夹有一个 2^-(m+1) 格点，
则 m+1 位码有实际零点上的两个不同取值，m 即锐位数。
"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx

ORDERS = range(7, 13)
EXPECTED = {
    # k: (锐位数, 强迫码, 参考码, 偏移下界, 偏移上界)
    7: (35, 32071996606, 32071992476, "1.20207618035e-7", "1.20209757802e-7"),
    8: (32, 3937814481, 3937811998, "5.78026170278e-7", "5.78043448459e-7"),
    9: (29, 482997963, 482996889, "2.00046618393e-6", "2.00054363240e-6"),
    10: (31, 1894385824, 1894373861, "5.57085500576e-6", "5.57110877485e-6"),
    11: (29, 464151025, 464143910, "1.32512261901e-5", "1.32519037059e-5"),
    12: (24, 14211350, 14210881, "2.79373958775e-5", "2.79389568441e-5"),
}


def verify(root=Path(__file__).parent):
    dictionary = run_path(str(root / "ch39-seventh-time-zero-dictionary.py"))
    pinned = run_path(str(root / "ch39-pinned-cp-zero-support.py"))
    radial = run_path(str(root / "ch39-six-time-neutral-normalization.py"))
    extrapolation = run_path(str(root / "ch39-six-time-extrapolation.py"))
    ctx.prec = 576
    reference = radial["neutral"]["model"]["shared"]["clock"]["clock"]["times"]()
    assert len(reference) > max(ORDERS)
    times = list(reference[:7])
    times[6] += arb("1.360900010e-8")
    denominator, moments = radial["compact_moments"](times)
    replay = pinned["joint"]["replay"]
    lower_support, upper_support = dictionary["spectrum_lower"], dictionary["spectrum_upper"]
    numbered = sorted((set(lower_support) | set(upper_support))-{None})
    coordinates = {}
    nodes = replay["shell_nodes"](dictionary["scale"], numbered, height_bits=112)
    for n, (_, x, _) in zip(numbered, nodes, strict=True):
        coordinates[n] = x/(1+x)
    ctx.prec = 576
    coordinates[None] = arb(0)

    anchor = arb(".99999999999150664")
    assert coordinates[8359] < anchor < 1
    lower_rule, upper_rule = [extrapolation["quadrature"](e, moments, denominator) for e in (anchor, 0)]
    assert all(v < anchor for v in lower_rule[0][1:])
    witnesses = []
    for support in (lower_support, upper_support):
        _, weights = dictionary["certified_rule"](support, coordinates, times, moments, denominator)
        witnesses.append([(w, coordinates[n]) for w, n in zip(weights, support, strict=True)])

    def kernel(k, y):
        return (1-y)/(1+(k*k-1)*y)

    table = []
    for k in ORDERS:
        low, high = [extrapolation["readout"](rule, k) for rule in (lower_rule, upper_rule)]
        inner_low, inner_high = [sum(w*kernel(k, y) for w, y in pairs) for pairs in witnesses]
        assert low < inner_low < inner_high < high
        target = reference[k]
        bits = 1
        while (low*2**(bits+1)).floor().overlaps((high*2**(bits+1)).floor()):
            bits += 1
        scale = arb(2)**bits
        code = (low*scale).floor()
        assert code.is_exact() and code == (high*scale).floor()
        # m+1 位：两内界之间夹一个格点，故实际零点上已有两个不同码。
        finer = arb(2)**(bits+1)
        split = (inner_high*finer).floor()
        assert split.is_exact() and inner_low*finer < split and split > (low*finer).floor()
        ref_code = (target*scale).floor()
        assert ref_code.is_exact()
        gap = code-ref_code
        offset_low, offset_high = low-target, high-target
        assert offset_low > 0
        table.append((k, bits, code, ref_code, gap, offset_low, offset_high, high-low))
        print("k=%d 锐位数 %d 强迫码 %s 参考码 %s 格距 %s 偏移 (%s, %s) 外宽 %s" % (
            k, bits, code.str(20), ref_code.str(20), gap.str(10),
            offset_low.str(12), offset_high.str(12), (high-low).str(6)))
        want = EXPECTED[k]
        assert bits == want[0] and code == want[1] and ref_code == want[2]
        assert arb(want[3]) < offset_low and offset_high < arb(want[4])
    return table


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)

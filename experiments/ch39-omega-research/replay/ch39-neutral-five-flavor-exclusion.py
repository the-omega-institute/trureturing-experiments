"""五项时间与中性壳响应对平均夸克跃迁和 CP 不变量的约束。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_poly, ctx

root = Path(__file__).parent
cp = run_path(str(root / "ch39-cp-ensemble-test.py"))
neutral = run_path(str(root / "ch39-neutral-clock-response.py"))


def verify():
    ctx.prec, ctx.cap = 512, 3
    times = neutral["model"]["shared"]["clock"]["clock"]["times"]()
    cut = arb(117387924448)
    assert 4000000*cut.sqrt() > arb("1370476848096.31")
    print("精确保留两个夸克概率所需尺度", (4000000*cut.sqrt()).str(19))
    spectrum = neutral["neutral_polynomial"]()
    assert all(v > 0 for v in spectrum(arb_poly([-2,1])))
    assert neutral["shell"]["soft"]["point"]["study"]["offset"] > 1
    max_mean = 1+18*(1-times[1])
    first_efficiency = 19*cut/(36+cut)
    first_flavor = arb("4.65271e-11")*first_efficiency/max_mean
    phase, c0 = cp["flavor"]["singlet"](11661)
    _, first = cp["flavor"]["singlet"](8359)
    variance = first_flavor*(first-c0)**2
    data = cp["two"]["spectral_data"]()
    u, v = cp["two"]["probabilities"](c0, data)
    mass = -u[1]/v[1]
    assert arb("1.16626614") < mass < arb("1.16626615")
    # 全相位支撑线的严格曲率来自两个跃迁概率及 CP 的既有证书。
    pair_gap = arb("3.5297e-8")*variance
    matrix = cp["flavor"]["reference_matrix"]()
    _, _, cp0 = cp["matrix_readings"](matrix,phase)
    cp_loss = arb(".0006")*variance/cp0
    assert pair_gap > arb("1.17256e-21") and cp_loss > arb("6.39464e-13")
    print("接近最小尺度时的味首点份额和均方相位偏差下界",
          first_flavor.str(13),variance.str(13))
    print("两概率线性组合的绝对超额下界",pair_gap.str(13))
    print("平均 |Vub|² 精确时，平均 J 的相对亏额下界",cp_loss.str(13))


if __name__ == "__main__":
    verify()

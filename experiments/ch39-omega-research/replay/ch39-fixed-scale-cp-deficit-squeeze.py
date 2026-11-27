"""固定共同尺度 2903600000 时留出 CP 亏额下确界的双侧夹逼：首枚项的锐化下界与新的八枚实际零点见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
window = run_path(str(root / "ch39-fixed-scale-cp-deficit-window.py"))
floor_rule, quark = window["floor_rule"], window["quark"]
replay, cp, flavor, model = quark["replay"], quark["cp"], quark["flavor"], quark["model"]
neutral = floor_rule["neutral"]

SCALE = 2_903_600_000
X0 = arb(117387924448)
# 簇的次序：首枚、x≈4.93e-3 三枚、x≈1.213e-3 两枚、x≈8.24e-9 两枚。
INDICES = (8359, 142122172741, 142213179754, 142294792965,
           295843057250, 296037909920, 143851497595768, 143935311868267)
BANDS = (("8462", "8463"), ("4.1328e10", "4.1329e10"), ("4.1353e10", "4.1354e10"),
         ("4.1376e10", "4.1377e10"), ("8.3325e10", "8.3326e10"), ("8.3378e10", "8.3379e10"),
         ("3.1985e13", "3.1986e13"), ("3.2003e13", "3.2004e13"))


def sharpened_floor():
    """D ≥ φ(A)·η(x₁)/E⁺·G(c₈₃₅₉)/J₀，其中 E⁺ 是 E_p η 的上界。"""
    ctx.prec = 1024
    scale = arb(SCALE)
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    cap = ((scale/second)**2).upper()
    _, _, _, contact, time_den, _, _, _ = floor_rule["envelope"](cap)
    times = model["shared"]["clock"]["clock"]["times"]()
    moments = neutral["contrast"]["moments_five"](times)
    # 接触多项式 B 在 0 与两个内节点处为零且在 x≥0 非负：E[B/Q] 由前五项时间精确给出。
    expectation = sum(k*v for k, v in zip(contact.coeffs(), moments, strict=True))
    x1 = (scale/first)**2
    floor = window["first_floor"](scale)
    ctx.prec = 1024
    at_first = contact(x1)/time_den(x1)
    ceiling = expectation/at_first
    rest = expectation-floor*at_first
    assert rest > 0 and floor < ceiling
    # 其余支撑 x≤cap；x≥X 部分的算术质量受 ∫_rest B/Q ≤ rest 控制。
    threshold = arb("0.05")
    low, smallest, cells = threshold, None, 0
    while low < cap:
        high = (low*arb("1.02")).min(cap)
        cell = low.union(high)
        value = (contact(cell)/time_den(cell)).lower()
        smallest = value if smallest is None else smallest.min(value)
        low, cells = high, cells+1
    assert smallest > 0
    heavy = rest/smallest
    poly = neutral["neutral_polynomial"]()
    poly = poly*(1/poly[0])
    slope = poly.derivative()
    # η-1 = xP'/P ≤ β·x/(1+x) 于 [0,X]；η-1 ≤ 18 处处成立。
    beta = None
    for j in range(400):
        cell = (threshold*j/400).union(threshold*(j+1)/400)
        value = ((1+cell)*slope(cell)/poly(cell)).upper()
        beta = value if beta is None else beta.max(value)
    # 真实响应 η_γ 用调谐辐边权 1-δ_γ 的拉普拉斯量 L_γ。L_γ = L₀ - δ_γ Σ_tuned L_e，
    # Σ_tuned L_e ⪯ L₀，故 (1-δ_γ)L₀ ⪯ L_γ ⪯ L₀，零模同为常向量；于是
    # η₀ ≤ η_γ 且 η_γ-1 ≤ (η₀-1)/(1-δ_γ)。δ_γ 随高度递减，其余支撑高于 γ₈₃₅₉。
    offset = neutral["shell"]["soft"]["point"]["study"]["offset"]
    assert offset > 0
    delta = arb(23)/25*arb(2)**-20/(first**2+offset)
    assert 0 < delta < arb("1.23e-14")
    effective = beta/(1-delta)
    assert effective < arb("1.35390")
    # 首枚取 η₀(x₁) ≤ η_γ(x₁)：ν₁ 对 η₁ 与 p₁ 都单调增，分母取 p₁ 的上界仍为下界。
    eta1 = 1+x1*slope(x1)/poly(x1)
    upper = 1+ceiling*(eta1-1)+effective*(1-times[1])+18*heavy
    factor, _ = window["forced_factor"]()
    ctx.prec = 1024
    gap = factor*(1+18*(1-times[1]))/(19*X0/(36+X0))
    bound = floor*eta1/upper*gap
    return dict(floor=floor, ceiling=ceiling, rest=rest, smallest=smallest, cells=cells,
                heavy=heavy, beta=beta, delta=delta, effective=effective, upper=upper, eta1=eta1, gap=gap, bound=bound,
                old=floor*factor)


def records(numbers):
    """逐枚重算三窗、真实中性响应、味读数与六面余弦；高度由 acb.zeta_zero 以 448 位给出。"""
    ctx.prec, ctx.cap = 768, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    matrix = flavor["reference_matrix"]()
    phase0, c0 = flavor["singlet"](11661)
    prob0, _, cp0 = cp["matrix_readings"](matrix, phase0)
    rows, cosines, heights = [], [], []
    for height, x, eta in replay["shell_nodes"](SCALE, numbers, height_bits=448):
        ctx.prec = 1024
        phase, cosine = flavor["singlet_at"](flavor["geom"]["phases"](height))
        probabilities, _, cp_value = cp["matrix_readings"](matrix, phase)
        times_at = [model["time_kernel"](x, j) for j in range(1, 6)]
        column = [1, *times_at, eta*model["gradient"](x, poles, slopes),
                  eta*(probabilities[6]-prob0[6])]
        rows.append((column, x, eta, probabilities, cp_value))
        cosines.append(cosine)
        heights.append(height)
    assert all(arb(lo) < h < arb(hi) for h, (lo, hi) in zip(heights, BANDS, strict=True))
    return rows, cosines, c0, times, invariants, ref, prob0, cp0



def verify():
    # 下界：对一切三窗支撑与可行正权，D ≥ φ·η(x₁)/E⁺·G(c₈₃₅₉)/J₀。
    low = sharpened_floor()
    assert arb("4.65277e-11") < low["floor"] < arb("4.65278e-11")
    assert arb("4.66880e-11") < low["ceiling"] < arb("4.66881e-11")
    assert arb("1.1129e-17") < low["rest"] < arb("1.1131e-17")
    assert low["smallest"] > arb("3.1624e-8") and low["heavy"] < arb("3.52e-10")
    assert arb("1.35389") < low["beta"] < arb("1.35390")
    assert low["cells"] == 1405
    assert low["delta"] < arb("1.23e-14") and low["effective"] < arb("1.35390")
    assert arb("1.0021698") < low["upper"] < arb("1.0021699")
    assert arb("0.0444239") < low["gap"] < arb("0.0444240")
    bound = low["bound"]
    assert arb("3.91869e-11") < bound < arb("3.91870e-11")
    assert arb("3.81708e-11") < low["old"] < arb("3.81709e-11")
    assert arb("1.0266") < bound/low["old"] < arb("1.0267")
    ctx.prec = 1024
    times = model["shared"]["clock"]["clock"]["times"]()
    assert arb("1.0288475") < 1+18*(1-times[1]) < arb("1.0288476")
    print("E_p η 的上界", low["upper"].str(10), "；旧上界 1+18(1-T1)")
    print("A=2903600000 上 CP 亏额的强制下界", bound.str(12))

    # 上界：同一尺度上的八枚实际三窗零点见证。
    rows, cosines, c0, times, invariants, ref, prob0, cp0 = records(INDICES)
    ctx.prec = 1024
    weights = window["solve"](rows, times)
    assert all(w > 0 for w in weights)
    assert low["floor"] < weights[0] < low["ceiling"]
    xs = [row[1] for row in rows]
    assert all(arb("4.92e-3") < x < arb("4.94e-3") for x in xs[1:4])
    assert all(arb("1.212e-3") < x < arb("1.215e-3") for x in xs[4:6])
    assert all(arb("8.23e-9") < x < arb("8.25e-9") for x in xs[6:])
    ub, cb, loss = quark["compare"](rows, weights, times, invariants, ref, prob0, cp0)
    assert ub.contains(0)
    sixth = sum(w*model["time_kernel"](row[1], 6) for w, row in zip(weights, rows, strict=True))-times[6]
    assert arb("1.360905797e-8") < sixth < arb("1.360905798e-8")
    assert arb("9.3711639e-10") < loss < arb("9.3711640e-10") and loss.rad() < arb("1e-19")
    assert arb("3.35355e-20") < cb < arb("3.35356e-20")
    distances = [abs(c-c0) for c in cosines[1:]]
    assert all(d < arb("2.28e-6") for d in distances) and abs(cosines[0]-c0) > arb(".0062")
    raw = [w*row[2] for w, row in zip(weights, rows, strict=True)]
    nu = [v/sum(raw) for v in raw]
    _, tangent = window["forced_factor"]()
    ctx.prec = 1024
    first_share = nu[0]*(1-rows[0][4]/cp0+tangent*(rows[0][3][6]-prob0[6]))
    second = sum(v*(c-c0)**2 for v, c in zip(nu, cosines, strict=True))
    assert arb("3.9322e-11") < first_share < arb("3.9323e-11")
    assert arb("1.0034") < first_share/bound < arb("1.0035")
    assert arb("1.1061e-12") < second < arb("1.1062e-12")
    # 与比值律 1.91919011e10 ΔC < D < 6.38699834e10 ΔC 及 581.861 ≤ D/E(c-c0)² ≤ 1936.427 一致。
    assert arb("1.91919011e10")*cb < loss < arb("6.38699834e10")*cb
    assert arb("581.861") < loss/second < arb("1936.427")
    assert arb("847.21") < loss/second < arb("847.22")
    assert arb("2.7944e10") < loss/cb < arb("2.7945e10")
    # 两侧之比与相对已有见证的改进。
    assert arb("23.91") < loss/bound < arb("23.92")
    assert arb("6.9126e-7")/loss > arb("737")
    assert arb("2.5021e-17")/cb > arb("745")
    assert sixth > arb("1.3301556e-8")
    assert arb("8.97e-10") < loss-first_share < arb("8.98e-10")
    print("八枚见证：CP 亏额", loss.str(10), "|Vcb|² 偏移", cb.str(8),
          "其余七枚最大 |c-c0|", max(float(d.upper()) for d in distances))
    print("首枚贡献", first_share.str(8), "味相位二阶矩", second.str(8), "两侧之比", (loss/bound).str(6))


if __name__ == "__main__":
    verify()

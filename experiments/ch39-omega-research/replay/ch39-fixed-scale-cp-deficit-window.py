"""固定共同尺度时留出 CP 亏额的强制下界、尺度包络与八枚实际零点见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
quark = run_path(str(root / "ch39-five-neutral-quark-probability.py"))
floor_rule = run_path(str(root / "ch39-neutral-five-scale-floor.py"))
band = run_path(str(root / "ch39-shell-cp-charm-band.py"))
replay, cp, flavor, model = quark["replay"], quark["cp"], quark["flavor"], quark["model"]
two = cp["two"]

SCALE = 2_903_600_000
# 每簇选取六面余弦最接近参考值、且逐枚通过三窗的实际零点。
indices = (8359, 142276819129, 142268853888, 295884733682, 295676934149,
           13219443393270, 143921301014041, 143921296609177)
X_FLOOR = arb(117387924448)


def forced_factor():
    """首枚味概率与 CP 切线余项的换算因子 K：D ≥ K · p(8359)。"""
    ctx.prec, ctx.cap = 384, 3
    data = two["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())
    reference = flavor["singlet"](11661)[1]
    up, _ = two["probabilities"](reference, data)
    current = cp["cp_value"](reference, polys, jet=True)
    slope = current[1]/up[1]
    first = flavor["singlet"](8359)[1]
    u, _ = two["probabilities"](first, data)
    gap = current[0]+slope*(u[0]-up[0])-cp["cp_value"](first, polys)
    times = model["shared"]["clock"]["clock"]["times"]()
    efficiency = 19*X_FLOOR/(36+X_FLOOR)
    factor = efficiency/(1+18*(1-times[1]))*gap/current[0]
    # 与既有全壳证书一致：首枚份额 4.65271e-11 对应 3.81702e-11。
    assert factor*arb("4.65271e-11") > arb("3.81702e-11")
    # 平均 U 精确时 G/J0 = H + (s/J0)(U-U0)，H = 1-J/J0。
    return factor, slope/current[0]


def first_floor(top):
    """尺度不超过 top 时，五项时间与参考味驻点强制首枚承担的算术份额下界。"""
    ctx.prec = 1024
    first, second = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    cap = ((top/second)**2).upper()
    deficit, high, den, _, _, _, _, _ = floor_rule["envelope"](cap)
    first_max = ((top/first)**2).upper()
    strip = (first_max-X_FLOOR)/128
    reward = None
    for j in range(128):
        cell = (X_FLOOR+j*strip).union(X_FLOOR+(j+1)*strip)
        value = ((cell-cap)*high(cell)/den(cell)).upper()
        reward = value if reward is None else reward.max(value)
    assert reward > 0 and deficit < 0
    return -deficit/reward


def records_at(scale, numbers):
    """沿用八点脚本的逐枚列，同时保留每枚的六面余弦。"""
    ctx.prec, ctx.cap = 768, 3
    times = model["shared"]["clock"]["clock"]["times"]()
    invariants, ref, poles, slopes = model["shared"]["flavor_data"]()
    matrix = flavor["reference_matrix"]()
    phase0, c0 = flavor["singlet"](11661)
    prob0, _, cp0 = cp["matrix_readings"](matrix, phase0)
    rows, cosines, heights = [], [], []
    # 两枚高处零点的时间列只在 x 的一阶上可分，故高度与列均取高精度。
    for height, x, eta in replay["shell_nodes"](scale, numbers, height_bits=448):
        ctx.prec = 1024
        phase, cosine = flavor["singlet_at"](flavor["geom"]["phases"](height))
        probabilities, _, cp_value = cp["matrix_readings"](matrix, phase)
        times_at = [model["time_kernel"](x, j) for j in range(1, 6)]
        column = [1, *times_at, eta*model["gradient"](x, poles, slopes),
                  eta*(probabilities[6]-prob0[6])]
        rows.append((column, x, eta, probabilities, cp_value))
        cosines.append(cosine)
        heights.append(height)
    bands = (("8462", "8463"), ("4.137e10", "4.138e10"), ("4.136e10", "4.137e10"),
             ("8.333e10", "8.334e10"), ("8.328e10", "8.329e10"), ("3.2e12", "3.2001e12"),
             ("3.2e13", "3.2001e13"), ("3.2e13", "3.2001e13"))
    assert all(arb(lo) < h < arb(hi) for h, (lo, hi) in zip(heights, bands, strict=True))
    return rows, cosines, c0, times, invariants, ref, prob0, cp0


def solve(rows, times):
    """八阶方程用预条件球算术求解；病态来自两枚 x 几乎相同的高处零点。"""
    matrix = arb_mat(list(zip(*(row[0] for row in rows), strict=True)))
    target = arb_mat([[v] for v in [*times[:6], 0, 0]])
    weights = matrix.solve(target, algorithm="precond")
    assert all(v.contains(0) for v in (matrix*weights-target).entries())
    return weights.entries()


def verify():
    factor, tangent = forced_factor()
    assert arb("43840.494666") < tangent < arb("43840.494667")
    assert arb(".8203") < factor < arb(".8204")

    # 下界：任何尺度 ≤ top、任意三窗支撑与正权，平均 |Vub|² 精确时 D ≥ K·floor(top)。
    at_scale = first_floor(arb(SCALE))*factor
    assert arb("3.81708e-11") < at_scale < arb("3.81709e-11")
    print("A=2903600000 及更低尺度上 CP 亏额的强制下界", at_scale.str(12))
    assert first_floor(arb(2_903_616_618)) > arb("4.65271e-11")
    assert (arb(SCALE)-arb("2899267716.58"))/arb("2899267716.58") < arb(".0015")
    # 指定 ΔC=2e-12 时由比值律强制的亏额。
    assert arb("1.91919011e10")*arb("2e-12") > arb(".0383838")
    envelope = (("3.5e-11", 3_007_568_000), ("3e-11", 3_169_126_000), ("2e-11", 3_481_405_000),
                ("1e-11", 3_777_881_000), ("1e-12", 4_031_418_000), ("1e-13", 4_056_112_000))
    for epsilon, top in envelope:
        # 在 top+1 处同样成立，故 A_eps ≥ top+1 > top。
        for limit in (top, top+1):
            bound = first_floor(arb(limit))*factor
            assert bound > arb(epsilon)
        print("尺度不超过", top+1, "时 D 至少", bound.str(10), ">", epsilon)
    assert arb(4_056_112_000) < arb("4058848884.28")

    # 上界：同一尺度上的八枚实际零点见证。
    rows, cosines, c0, times, invariants, ref, prob0, cp0 = records_at(SCALE, indices)
    ctx.prec = 1024
    weights = solve(rows, times)
    assert all(w > 0 for w in weights)
    assert arb("4.65271e-11") < weights[0] < arb("4.66881e-11")
    xs = [row[1] for row in rows]
    assert all(arb("4.9e-3") < x < arb("4.95e-3") for x in xs[1:3])
    assert all(arb("1.2e-3") < x < arb("1.22e-3") for x in xs[3:5])
    assert arb("8.2e-7") < xs[5] < arb("8.25e-7")
    assert all(arb("8.2e-9") < x < arb("8.25e-9") for x in xs[6:])
    assert arb("5.8e-8") < abs(xs[6]/xs[7]-1) < arb("6.0e-8")
    ub, cb, loss = quark["compare"](rows, weights, times, invariants, ref, prob0, cp0)
    assert ub.contains(0)
    sixth = sum(w*model["time_kernel"](row[1], 6) for w, row in zip(weights, rows, strict=True))-times[6]
    assert arb("1.360900045e-8") < sixth < arb("1.360900046e-8")
    assert arb("6.9126e-7") < loss < arb("6.9127e-7")
    assert arb("2.5021e-17") < cb < arb("2.5022e-17")
    distances = [abs(c-c0) for c in cosines[1:]]
    assert all(d < arb("4.9e-5") for d in distances)
    assert abs(cosines[0]-c0) > arb(".0062")
    raw = [w*row[2] for w, row in zip(weights, rows, strict=True)]
    nu = [v/sum(raw) for v in raw]
    first_share = nu[0]*(1-rows[0][4]/cp0+tangent*(rows[0][3][6]-prob0[6]))
    second_moment = sum(v*(c-c0)**2 for v, c in zip(nu, cosines, strict=True))
    assert arb("3.93e-11") < first_share < arb("3.94e-11")
    print("味相位二阶矩", second_moment.str(10))
    assert arb("8.0e-10") < second_moment < arb("8.5e-10")
    assert arb("581.861") < loss/second_moment < arb("1936.427")
    # 与比值律 1.91919011e10 ΔC < D < 6.38699834e10 ΔC 一致。
    assert arb("1.91919011e10")*cb < loss < arb("6.38699834e10")*cb
    assert arb("2.762e10") < loss/cb < arb("2.763e10")
    print("八枚见证：CP 亏额、|Vcb|² 偏移、其余七枚的最大 |c-c0|",
          loss.str(10), cb.str(8), max(float(d.upper()) for d in distances))
    print("首枚对亏额的贡献与味相位二阶矩", first_share.str(6), second_moment.str(6))
    assert loss/at_scale < arb("18200") and loss/at_scale > arb("18100")
    assert arb("1.02") < first_share/at_scale < arb("1.04")
    assert loss.rad() < arb("2e-17")
    assert arb("1.90254596927e-12")/cb > arb("76000")
    assert arb(".06083301837")/loss > arb("88000")


if __name__ == "__main__":
    verify()

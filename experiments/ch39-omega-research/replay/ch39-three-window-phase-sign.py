"""在两圈阈值和偶荷重模脱耦的原处方下，三窗能否留下负的六面味坐标。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
window = run_path(str(root / "ch39-common-step-zero-window.py"))
physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
density = run_path(str(root / "ch39-curvature-window-density-cap.py"))
charges = (1, 3)
frequencies = range(-3, 4)

# 端点以千分之一为单位，每组三元组的末项为谱重数。方框编号给出靠近半周的钟。
bands = {
    (1, 0): [
        (2176, 2216, 1), (7798, 7829, 1), (8004, 8038, 1), (10151, 10202, 1),
        (10259, 10308, 1), (10987, 11008, 1), (11468, 11507, 1), (14979, 15013, 3),
        (15632, 15668, 2), (17018, 17046, 1), (17990, 18029, 1), (18111, 18156, 1),
        (19070, 19097, 2), (19176, 19198, 1), (22189, 22234, 1),
    ],
    (1, 1): [
        (2177, 2217, 1), (7767, 7796, 1), (8036, 8070, 1), (10173, 10224, 1),
        (10236, 10286, 1), (10983, 11004, 1), (11472, 11510, 1), (14979, 15013, 3),
        (15631, 15669, 2), (17018, 17046, 1), (17982, 18021, 1), (18119, 18164, 1),
        (19070, 19100, 2), (19173, 19195, 1), (22189, 22234, 1),
    ],
    (3, 0): [
        (2145, 2326, 1), (7542, 7619, 1), (8094, 8242, 1), (10061, 10347, 1),
        (10373, 10596, 1), (10967, 10995, 1), (11371, 11536, 1), (14930, 14951, 1),
        (14985, 15025, 2), (15611, 15701, 2), (16961, 17041, 1), (17793, 17925, 1),
        (18020, 18279, 1), (19000, 19029, 1), (19073, 19102, 1), (19293, 19322, 1),
        (22191, 22394, 1),
    ],
    (3, 1): [
        (2155, 2336, 1), (7458, 7532, 1), (8168, 8329, 1), (10119, 10548, 2),
        (10927, 10961, 1), (11401, 11558, 1), (14930, 14951, 1), (14985, 15025, 2),
        (15610, 15705, 2), (16960, 17041, 1), (17772, 17904, 1), (18045, 18302, 1),
        (19005, 19033, 1), (19073, 19102, 1), (19286, 19315, 1), (22191, 22394, 1),
    ],
}


def schur_coefficients(charge, level):
    cut = window["audit"]["cut_matrix"](charge)
    for i in range(19):
        cut[i, i] -= level
    core, inner = window["core"], window["inner"]
    sub = window["sub"]
    inverse = sub(cut, core, core).inv()
    fixed = sub(cut, inner, inner)
    cross = sub(cut, core, inner)
    edges = [(core.index(i), inner.index(j), d, -u**charge)
             for i, j, d, u in window["clock_edges"]]
    roots = [acb(0, 2*arb.pi()*i/7).exp() for i in range(7)]
    samples = acb_mat(7, 7)
    for i, z in enumerate(roots):
        for j, w in enumerate(roots):
            block = acb_mat(cross)
            for a, b, d, value in edges:
                block[a, b] = value * (z if d == 724 else w)
            samples[i, j] = (fixed - block.conjugate().transpose()*inverse*block).det().real
    transform = acb_mat([[roots[-k*i % 7] for i in range(7)] for k in frequencies])
    return transform * samples * transform.transpose() / 49


def determinant_sign_on_cell(coefficients, charge, corner):
    radius = arb.pi()*charge/16
    second = sum((a*a+b*b)*abs(coefficients[i, j])
                 for i, a in enumerate(frequencies) for j, b in enumerate(frequencies))

    @cache
    def powers(n, index, clock):
        center = arb.pi() if corner == clock else arb(0)
        phase = acb(0, center + radius*(2*arb(index)/n-1)).exp()
        return [phase**k for k in frequencies]

    @cache
    def value(n, i, j):
        row = acb_mat([powers(n, i, 0)])
        column = acb_mat([powers(n, j, 1)]).transpose()
        return (row * coefficients * column)[0, 0].real

    def cell(n, i, j):
        values = [value(n, i+a, j+b) for a in (0, 1) for b in (0, 1)]
        error = (2*radius/n)**2 * second/8
        return all(v > error for v in values) or all(v < -error for v in values)

    return cell


def certify_no_crossing(coefficients, charge, corner):
    definite = determinant_sign_on_cell(coefficients, charge, corner)
    pending = [(1, 0, 0)]
    while pending:
        n, i, j = pending.pop()
        if definite(n, i, j):
            continue
        assert n < 16
        pending.extend((2*n, 2*i+a, 2*j+b) for a in (0, 1) for b in (0, 1))


def certified_thresholds(model, corner):
    unit = arb("2.435e18")/724**2
    delta = arb(23)/25*arb(2)**-20/14**2
    clocks = {724: -1 if corner == 0 else 1, 804: -1 if corner == 1 else 1}
    knots = {}
    for charge in charges:
        count = 0
        masses = []
        for low, high, multiplicity in bands[charge, corner]:
            for edge, expected in ((low, count), (high, count+multiplicity)):
                level = arb(edge)/1000
                assert window["supplier"](charge, level)(clocks)+1 == expected
                certify_no_crossing(schur_coefficients(charge, level), charge, corner)
            count += multiplicity
            eig = (arb(low)/1000 - 4*delta).union(arb(high)/1000)
            center = (15-delta/3).union(15)
            masses.extend([eig*abs(eig-center)] * multiplicity)
        assert count == 19
        ordered = physical["joint"]["chiral"]["ordered"](masses)
        knots[charge] = [(mass*unit/model["mz"]).log() for mass in ordered[3:]]
        assert all(t > 25 for t in knots[charge])
    return knots


def dot(left, right):
    return sum(x*y for x, y in zip(left, right, strict=True))


def gauge_test(initial, knots, fields, alpha):
    gauge = physical["joint"]["gauge"]["sm"]
    beta = gauge["BETA"]
    one = {q: row[1] for q, row in zip(charges, fields, strict=True)}
    two = {q: row[2] for q, row in zip(charges, fields, strict=True)}
    start = [*initial[:2], 1/alpha]
    correction = [arb(0) for _ in range(3)]
    logs = [arb(x).log() for x in (".98", ".016", ".010")]
    pi = arb.pi()
    step = arb(1)/256
    undecided = 0
    tested = 0
    margin = arb(100)
    # 覆盖 L 从 0 到 29.4；一圈部分保留阈值正部的精确表达式。
    for index in range((294*256+9)//10):
        lo, hi = step*index, step*(index+1)
        time = lo.union(hi)
        opened = {}
        elapsed = {}
        for q in charges:
            if hi <= 25:
                opened[q] = elapsed[q] = arb(0)
            else:
                opened[q] = arb(sum(t < lo for t in knots[q])).union(sum(not (t > hi) for t in knots[q]))
                elapsed[q] = sum((time-t).max(0) for t in knots[q])
        inverses = []
        for i in range(3):
            linear = beta[i]*time + sum(one[q][i]*elapsed[q] for q in charges)
            inverses.append(start[i] - linear/(2*pi) + correction[i] + arb(0, step*arb(".4")))
        assert all(x > 7 for x in inverses)
        couplings = [1/x for x in inverses]
        yukawa = [(2*(v+arb(0, step*arb(".12")))).exp() for v in logs]
        matrix = [[gauge["B"][i][j] + sum(two[q][i][j]*opened[q] for q in charges)
                   for j in range(3)] for i in range(3)]
        dx = [-dot(r, couplings)/(8*pi**2) + dot(d, yukawa)/(32*pi**3)
              for r, d in zip(matrix, gauge["D"])]
        trace = 3*(yukawa[0]+yukawa[1]) + yukawa[2]
        differences = (yukawa[0]-yukawa[1], yukawa[1]-yukawa[0], yukawa[2])
        dy = [(arb(3)/2*d+trace-4*pi*dot(row, couplings))/(16*pi**2)
              for d, row in zip(differences, gauge["C"])]
        assert all(abs(v) < arb(".4") for v in dx) and all(abs(v) < arb(".12") for v in dy)

        def residual(vector):
            running = dot(beta, vector)*time + sum(dot(one[q], vector)*elapsed[q] for q in charges)
            return (dot(start, vector) + dot(correction, vector)
                    + (step*dot(dx, vector)).union(0) - running/(2*pi))

        if residual((1, -1, 0)).contains(0):
            tested += 1
            tests = [residual((1-k, k, -1)) for k in (-1, -.5, 0, .5, 1, 1.5, 2)]
            separation = max(abs(v).lower() for v in tests)
            margin = margin.min(separation)
            if separation <= 0:
                undecided += 1
        for i in range(3):
            correction[i] += step*dx[i]
            logs[i] += step*dy[i]
    assert tested > 0
    return undecided, margin


def verify():
    model, initial = physical["model_inputs"]()
    ctx.prec = 160
    density["restrict_corners"]()
    gauge = physical["joint"]["gauge"]
    # 用原有场表产生荷一费米子与荷三费米子加标量的标准一圈、两圈系数。
    fields = gauge["two"]["thresholds"](model, {1: [0, 0, 0, arb(1)], 3: [0, 0, 0, arb(2)]}, False)
    assert [(row[3], row[4]) for row in fields] == [(1, 0), (1, 1)]
    # 利用同一中性谱及辐边扰动约束偶重模仍脱耦时的最大能标。
    neutral = model["shell"]["shell_spectra"]((0,), zero_index=11431)[0]
    delta = arb(23)/25*arb(2)**-20/14**2
    assert any(v*abs(v-15) + (8*80+12)*delta < arb("3.59") for v in neutral)
    assert (32*arb("3.59")*arb("2.435e18")/724**2/model["mz"]).log() < arb("29.4")
    for corner in (0, 1):
        knots = certified_thresholds(model, corner)
        missed, margin = gauge_test(initial, knots, fields, arb(".117").union(arb(".119")))
        assert missed == 0 and margin > (arb(".072") if corner == 0 else arb(".066"))
        print("靠近半周的钟、两圈交点剩余格数、逆耦合线性对比余量", corner+1, missed, margin.str(10))
        if corner == 1:
            # 上一开窗实算强耦合约 0.12121，该对照应允许交点。
            control, _ = gauge_test(initial, knots, fields, arb(".1211").union(arb(".1213")))
            assert control > 0
    flavor = physical["joint"]["flavor"]
    a, b = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    assert a*(arb.pi()/8).cos() > abs(b)
    cap = arb(246020)/4**14/4
    assert cap < arb(".000229124")
    print("三窗相对上密度界", cap.str(14))


if __name__ == "__main__":
    verify()

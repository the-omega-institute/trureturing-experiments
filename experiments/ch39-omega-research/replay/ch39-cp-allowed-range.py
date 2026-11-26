"""固定算术时间、参考驻点与 |Vub|² 时，平均 J 可取哪些值。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


root = Path(__file__).parent
joint = run_path(str(root / "ch39-common-vub-weights.py"))
model = joint["model"]


def verify():
    physical = model["physical"]
    parameters, initial = physical["model_inputs"]()
    flavor = physical["joint"]["flavor"]
    chiral = physical["joint"]["chiral"]
    supplier = chiral["source"]["spectral_supplier"](parameters)
    ctx.prec, ctx.cap = 432, 3
    readout = flavor["reference_readout"]()
    ref = readout(flavor["singlet"](11661)[0])
    sides = []
    for n in (58214, 2529893):
        height = acb.zeta_zero(n).imag
        rows = supplier(height)
        layers = {q: chiral["ordered"]([v for v, charge in rows if charge == q])
                  for q in (0, *parameters["charges"])}
        _, alpha, curve = physical["layer_observables"](parameters, initial, layers)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
        sides.append(readout(flavor["singlet"](n)[0]))
    fraction = (ref[2]**2 - sides[1][2]**2)/(sides[0][2]**2-sides[1][2]**2)
    assert 0 < fraction < 1
    endpoint = 1 - (fraction*sides[0][3]+(1-fraction)*sides[1][3])/ref[3]
    assert arb(".12276863") < endpoint < arb(".12276865")
    print("相同 |Vub|² 的相位三角形下边界相对亏额", endpoint.str(12))

    indices = (8359, 1674503460, 1674725340, 6025261094, 6025461713,
               1674353576, 1675700399, 6025296660, 6026035086)
    for n in (1675700399, 6026035086):
        _, alpha, curve = physical["observables"](n)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
    ctx.prec = 432
    times = model["shared"]["clock"]["clock"]["times"]()
    _, _, lambdas, slopes = model["shared"]["flavor_data"]()
    columns, cp, fifth = [], [], []
    for n in indices:
        response = readout(flavor["singlet"](n)[0])
        x = (joint["scale"]/acb.zeta_zero(n).imag)**2
        columns.append([1, *[1/(1+j*j*x) for j in range(1, 5)], model["gradient"](x,lambdas,slopes), response[2]**2])
        cp.append(response[3])
        fifth.append(1/(1+25*x))
    target = arb_mat([[v] for v in [*times[:5], 0, ref[2]**2]])
    deviations, fifth_readings = [], []
    for support in combinations(range(len(indices)), 7):
        matrix = arb_mat(list(zip(*(columns[i] for i in support), strict=True)))
        weights = matrix.solve(target).entries()
        if all(w > 0 for w in weights):
            deviations.append(1-sum(w*cp[i] for w, i in zip(weights, support))/ref[3])
            fifth_readings.append(sum(w*fifth[i] for w, i in zip(weights, support)))
        else:
            assert any(w < 0 for w in weights)
    assert len(deviations) == 4
    ordered = sorted(deviations, key=float)
    assert all(a < b for a, b in zip(ordered, ordered[1:]))
    assert arb(".01409299") < ordered[0] < arb(".01409300")
    assert arb(".02131550") < ordered[-1] < arb(".02131551")
    assert all(arb("2.49764954e-7") < times[5]-value < arb("2.49764956e-7") for value in fifth_readings)
    for bits, target_code, model_codes in (
            (21, 2020141, {2020141}),
            (22, 4040283, {4040282})):
        assert (times[5]*2**bits).floor().unique_fmpz() == target_code
        assert {(value*2**bits).floor().unique_fmpz() for value in fifth_readings} == model_codes
    codes = {(value*2**53).floor().unique_fmpz() for value in fifth_readings}
    assert codes == {8676439645648768, 8676439645648769}
    assert len({code//2 for code in codes}) == 1
    print("第五个时间的 22 位参考码与九点模型码", 4040283, 4040282)
    print("九枚可行域的顶点数与 J 相对亏额区间的两个端点", len(ordered), ordered[0].str(13), ordered[-1].str(13))


if __name__ == "__main__":
    verify()

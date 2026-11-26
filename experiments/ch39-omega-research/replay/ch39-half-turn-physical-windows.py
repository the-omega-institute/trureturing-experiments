"""一只壳钟转过半周时，检查四个物理分支对两组开相位窗的反应。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
cp = run_path(str(root / "ch39-cp-ensemble-test.py"))


def verify():
    params, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    supplier = chiral["source"]["spectral_supplier"](params)
    matrix = cp["flavor"]["reference_matrix"]()
    phase, _ = cp["flavor"]["singlet"](11661)
    reference, _, current = cp["matrix_readings"](matrix,phase)
    ctx.prec, ctx.cap = 384, 3
    spectral = cp["two"]["spectral_data"]()
    cp_polys = cp["cp_data"](spectral,matrix)
    radius = arb(2)**-32
    alpha, curves, up = [], [], []
    window_layers = []
    for second in (2,66):
        clocks = {d: acb(0,arb.pi()*r/64).exp() for d,r in zip((724,804),(1,second),strict=True)}
        rows = supplier(arb(10)**13, clocks=clocks)
        layers = {}
        for q in (0,*params["charges"]):
            sorted_masses = chiral["ordered"]([v for v,t in rows if q == t])
            layers[q] = [v+arb(0,4*80*q*radius) for v in sorted_masses]
        _, coupling, curvature = physical["layer_observables"](params,initial,layers)
        window_layers.append(layers)
        broad = {d:z*acb(0,arb(0,radius)).exp() for d,z in clocks.items()}
        _, coordinate = cp["flavor"]["singlet_at"](broad)
        pu, pc = cp["two"]["probabilities"](coordinate,spectral)
        value = cp["cp_value"](coordinate,cp_polys)
        alpha.append(coupling)
        curves.append(curvature)
        up.append(pu[0]/reference[6])
        print("第二钟角的 π/64 倍数、强耦合、曲率比、U/U₀、C/C₀、J/J₀", second,
              *(v.str(12) for v in (coupling,curvature,up[-1],pc[0]/reference[7],value/current)))
    assert arb(".1184") < alpha[0] < arb(".1185")
    assert arb(".1211") < alpha[1] < arb(".1213")
    assert alpha[1]-alpha[0] > arb(".0027")
    assert all(arb("1.00072") < r < arb("1.00073") for r in curves)
    assert arb(".34719") < up[0] < arb(".34720") and 131 < up[1] < 132

    # 在原窗与半周窗之间逐项替换奇数荷，检验两部分阈值效应的方向。
    assert {q for q in params["charges"] if q % 2} == {1, 3}
    mixed = {(0, 0): alpha[0], (1, 1): alpha[1]}
    for charge_one, charge_three in ((0, 1), (1, 0)):
        layers = dict(window_layers[0])
        layers[1] = window_layers[charge_one][1]
        layers[3] = window_layers[charge_three][3]
        _, coupling, _ = physical["layer_observables"](params, initial, layers)
        mixed[charge_one, charge_three] = coupling
    assert arb(".11677") < mixed[0, 1] < arb(".11679")
    assert arb(".12298") < mixed[1, 0] < arb(".12300")
    for other in (0, 1):
        assert mixed[1, other] - mixed[0, other] > arb(".0044")
        assert mixed[other, 0] - mixed[other, 1] > arb(".00164")
    for choice, coupling in sorted(mixed.items()):
        print("荷 1、荷 3 是否取半周窗，反推强耦合", *choice, coupling.str(12))


if __name__ == "__main__":
    verify()

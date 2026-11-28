"""固定夸克复源及两标尺时，任意高的实际零点对满足十项窗且轨道配重费用低于中央。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def verify():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay/f"ch39-{name}.py"))

    model = load("two-zero-quark-window")
    # 原记录直接在各荷的两钟角盒及全部径向软化区间上验证六质量、三角与 CP。
    assert len(model["values"]) == 10 and all(abs(v-r) < w/8
        for v,r,w in zip(model["values"],model["ref"],model["tol"],strict=True))
    shell, pump = load("three-ray-chiral-boundary"), load("native-light-projector-pump")
    geometry = load("geometric-selector-mixture")
    optimal = load("minimum-orbit-specific-probability")["setup"]()
    data = json.loads((replay/"ch39-nine-orbit-cost-bases.json").read_text())
    basis = next(e["support"] for e in data if e["offsets"] == [0,0])
    center = optimal((0,0),basis)
    geometry["geometry"](shell,pump)
    adj = pump["adjoint"]
    delta = arb(23)/25*arb(2)**-20/(arb(500000)**2+shell["soft"]["point"]["study"]["offset"])
    assert 0 < delta and model["deficit"].contains(0) and model["deficit"].contains(delta)
    energy_radius = arb("2e-10")
    assert 18*model["radius"]+6*delta+arb("0.5e-46") < energy_radius
    cases = (
        (("-.0013849","-.00091671",".00021066","-.00021829"), ".121", "3.05e-7"),
        (("-.00411497","-.00122105",".00131647","-.00115155"), ".901", "9.10e-7"),
    )
    # 在直积编号的计数中，两荷分别计数后相乘，无需不同零点的独立性假设。
    pair_fraction = (model["radius"]/arb.pi())**4
    assert pair_fraction > arb("1.0265e-46")
    assert abs(arb(cases[0][0][0])-arb(cases[1][0][0])) > 2*model["radius"]
    for shifts, window, lower in cases:
        roots, outputs = [], []
        for i,(q,x,y,scale) in enumerate(model["PHASE_BOXES"]):
            a = arb(x)+arb(shifts[2*i])
            b = arb(y)+arb(shifts[2*i+1])
            # prepare 在两角各自半宽 1e-11 及完整径向亏额区间上执行。
            y_matrix = model["prepare"](q,a,b)*model["K"]*arb(scale)
            outputs.append(model["eig"](adj(y_matrix)*y_matrix))
            base, links = pump["clock_matrices"](shell,q)
            off = acb(0,a).exp()*links[0]+acb(0,b).exp()*links[1]
            laplacian = base+off+adj(off)
            for frame in pump["rotation_bases"](shell):
                roots.append(geometry["minimum"](adj(frame)*laplacian*frame))
        (up,u),(down,d) = outputs
        values = [v.sqrt() for v in up+down]+model["mixing_angles"](u,d)
        error = max((abs(v-r)/w).upper() for v,r,w in
                    zip(values,model["ref"],model["tol"],strict=True))
        cost = optimal(None,basis,roots=roots,root_radius=energy_radius)
        improvement = center-cost
        assert error < arb(window) and improvement > arb(lower)
        print("十项最大相对半宽与中央费用减替代费用",error.str(12),improvement.str(10))



if __name__ == "__main__":
    verify()

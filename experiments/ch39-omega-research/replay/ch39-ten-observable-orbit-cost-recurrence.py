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
    eps = arb("1e-22")
    energy_radius = arb("1e-20")
    delta = arb(23)/25*arb(2)**-20/(arb(10**9)**2+shell["soft"]["point"]["study"]["offset"])
    assert 0 < delta and delta < arb("4e-18")
    assert model["deficit"].contains(arb(0).union(delta))
    # 每个参考角周围的小开方形在旧十项盒内。单位钟相位变化不超过角变化。
    assert arb("6e-12")+eps < model["radius"]
    assert 18*eps+6*delta+arb("0.5e-46") < energy_radius
    shifts = {-4:(-1,-1),2:(1,-1)}
    roots = []
    for q,x,y,_ in model["PHASE_BOXES"]:
        base, links = pump["clock_matrices"](shell,q)
        angles = [arb(a)+v*arb("6e-12") for a,v in zip((x,y),shifts[q],strict=True)]
        off = sum((acb(0,a).exp()*link for a,link in zip(angles,links,strict=True)),0*base)
        laplacian = base+off+adj(off)
        for frame in pump["rotation_bases"](shell):
            roots.append(geometry["minimum"](adj(frame)*laplacian*frame))
    cost = optimal(None,basis,roots=roots,root_radius=energy_radius)
    improvement = center-cost
    assert improvement > arb("3.71e-15")
    print("两个开相位盒中，中央最小费用减去替代最小费用",improvement.str(12))


if __name__ == "__main__":
    verify()

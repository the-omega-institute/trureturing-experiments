"""十九点壳同一味概率下，两个留出读数的全相位定量关系。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

previous = run_path(str(Path(__file__).with_name("ch39-neutral-five-flavor-exclusion.py")))
cp = previous["cp"]
two = cp["two"]


def verify():
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    axes = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    radius = arb(".998811")
    assert sum(v*v for v in axes) < radius**2
    data = two["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())
    reference = flavor["singlet"](11661)[1]
    up, charm = two["probabilities"](reference, data)
    current = cp["cp_value"](reference, polys, jet=True)
    slope = current[1]/up[1]
    k = -up[1]/charm[1]
    def gap(location):
        u, _ = two["probabilities"](location, data)
        return current[0]+slope*(u[0]-up[0])-cp["cp_value"](location, polys)

    # 在参考点处切线余项有二重零，用 Taylor 二次系数覆盖邻域。
    lower, upper = arb(".00063"), arb(".061")
    nearby = reference+arb(0, arb(".00099"))
    for i in range(-16, 16):
        region = (reference+arb(i)/16000).union(reference+arb(i+1)/16000)
        u, _ = two["probabilities"](region, data)
        curvature = slope*u[2]-cp["cp_value"](region, polys, jet=True)[2]
        assert lower < curvature < upper
    pending = [(-radius, radius)]
    while pending:
        left, right = pending.pop()
        region = left.union(right)
        if nearby.contains(region):
            continue
        square = (region-reference)**2
        value = gap(region)
        if value > lower*square and value < upper*square:
            continue
        assert right-left > arb("1e-10")
        middle = (left+right)/2
        pending.extend(((left,middle), (middle,right)))

    # 复用两概率证书的全域二阶 Taylor 系数上下界，换算为可观测单位。
    low_ratio = arb("3.5297e-8")*current[0]/(k*upper)
    high_ratio = arb("3.5365e-8")*current[0]/(k*lower)
    assert low_ratio > arb("1.54e-11") and high_ratio < arb("1.51e-9")
    first = flavor["singlet"](8359)[1]
    times = previous["neutral"]["model"]["shared"]["clock"]["clock"]["times"]()
    x_floor = arb(117387924448)
    flavor_floor = arb("4.65271e-11")*19*x_floor/((36+x_floor)*(1+18*(1-times[1])))
    first_loss = flavor_floor*gap(first)/current[0]
    assert first_loss > arb("3.81702e-11")
    print("两个留出读数的倍率界", low_ratio.str(14), high_ratio.str(14))
    print("首枚的切线亏额和强制 CP 相对亏损下界", gap(first).str(14), first_loss.str(14))


if __name__ == "__main__":
    verify()

"""七枚实际三窗零点与一次中性壳响应的五项精确时间。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
neutral = run_path(str(root / "ch39-neutral-clock-response.py"))
physical = run_path(str(root / "ch39-first-three-mode-observables.py"))
indices = (8359, 143921293989490, 3608485756267, 1742020063894,
           13219432687559, 3264638580330, 1736367964191)
scale = 32_000_000_000


def shell_nodes(scale, indices, *, height_bits=112):
    params, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    supply_spectrum = chiral["source"]["spectral_supplier"](params)
    for index in indices:
        ctx.prec = height_bits
        height = acb.zeta_zero(index).imag
        ctx.prec = 224
        if index == 8359:
            _, alpha, curvature = physical["observables"](index)
        else:
            rows = supply_spectrum(height)
            levels = {q: chiral["ordered"]([v for v, charge in rows if charge == q])
                      for q in (0, *params["charges"])}
            _, alpha, curvature = physical["layer_observables"](params, initial, levels)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curvature < arb("1.002")

        ctx.prec = 512
        x = (scale/height)**2
        shifted = neutral["shell"]["laplacian"](0, height)
        for i in range(19):
            shifted[i, i] += x
        response = (x*shifted.inv().trace()).real
        assert 1 < response < 19
        print(index, "强耦合、四层曲率比、中性响应",
              alpha.str(10), curvature.str(12), response.str(12), flush=True)
        yield height, x, response


def verify_at(scale, indices):
    ctx.prec, ctx.cap = 512, 3
    model = neutral["model"]
    times = model["shared"]["clock"]["clock"]["times"]()
    invariants, reference, poles, slopes = model["shared"]["flavor_data"]()
    nodes, efficiencies, columns = [], [], []
    for _, x, response in shell_nodes(scale, indices):
        time_column = [model["time_kernel"](x, j) for j in range(1, 6)]
        columns.append([1, *time_column, response*model["gradient"](x, poles, slopes)])
        nodes.append(x)
        efficiencies.append(response)

    matrix = arb_mat(list(zip(*columns, strict=True)))
    target = arb_mat([[v] for v in [*times[:6], 0]])
    weights = matrix.solve(target)
    assert all(w > 0 for w in weights.entries())
    assert all(v.contains(0) for v in (matrix*weights-target).entries())
    flavor = [w*eta for w, eta in zip(weights.entries(), efficiencies, strict=True)]
    normalizer = sum(flavor)
    neutral["verify_maximum"](nodes, [v/normalizer for v in flavor], invariants, reference)
    sixth = sum(w*model["time_kernel"](x, 6) for w, x in zip(weights.entries(), nodes, strict=True))
    difference = sixth - times[6]
    print("算术概率", *[w.str(12) for w in weights.entries()])
    print("第六项与参考读数之差", difference.str(15))

    return difference


def verify():
    difference = verify_at(scale, indices)
    assert arb("1.33072788423e-8") < difference < arb("1.33072788426e-8")


if __name__ == "__main__":
    verify()

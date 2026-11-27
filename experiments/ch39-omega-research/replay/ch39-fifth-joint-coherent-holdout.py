"""在第六与第七壳直接读取第五壳完整共同见证的固定十九站算子。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def exclude_radius(amplitudes, weights, field, k, w, cap, delta):
    # 对所有算子范数距离至多 delta 的候选，在相同谱线上包围热矩导数。
    lip_k, lip_w = arb(0), arb(0)
    for amplitude, c in zip(amplitudes, weights, strict=True):
        left = max(arb(0), (amplitude.lower()-delta).lower())
        right = min(cap.upper(), (amplitude.upper()+delta).upper())
        sk, sw = arb(0), arb(0)
        for g, alpha, z in field:
            if z.is_zero():
                kernel = arb(1)
            else:
                lower = (-z/left**2).exp() if left > 0 else arb(0)
                kernel = lower.union((-z/right**2).exp())
            sk += alpha*kernel
            sw += g*kernel
        lip_k += abs(c)*right*abs(sk)/arb.pi()
        lip_w += abs(c)*right**3*abs(sw)/(2*arb.pi())
    assert abs(k-1) > delta*lip_k or abs(w-arb("1e-122")) > delta*lip_w


def colour_probe(root, number, bands, weights):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    gauge = lib("all-field-signed-crossing")
    data = gauge["joint_system"](root, "1.176397e14", number=number, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, number)
    ctx.prec = 192
    mp = arb("1.22089e19")
    rho = (arb("246.22") / mp)**2 / 2
    species, _ = lib("full-neutral-free-cutoff")["response_rows"](
        data, masses, signs, rho, arb(".42995401872")
    )
    p0, p = data["p0"], (arb(560000) / mp)**2
    log = (p / p0).log()
    rule = gauge["angular_rule"](15, 42)
    total = data["initial"][2]
    for ends_x, c in zip(bands, weights, strict=True):
        for kind, z, (_, _, beta) in species:
            if beta == 0:
                continue
            # 积分中的 exp(-(z+p h)/x) 随 x 单调增。按每个场的两个端点包围。
            ends = [gauge["angular_exposure"](kind, z, e, p0, p, log, rule)
                    for e in ends_x]
            total -= c * beta * ends[0].union(ends[1]) / (4 * arb.pi())
    return total


def readout(root, fifth_basis, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    ctx.prec = 320
    numbers = (2_306_888_100_091, 2_306_888_100_092)
    _, heat, fields = lib("three-index-free-atoms")["broken_fields"](root, numbers)
    shell = lib("three-ray-chiral-boundary")
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    native = lib("native-rh-ray-pvm")
    low, heavy, highest = levels
    assert 0 < low < heavy < highest
    gram = columns.conjugate().transpose() * columns
    assert all(gram[i, j].contains(int(i == j)) for i in range(2) for j in range(2))
    output = []
    for number, field, radius in zip(numbers, fields, (".02", ".04"), strict=True):
        _, nodes, pivots, _ = next(row for row in native["shell_weights"](root, number) if row[0] == 6)
        weights = [w.real for w in pivots[2]]
        assert sum(weights).contains(1)
        energies, basis = spectrum(shell, 6, number)
        assert all(u.overlaps(v) for u, v in zip(energies, nodes, strict=True))
        projections = (basis.conjugate().transpose() * fifth_basis) * columns
        amplitudes = [low + (heavy-low)*abs(projections[i, 0])**2
                      + (highest-low)*abs(projections[i, 1])**2 for i in range(19)]
        assert (sum(amplitudes) - (17*low + heavy + highest)).contains(0)
        bands = [(max(a.lower(), low.lower())**2, min(a.upper(), highest.upper())**2)
                 for a in amplitudes]
        assert all(l > 0 and r > 0 for l, r in bands)
        # 逐物种的 x E₂(z/x) 和 x² E₃(z/x) 在正轴递增。分开保存两个端点。
        columns_kw = []
        for left, right in bands:
            kc, wc = arb(0), arb(0)
            for g, alpha, z in field:
                k_band = (left*heat(2, z/left)).union(right*heat(2, z/right))
                w_band = (left**2*heat(3, z/left)).union(right**2*heat(3, z/right))
                kc += alpha*k_band
                wc -= g*w_band
            columns_kw.append((kc/(2*arb.pi()), wc/(4*arb.pi())))
        k, w = [sum(c*v for c, v in zip(weights, values, strict=True))
                for values in zip(*columns_kw, strict=True)]
        exclude_radius(amplitudes, weights, field, k, w, highest, arb(radius))
        print("同时恢复引力双矩需要的算子范数改变量严格大于", radius, flush=True)
        colour = colour_probe(root, number, bands, weights)
        print("零点编号、最小平方截止、牛顿项、真空项、56 万 GeV 色项",
              number, min(l.lower() for l, r in bands).str(14),
              *[v.str(16) for v in (k, w, colour)], flush=True)
        output.append((k, w, colour))
    published = ((("2.28104", "2.28106"), ("-1.398475", "-1.398471"), ("20.454856", "20.454858")),
                 (("5.19898", "5.19900"), ("-3.35094", "-3.35092"), ("18.343919", "18.343921")))
    for values, bounds in zip(output, published, strict=True):
        assert all(arb(l) < v < arb(r) for v, (l, r) in zip(values, bounds, strict=True))
    return output


def verify(root):
    joint = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    operator = joint["verify"](root, return_operator="spectral")
    return readout(root, *operator)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

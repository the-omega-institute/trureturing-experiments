"""允许跨越平凡字符子空间时，一次同谱输运可连续满足七壳的弱差。"""

from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def experiment(root, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    geo = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    frames, _ = lib("cyclic-two-characters-nine")["characters"](geo["geo"]["sites"])
    gram = native["gram_cholesky"](root)
    previous = lib("six-cyclic-weak-horizon")
    adj, squared = previous["adj"], previous["squared"]
    centers = list(map(arb, (*previous["Q"], "506.579672393")))
    ctx.prec = 256
    bases = []
    for j in range(7):
        _, v = lib("adjacent-pvm-transport")["spectrum"](geo, 6, 2_306_888_100_092 + j)
        weights = [
            z.real for z in native["readout"](v, range(19), [14, 8, 7, 6], gram[3], 3)
        ]
        bases.append((adj(v), weights))
    b = acb_mat([[basis[0, k] for k in range(19)] for basis, _ in bases]) * frames[0]
    inv = b.inv()
    # 选取 B 的实际相位；同一复相位既出现在 B 中也出现在指定幅度中。
    # 以坐标框中每一行最强分量将低线统一定相，避免依赖本征向量例程的相位。
    phases = []
    for j in range(7):
        k = max(range(7), key=lambda k: abs(b[j, k]).lower())
        phases.append((-1) ** j * b[j, k] / abs(b[j, k]))
    width = arb(".00005")
    mp, mz = arb("1.22089e19"), arb("91.1876")
    floor, first, second = levels
    y = acb_mat(
        [
            [phase * ((q + arb(0, width.upper())) / mp - floor).sqrt()]
            for phase, q in zip(phases, centers, strict=True)
        ]
    )
    lift = inv * y
    mass = squared(lift)
    assert mass < first - floor
    one = acb_mat([[arb(1 if k == 0 else 0)] for k in range(6)])
    f = frames[1] * one
    g = frames[2] * one
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0 = (mz / mp) ** 2
    p1 = (arb("1.176397e14") / mp) ** 2
    assert 0 < mass < arb(".06")
    minimum = arb(1)
    margin = arb(1)
    numbers = range(2_306_888_100_092, 2_306_888_100_099)
    fields = lib("three-index-free-atoms")["broken_fields"](root, numbers)[2]
    vacuum_bounds = []
    for j, ((basis, w), field) in enumerate(zip(bases, fields, strict=True)):
        cutoffs = []
        for i in range(1, 19):
            row = acb_mat([[basis[i, k] for k in range(19)]])
            membership = [
                k for k, E in enumerate(frames) if squared(row * E) > arb(".9")
            ]
            assert len(membership) == 1
            if membership[0] == 0:
                a = floor + abs((row * frames[0] * lift)[0, 0]) ** 2
            elif membership[0] == 1:
                a = floor + (first - floor - mass) * abs((row * f)[0, 0]) ** 2
            else:
                a = floor + (second - floor) * abs((row * g)[0, 0]) ** 2
            cutoffs.append(a)
        low = min(a.lower() for a in cutoffs)
        all_cutoffs = [(centers[j] + arb(0, width.upper())) / mp, *cutoffs]
        coefficient = -sum(g for g, _, _ in field) / (8 * arb.pi())
        mass_bound = sum(abs(g) * z for g, _, z in field) / (4 * arb.pi())
        vacuum = sum(
            wi * (coefficient * a**4 + arb(0, (mass_bound * a**2).upper()))
            for wi, a in zip(w, all_cutoffs, strict=True)
        )
        assert vacuum < arb("-.55")
        vacuum_bounds.append(vacuum.upper())
        assert (
            squared(acb_mat([[basis[0, k] for k in range(19)]]) * frames[0]) > arb(".9")
            and low**2 > p0
        )
        tail = (large(low**2, p0, p1) - ref).abs_upper() * sum(
            abs(wi).upper() for wi in w[1:]
        )
        faces = [
            ref
            + w[0] * (weak(((centers[j] + sign * width) / mz) ** 2) - ref)
            + arb(0, tail)
            for sign in (-1, 1)
        ]
        assert faces[0] > 0 > faces[1], (j, faces)
        minimum = min(minimum, low)
        margin = min(margin, faces[0].lower(), (-faces[1]).lower())
    print("整个七壳校准盒的真空项统一上界", max(vacuum_bounds).str(12), flush=True)
    assert minimum > arb("1e-8") and margin > arb("1e-8")
    print(
        "七壳共同弱差的投影范数平方、其他截止下界及面余量",
        mass.str(12),
        minimum.str(12),
        margin.str(12),
        flush=True,
    )


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, source["verify"](root, return_operator="spectral")[2])

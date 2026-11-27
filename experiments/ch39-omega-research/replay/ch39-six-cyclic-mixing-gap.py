"""保持平凡字符子空间的六壳共同解必须有非零的定量字符混合。"""

from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def load_replay(root, name):
    return run_path(str(root / ("ch39-" + name + ".py")))


def leakage_error(levels, traces, shells):
    main = 1 if traces[1] > traces[2] else 2
    tiny = 3 - main
    assert traces[main] > traces[tiny]
    beta, alpha = traces[main].sqrt(), traces[tiny].sqrt()
    t = arb(".000033")
    # 以 U 的两个六维对角块的酉极分解近似。平凡块按原样保留。
    displacement = beta * t**2 / (1 + (1 - t**2).sqrt()) + alpha * t
    cutoff_l1 = (2 * beta + displacement) * displacement
    large_cutoff = levels[0] + (beta + displacement) ** 2
    small_cutoff = levels[0] + (alpha + t * beta) ** 2
    errors = []
    for frames, weights, field in shells:
        groups = [
            [
                i
                for i in range(19)
                if sum(abs(f[i, k]) ** 2 for k in range(f.ncols())) > arb(".9")
            ]
            for f in frames
        ]
        assert list(map(len, groups)) == [7, 6, 6]
        assert len(set().union(*map(set, groups))) == 19
        main_weight = max(abs(weights[i]).upper() for i in groups[main])
        small_weight = sum(abs(weights[i]) for i in groups[tiny])
        counts, coefficients, masses = zip(*field, strict=True)
        # |1-exp(-x)| ≤ min(1,x) ≤ sqrt(x)。
        newton_slope = (
            large_cutoff * abs(sum(coefficients))
            + sum(abs(c) * m.sqrt() for c, m in zip(coefficients, masses, strict=True))
        ) / arb.pi()
        vacuum_slope = (
            large_cutoff**3 * abs(sum(counts))
            + large_cutoff
            * sum(abs(g) * m for g, m in zip(counts, masses, strict=True))
        ) / (2 * arb.pi())
        newton = (
            main_weight * newton_slope * cutoff_l1
            + small_weight * small_cutoff**2 * sum(map(abs, coefficients)) / arb.pi()
        )
        vacuum = (
            main_weight * vacuum_slope * cutoff_l1
            + small_weight * small_cutoff** 4 * sum(map(abs, counts)) / (4 * arb.pi())
        )
        errors.append((newton.upper(), vacuum.upper()))
    return [max(pair[i] for pair in errors) for i in (0, 1)]


def experiment(root, fifth, columns, levels):
    graph = load_replay(root, "three-ray-chiral-boundary")
    frames, _ = load_replay(root, "cyclic-two-characters-nine")["characters"](
        graph["geo"]["sites"]
    )
    sieve = load_replay(root, "six-cyclic-line-sieve")
    adj = sieve["adj"]
    floor, heavy, top = levels
    source = (
        fifth
        * columns
        * acb_mat([[(heavy - floor).sqrt(), 0], [0, (top - floor).sqrt()]])
    )
    traces = [sieve["squared"](adj(f) * source) for f in frames]
    native = load_replay(root, "native-rh-ray-pvm")
    gram = native["gram_cholesky"](root)
    numbers = range(2_306_888_100_092, 2_306_888_100_098)
    fields = load_replay(root, "three-index-free-atoms")["broken_fields"](
        root, numbers
    )[2]
    shells = []
    for number, field in zip(numbers, fields, strict=True):
        _, basis = load_replay(root, "adjacent-pvm-transport")["spectrum"](
            graph, 6, number
        )
        weights = [
            w.real
            for w in native["readout"](basis, range(19), [14, 8, 7, 6], gram[3], 3)
        ]
        shells.append(([adj(basis) * f for f in frames], weights, field))
    newton, vacuum = leakage_error(levels, traces, shells)
    assert newton < arb("2.20e-7") and vacuum < arb("9.42e-7")
    print(
        "混合范数不大于 0.000033 时比较输运的牛顿与真空绝对误差上界",
        newton.str(12),
        vacuum.str(12),
        flush=True,
    )


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = load_replay(root, "fifth-full-positive-joint")["verify"](
        root, return_operator="spectral"
    )
    experiment(root, *source)

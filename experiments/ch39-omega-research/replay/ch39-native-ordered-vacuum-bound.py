"""十九站固定实交叉权的保序截止，在牛顿标定后对真空项的上界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    spectra = run_path(str(root / "ch39-three-index-free-atoms.py"))
    _, heat, (fields,) = spectra["broken_fields"](root, [2_306_888_100_085])
    ctx.prec = 256
    alpha = -sum(a for _, a, _ in fields) / (2 * arb.pi())
    beta = -sum(g for g, _, _ in fields) / (8 * arb.pi())
    delta_k = sum(abs(a) * (1 - heat(2, z)) for _, a, z in fields) / (2 * arb.pi())
    # W'(x) = -x sum(g exp(-z/x))/(4 pi) 且 x(1-exp(-z/x)) <= z。
    derivative_error = sum(abs(g) * z for g, _, z in fields) / (4 * arb.pi())
    assert alpha > 0 and beta > 0
    excluded = {}
    for charge, _, direct, mirror in native["shell_weights"](root):
        for pivot, (source, reflected) in enumerate(zip(direct, mirror, strict=True), 1):
            c = [v.real for v in source]
            assert all(w.overlaps(z.real) for w, z in zip(c, reflected, strict=True))
            tv = sum(abs(w) for w in c)
            a = (1 - tv * delta_k) / alpha
            tails = [sum(c[j:]) for j in range(1, len(c))]
            depth = max(-t.lower() for t in tails)
            ceiling = max(arb(0), max(t.upper() for t in tails))
            assert a > 0 and depth > 0
            # 有序阶梯 T 满足 ∫|T|=2∫T₊−∫T，其中 −∫T 也被牛顿标定所限。
            absolute_tail_integral = 2 * (1 + ceiling) + (1 + tv * delta_k) / alpha
            if ceiling <= 1:
                # 负尾占最小 u、正尾占最大 u 时二阶矩最高。这也覆盖首层以上的零尾。
                second_moment = ceiling - (ceiling + a)**2 / (depth + ceiling)
            else:
                second_moment = ceiling - a**2 / depth
            bound = beta * second_moment + derivative_error * absolute_tail_integral
            print(charge, pivot, "后缀上界、最大深度、W 上界",
                  *[v.upper().str(12) for v in (ceiling, depth, bound)], flush=True)
            if bound < 0:
                excluded[charge, pivot] = bound
            # 两位置若高端总权非负则一阶矩非负，若为负则非平凡后缀全非正。
            assert tv * delta_k < 1
            two_site_bound = -beta * a**2/depth + derivative_error * (2 + (1+tv*delta_k)/alpha)
            assert two_site_bound < 0
    expected = {(q, 2) for q in (1, -4, 2, -3, 6)} | {
        (1, 1), (2, 1), (-3, 1), (-4, 3), (2, 3), (6, 3)
    }
    assert set(excluded) == expected
    assert excluded[2, 2] < arb("-.090") and excluded[2, 1] < 0
    return excluded


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

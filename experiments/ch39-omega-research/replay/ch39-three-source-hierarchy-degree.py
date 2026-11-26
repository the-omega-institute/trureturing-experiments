"""三份 Gram 源上，非负系数多项式重现最轻与最重上型质量比的最低次数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify(root):
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    ray = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ctx.prec = 320
    frame, _ = flavor["family_space"](native["gram_cholesky"](root))
    sites = [ray["soft"]["pole"], *ray["boundary_sites"]()]
    assert sites == [14, 8, 7, 6]
    height = acb.zeta_zero(2_306_888_100_085).imag
    up = ray["laplacian"](-4, height)
    levels, vectors = up.eig(right=True, algorithm="rump")
    assert all(value.imag.contains(0) for value in levels)
    order = sorted(range(19), key=lambda i: levels[i].real.mid())
    spectrum = [levels[i].real for i in order]
    assert spectrum[0] > 0 and all(a < b for a, b in zip(spectrum, spectrum[1:]))

    # 每个十九站谱投影先取四个站位，再压缩至同一个三源空间。
    adj = flavor["adj"]
    projections = []
    for index in order:
        column = acb_mat([[vectors[site, index]] for site in sites])
        vector = adj(frame) * column
        norm_sq = sum(abs(vectors[k, index])**2 for k in range(19))
        assert norm_sq > 0
        projections.append((vector * adj(vector)) * (1 / norm_sq))
    zero = acb_mat(3, 3)
    assert (sum(projections, zero) - flavor["identity"](3)).contains(zero)

    # 取同一单位 L = λmax(H₋₄)，最高谱点的归一值严格等于一。
    normalized = [level / spectrum[-1] for level in spectrum[:-1]] + [arb(1)]
    weights = list(normalized)
    mass_ratios = [arb(1)]
    for degree in range(1, 1953):
        matrix = sum((proj * weight for proj, weight in zip(projections, weights, strict=True)), zero)
        eigenvalues = matrix.eig(algorithm="rump")
        assert all(value.imag.contains(0) for value in eigenvalues)
        low, middle, high = sorted((value.real for value in eigenvalues), key=lambda value: value.mid())
        assert 0 < low < middle < high
        mass_ratios.append((low / high).sqrt())
        weights = [weight * step for weight, step in zip(weights, normalized, strict=True)]

    lower = arb(".00000601") / arb(".8659")
    center = arb(".00000615") / arb(".8616")
    upper = arb(".00000629") / arb(".8573")
    assert min(v.lower() for v in mass_ratios[:3]) > arb(".90760897")
    assert min(v.lower() for v in mass_ratios[:1947]) > arb("7.37949e-6") > upper
    assert lower < center < mass_ratios[1947] < upper
    assert min(v.lower() for v in mass_ratios[:1952]) > arb("7.16467e-6") > center
    assert mass_ratios[1952] < center

    print("参考窗下端、中心与上端", lower.str(15), center.str(15), upper.str(15))
    for degree in (1, 2, 1946, 1947, 1951, 1952):
        print("单项式次数与上型最轻最重质量比", degree, mass_ratios[degree].str(15))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

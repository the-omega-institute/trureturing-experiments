"""冻结三壳训练值，求所有随距离非增的非负加权平均在外侧夸克谱上的必要范围。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def prefix_bounds(energy, records):
    ranked = sorted(((abs(energy - x), y) for x, y in records), key=lambda row: row[0].mid())
    assert ranked[0][0] > 0
    assert all(a[0] < b[0] for a, b in zip(ranked, ranked[1:]))
    means = []
    total = arb(0)
    for k, (_, y) in enumerate(ranked, start=1):
        total += y
        means.append(total/k)
    low = min(y.lower() for y in means).union(min(y.upper() for y in means))
    high = max(y.lower() for y in means).union(max(y.upper() for y in means))
    return low, high


def charge_records(records, center, ordered, shell):
    labels = []
    for number in range(center - 1, center + 2):
        height = acb.zeta_zero(number).imag
        for charge in (-4, 2, 6):
            labels.extend((t, charge) for t in ordered(shell, charge, height))
    labels.sort(key=lambda item: item[0].mid())
    assert len(records) == 171
    paired = list(zip(labels, records, strict=True))
    assert all(t.overlaps(v) for (t, _), (v, _) in paired)
    assert all(a[0] < b[0] for a, b in zip(labels, labels[1:]))
    same = {q: [row for (_, label), row in paired if label == q] for q in (-4, 2)}
    assert all(len(rows) == 57 for rows in same.values())
    return same


def mass_bounds(projectors, bounds):
    matrices = [
        sum((p * w for p, w in zip(projectors, endpoint, strict=True)), acb_mat(4, 4))
        for endpoint in zip(*bounds, strict=True)
    ]
    spectra = []
    for matrix in matrices:
        eigenvalues = matrix.eig(multiple=True, algorithm="rump")
        assert len(eigenvalues) == 4 and all(z.imag.contains(0) for z in eigenvalues)
        left = sorted(z.real.lower() for z in eigenvalues)
        right = sorted(z.real.upper() for z in eigenvalues)
        spectra.append([a.union(b) for a, b in zip(left, right, strict=True)])
    return list(zip(spectra[0][:3], spectra[1][1:], strict=True))


def free_fraction(low, high, target):
    a, b = map(arb, target)
    if low > b:
        return (1 - b / low).lower()
    if high < a:
        return ((a - high) / (1 - high)).lower()
    return arb(0)


def verify(root):
    model = run_path(str(root / "ch39-three-shell-outer-holdout.py"))
    training = run_path(str(root / "ch39-three-shell-two-space-readout.py"))
    records, _ = training["verify"](root, return_data=True)
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))["native"]
    center = model["CENTER"]
    same = charge_records(records, center, ordered, shell)
    # 每个整数都是该荷中单独就被排除的质量升序位置（从一数起）。
    excluded = {
        (-2, -4, "all"): {3}, (-2, 2, "all"): {2, 3},
        (2, -4, "all"): {2, 3}, (2, 2, "all"): {1, 2},
        (-2, -4, "same"): {1, 3}, (2, -4, "same"): {1, 2, 3},
        (-2, -4, "mixed"): {3}, (2, -4, "mixed"): {3},
    }
    for offset in model["OFFSETS"]:
        number = center + offset
        projectors = native(root, number)
        ctx.prec = 320
        height = acb.zeta_zero(number).imag
        for charge, rows, reference in zip((-4, 2), projectors, model["WINDOWS"], strict=True):
            nodes = ordered(shell, charge, height)
            all_bounds = [prefix_bounds(t, records) for t in nodes]
            same_bounds = [prefix_bounds(t, same[charge]) for t in nodes]
            # 各层还可独立混合两个平均值，取两个区间的凸包。
            mixed_bounds = [(a.min(c), b.max(d)) for (a, b), (c, d)
                            in zip(all_bounds, same_bounds, strict=True)]
            for pool, bounds in (("all", all_bounds), ("same", same_bounds), ("mixed", mixed_bounds)):
                masses = mass_bounds(rows, bounds)
                failed = set()
                for i, ((low, high), target) in enumerate(zip(masses, reference, strict=True), 1):
                    fraction = free_fraction(low, high, target)
                    if fraction > 0:
                        failed.add(i)
                    else:
                        # 真实可达区间的两端均由选取三个本征向量张成的空间取得。
                        assert low < arb(target[0]) and high > arb(target[1])
                    print(number, pool, charge, i, low.str(13), high.str(13),
                          fraction.str(13), flush=True)
                    if pool == "all" and charge == 2 and (offset, i) == (-2, 3):
                        assert low > arb(".0722179") and fraction > arb(".80475")
                    if pool == "all" and charge == 2 and (offset, i) == (2, 1):
                        assert low > arb(".0001396976") and fraction > arb(".90193")
                    if pool == "all" and charge == -4 and i == 3:
                        assert high < arb(".80930" if offset == -2 else ".75824")
                    if pool == "same" and charge == -4 and i == 1:
                        assert fraction > arb(".99315" if offset == -2 else ".99359")
                    if pool == "mixed" and charge == -4 and i == 3:
                        left, right, margin = (".8272257", ".8272258", ".17406") if offset == -2 else (
                            ".7856137", ".7856138", ".33437")
                        assert arb(left) < high < arb(right) and fraction > arb(margin)
                assert excluded.get((offset, charge, pool), set()) == failed


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

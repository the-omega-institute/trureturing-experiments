"""固定第五壳源，多壳零弱差要求站位酉的支撑触及多少个三循环轨道。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def cutoff_bound(rows, source, support, floor):
    count = len(rows)
    restricted = [[r[k] for k in support] for r in rows]
    minors = [(-1)**j * acb_mat(restricted[:j]+restricted[j+1:]).det() for j in range(count)]
    denominator = sum(abs(z) for z in minors)
    assert denominator > 0
    # 支撑内的分量依行列式恒等式精确为零，先消去，再与原始两列相乘。
    outside = [k for k in range(7) if k not in support]
    combination = [sum(c*r[k] for c, r in zip(minors, rows, strict=True)) for k in outside]
    energies = [sum(z*source[k,h] for k,z in zip(outside, combination, strict=True)) for h in (0,1)]
    numerator = sum(max(arb(0), abs(z).lower())**2 for z in energies)
    return ((floor + numerator/denominator**2)*arb("1.22089e19")).lower()


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    candidates, shells, _ = lib("seven-unrestricted-vacuum-weak")["prepare"](root, levels)
    graph = lib("three-ray-chiral-boundary")
    frames, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    floor, first, second = levels
    scale = acb_mat([[(first-floor).sqrt(), 0], [0, (second-floor).sqrt()]])
    source = frames[0].conjugate().transpose()*(fifth*columns*scale)
    # 前八个窗另加第七至十二壳，共九个全覆盖。列出严格向下舍入的 GeV 下界。
    windows = [(0,2,"1.99e14"), (0,3,"1.66e13"), (0,4,"6.68e9"), (0,5,"5.7e5"),
               (0,6,"14000"), (1,2,"1.13e14"), (1,3,"9.23e11"), (1,4,"9.6e9"),
               (1,5,"6.2e6"), (0,7,None)]
    for start, size, threshold in windows:
        choices = candidates[start:start+size]
        raw = [s[0].tolist() for s in shells[start:start+size]]
        count = 0
        minimum = arb("1e19")
        unresolved = []
        for ids in product(*choices):
            rows = [r[i] for r,i in zip(raw, ids, strict=True)]
            for support in combinations(range(7), size-1):
                lower = cutoff_bound(rows, source, support, floor)
                assert lower.is_finite()
                minimum = min(minimum, lower)
                count += 1
                if not lower > arb(3000):
                    unresolved.append((ids, tuple(k for k in range(7) if k not in support)))
        if threshold is not None:
            assert not unresolved and minimum > arb(threshold)
        else:
            assert set(unresolved) == {((0,)*7, (k,)) for k in (3,4,6)}
        print("起始壳、壳数、组合、共同截止下界 GeV、未分辨数",
              start+7, size, count, minimum.str(10), len(unresolved), flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))

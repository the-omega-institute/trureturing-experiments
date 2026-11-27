"""在原生第三枢轴选出的同一组站位上限制五荷三维轻带。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def least_singular(rows):
    gram = rows.conjugate().transpose()*rows
    eigenvalues = gram.eig(algorithm="rump")
    assert len(eigenvalues) == 3
    assert all(z.is_finite() and z.imag.contains(0) and z.real > 0 for z in eigenvalues)
    return min(z.real.sqrt().lower() for z in eigenvalues)


def uniform_bound(graph, sites, pump):
    bases = list(pump["rotation_bases"](graph))
    reference = []
    for charge in (1, -4, 2, -3, 6):
        base = pump["clock_matrices"](graph, charge)[0]
        columns = []
        for frame in bases:
            block = frame.conjugate().transpose()*base*frame
            vector, _ = pump["first_line"](block)
            columns.append((frame*acb_mat([[z] for z in vector])).entries())
        rows = acb_mat([[v[i] for v in columns] for i in sites])
        reference.append(least_singular(rows))
    m = arb(".55446")
    assert min(reference) > m
    # 既有全环面证书给当前轻线对参考线的主角正弦上界 1/3。
    lower = m*(arb(8)/9).sqrt()-(1-m*m).sqrt()/3
    assert lower > arb(".24534")
    print("整环面参考空间、实际轻带的统一最小奇异值下界", m, lower.str(13), flush=True)


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    graph = lib("three-ray-chiral-boundary")
    sites = lib("terminal-port-support")["verify"](root)
    ctx.prec = 160
    pump = lib("native-light-projector-pump")
    uniform_bound(graph, sites, pump)
    frames = list(pump["rotation_bases"](graph))
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    results = []
    for shell in (0, *range(7, 20)):
        singular = []
        for charge in (1, -4, 2, -3, 6):
            levels, basis = spectrum(graph, charge, 2_306_888_100_085+shell)
            assert levels[2] < arb("8.685") < 9 < arb("9.6") < levels[3]
            if shell == 0:
                light = acb_mat([[basis[i, j] for j in range(3)] for i in range(19)])
                for frame in frames:
                    projection = frame.conjugate().transpose()*light
                    multiplicity = sum(abs(z)**2 for z in projection.entries())
                    assert arb(".9") < multiplicity < arb("1.1")
            rows = acb_mat([[basis[i, j] for j in range(3)] for i in sites])
            sigma = least_singular(rows)
            singular.append(sigma)
            results.append((sigma, shell, charge))
        if shell == 0:
            assert min(singular) > arb(".57369")
        print("壳偏移、按五荷顺序的最小奇异值下界", shell, *(v.str(13) for v in singular), flush=True)
    assert len(results) == 70
    least, shell, charge = min(results, key=lambda item: item[0])
    assert least > arb(".53219")
    print("七十组的共同下界及取到最小包围下端的壳偏移、荷", least.str(15), shell, charge)
    return least


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

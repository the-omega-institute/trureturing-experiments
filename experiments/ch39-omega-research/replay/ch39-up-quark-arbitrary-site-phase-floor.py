"""右邻到右外壳在任意十九站对角相位下的上夸克 Yukawa 下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CENTER = 2_306_888_100_085
ORBITS = ((0, 3, 4), (1, 2, 5), (6, 7, 8), (9, 15, 16), (10, 13, 17),
          (11, 12, 18), (14,))
DUALS = ((1905, 1919, -1136, -4724, -5113, 1919, -1001),
         (1779, -3316, -19, 239, 232, -96, -7),
         (472, 124, -10, 208, 211, -1491, -2))


def adj(m):
    return m.conjugate().transpose()


def positive(matrix):
    """已知矩阵为厄米，以顺序主子式严格证其正定。"""
    for size in range(1, matrix.nrows()+1):
        minor = acb_mat([[matrix[i, j] for j in range(size)] for i in range(size)]).det()
        assert minor.imag.contains(0) and minor.real > 0


def verify(root):
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    source = CENTER + 1
    _, responses, _ = training(root, ((source, ".0294"),), ".1", (),
                               ((6, 7, 10, 13), (6, 7, 10, 12)))
    weights = responses[source, -4]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    ctx.prec = 384
    _, previous = spectrum(shell, -4, source)
    _, current = spectrum(shell, -4, CENTER+2)
    weighted = acb_mat([[previous[k, i]*weights[i] for i in range(19)] for k in range(19)])
    kernel = weighted*adj(previous)
    sites = shell["geo"]["sites"]
    assert sorted(k for group in ORBITS for k in group) == list(range(19))
    for group in ORBITS:
        assert all(sites.index((sites[k][1], sites[k][2], sites[k][0])) in group for k in group)

    floors = [arb(0)]*19
    for row in DUALS:
        diagonal = [arb(0)]*19
        for coefficient, group in zip(row, ORBITS, strict=True):
            for k in group:
                diagonal[k] = arb(coefficient)/10**6
        shifted = acb_mat(kernel)
        for k in range(19):
            shifted[k, k] -= diagonal[k] + arb("4e-7")
        positive(shifted)
        for j in range(19):
            bound = sum(diagonal[k]*abs(current[k, j])**2 for k in range(19))
            floors[j] = max(floors[j], bound.lower())

    keep = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert keep == (14, 8, 7, 6)
    four = acb_mat([[current[k, j] for j in range(19)] for k in keep])
    weighted = acb_mat([[four[k, j]*floors[j] for j in range(19)] for k in range(4)])
    raw = weighted*adj(four)
    threshold = arb(".0000107940")
    shifted = acb_mat(raw)
    for k in range(4):
        shifted[k, k] -= threshold
    positive(shifted)
    assert threshold > arb(".00000629")
    spectrum_floor = sorted((z.real for z in raw.eig(algorithm="rump")), key=lambda x: x.mid())
    print("四站下界矩阵的本征值", [v.str(14) for v in spectrum_floor], flush=True)
    print("十九个任意站位相位及任意三平面上的严格下界", threshold.str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

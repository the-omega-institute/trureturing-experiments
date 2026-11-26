"""两荷共用一个实多项式，在三源空间实现给定的六夸克汤川强度和 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def columns(*vectors):
    return acb_mat([[v[j] for v in vectors] for j in range(len(vectors[0]))])


def adapted_frame(gram):
    normal = [arb(1)]
    for i in (1, 2, 3):
        normal.append(-sum(gram[i][k] * normal[k] for k in range(i)) / gram[i][i])
    length = sum(v*v for v in normal).sqrt()
    n = [v / length for v in normal]
    mean = sum(n[1:]) / 3
    symmetric = [n[0], mean, mean, mean]
    standard = [v-w for v, w in zip(n, symmetric)]
    a = sum(v*v for v in symmetric).sqrt()
    b = sum(v*v for v in standard).sqrt()
    assert a > 0 and b > 0
    sqrt3 = arb(3).sqrt()
    x = [sqrt3 * mean / a] + [-n[0] / (sqrt3*a)] * 3
    y = [arb(0)] + [(standard[j % 3 + 1] - standard[(j+1) % 3 + 1]) / (sqrt3*b)
                    for j in range(1, 4)]
    z = [b*t/a-a*w/b for t, w in zip(symmetric, standard)]
    frame = columns(x, y, z)
    assert (adj(frame)*frame).contains(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
    assert (acb_mat([normal]) * frame).contains(acb_mat(1, 3))
    return frame


def coordinates(matrix):
    assert matrix[0, 1].contains(0) and matrix[1, 2].real.contains(0)
    return [matrix[j, j].real for j in range(3)] + [
        matrix[0, 2].real, matrix[0, 2].imag, matrix[1, 2].imag]


def target_blocks(ckm, flavor):
    up = list(map(arb, (".00000615", ".00311", ".8616")))
    down = list(map(arb, (".0000135", ".000268", ".01401")))
    v = acb_mat(ckm["standard_matrix"](*ckm["centers"]))

    def diagonal(s):
        return acb_mat([[s[i]**2 if i == j else 0 for j in range(3)] for i in range(3)])

    a, b = adj(v) * diagonal(up) * v, diagonal(down)
    residual_norm = (abs(a[1, 0])**2 + abs(a[2, 0])**2).sqrt()
    assert residual_norm > 0
    z = [acb(0)] + [acb(0, 1) * a[j, 0] / residual_norm for j in (1, 2)]
    x = [acb(0), -z[2].conjugate(), z[1].conjugate()]
    unitary = columns(x, [acb(1), acb(0), acb(0)], z)
    assert (adj(unitary) * unitary).contains(diagonal([1, 1, 1]))
    targets = [adj(unitary) * matrix * unitary for matrix in (a, b)]
    values_u, pu = flavor["projections"](targets[0])
    values_d, pd = flavor["projections"](targets[1])
    for value, mass in zip(values_u + values_d, up + down, strict=True):
        assert value > 0 and value.overlaps(mass**2)
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(v[i,j])**2)
               for i in range(3) for j in range(3))
    jarlskog = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    reference = (v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag
    assert jarlskog > 0 and jarlskog.overlaps(reference)
    print("中心 CKM 的有符号 J", jarlskog.str(12))
    return targets


def verify(root):
    ray = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    ctx.prec = 320
    frame = adapted_frame(native["gram_cholesky"](root))
    sites = [ray["soft"]["pole"], *ray["boundary_sites"]()]
    assert sites == [14, 8, 7, 6]
    labels = ray["geo"]["sites"]
    rotate = [labels.index((y,z,x)) for x,y,z in labels]
    assert rotate[sites[0]] == sites[0] and set(rotate[j] for j in sites[1:]) == set(sites[1:])
    tuned = set(ray["soft"]["star"]["orbit"])
    assert {rotate[j] for j in tuned} == tuned
    # 带额外时钟相位的有向边，在坐标三循环下保持端点顺序。
    assert all(rotate[i] < rotate[j] for i,j,d,_ in ray["geo"]["edges"] if d in (724, 804))
    height = acb.zeta_zero(2_306_888_100_085).imag
    identity = acb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    blocks, spectra = [], []
    for charge in (-4, 2):
        matrix = ray["laplacian"](charge, height)
        scaled = (matrix - identity*18) * (arb(1)/8)
        levels = scaled.eig(algorithm="rump")
        assert all(v.imag.contains(0) for v in levels)
        spectrum = sorted((v.real for v in levels), key=lambda v: v.mid())
        assert all(a < b for a, b in zip(spectrum, spectrum[1:]))
        spectra.append(spectrum)
        power = identity
        rows = []
        for degree in range(12):
            block = acb_mat([[power[i, j] for j in sites] for i in sites])
            rows.append(coordinates(adj(frame) * block * frame))
            power = power * scaled
        blocks.append(arb_mat(rows).transpose())
    system = arb_mat(blocks[0].tolist()+blocks[1].tolist())
    determinant = system.det()
    assert arb("5.5e-40") < determinant < arb("6.5e-40")
    for block in blocks:
        assert arb_mat([row[:6] for row in block.tolist()]).det() > 0
    print("两荷共同多项式的十二维行列式", determinant.str(12))
    targets = target_blocks(ckm, flavor)
    target = arb_mat([[value] for matrix in targets for value in coordinates(matrix)])
    coeff = system.solve(target, algorithm="precond")
    assert (system*coeff-target).contains(arb_mat(12, 1))
    print("十二个系数", *(c.str(12) for c in coeff.entries()))
    for charge, spectrum in zip((-4, 2), spectra):
        weights = []
        for point in spectrum:
            value = arb(0)
            for c in reversed(coeff.entries()):
                value = value*point+c
            # 对偶系统保留同一组系数在谱点处相消的精度。
            powers = arb_mat([[point**n] for n in range(12)])
            dual = system.transpose().solve(powers, algorithm="precond")
            weights.append(value.intersection((dual.transpose() * target)[0, 0]))
        positive, negative = sum(w>0 for w in weights), sum(w<0 for w in weights)
        assert (positive, negative) == (11, 8)
        print("荷、正负谱权数、最小和最大谱权", charge, positive, negative,
              min(weights).str(12), max(weights).str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

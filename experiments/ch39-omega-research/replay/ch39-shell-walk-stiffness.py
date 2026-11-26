"""用分段谱矩判定任意步数的规范刚度，并计算有限角度的面权。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def spectral_order(spectra, weights, floor):
    nodes = sorted(((v, a) for q, a in weights.items() for v in spectra[q]),
                   key=lambda pair: float(pair[0]))
    assert all(left[0] < right[0] for left, right in zip(nodes, nodes[1:]))
    assert nodes[0][1] > 0 and nodes[-1][1] < 0
    # On the first and last intervals a single squared distance suffices,
    # using exact cancellation of moments 0,1,2 for the last interval.
    for i in range(1, len(nodes) - 2):
        left, right = nodes[i][0], nodes[i + 1][0]
        active = [(left - v, a) for v, a in nodes[:i]]
        first = sum(a * d * d for d, a in active)
        middle = first + (right - left) * sum(a * d for d, a in active)
        last = sum(a * (right - v) * (right - v) for v, a in nodes[:i + 1])
        assert min(first, middle, last) > arb(floor)
    # The three-mode band has equal masses (ranks). Its first positive-part
    # moment stays strictly positive by linear interpolation between knots.
    band = [(v, a) for v, a in nodes if v < 9]
    assert len(band) == 3 * len(weights) and sum(a for _, a in band) == 0
    assert all(sum(a * (u - v) for v, a in band[:i]) > 0
               for i, (u, _) in enumerate(band) if i > 0)


def check_three_characters(matrix, identity, orbits):
    # The cycle commutes with the magnetic Laplacian. A negative determinant
    # for each restricted Hermitian form gives at least one level below 9.
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    shifted = matrix - 9 * identity
    for character in range(3):
        active = [orbit for orbit in orbits if len(orbit) == 3 or character == 0]
        basis = acb_mat(19, len(active))
        for col, orbit in enumerate(active):
            for k, row in enumerate(orbit):
                basis[row, col] = omega ** (k * character)
        determinant = (basis.conjugate().transpose() * shifted * basis).det()
        assert determinant.imag.contains(0) and determinant.real < 0


def trace(matrix):
    return sum(matrix[i, i].real for i in range(matrix.nrows()))


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory / "ch39-three-ray-chiral-boundary.py"))
    reps = run_path(str(directory / "ch39-native-trace-plaquette-jets.py"))
    ctx.prec = 256
    height = acb.zeta_zero(11431).imag
    identity = acb_mat(19, 19)
    for i in range(19):
        identity[i, i] = 1
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y, z, x)) for x, y, z in sites]
    orbits = []
    for i, j in enumerate(rotation):
        if i == j:
            orbits.append([i])
        elif i == min(i, j, rotation[j]):
            orbits.append([i, j, rotation[j]])
    assert sorted(i for orbit in orbits for i in orbit) == list(range(19))
    spectra, walks = {}, {}
    for q in (0, 1, 2, 3, 4, 6):
        matrix = shell["laplacian"](q, height)
        walks[q] = identity - matrix / 36
        if q:
            eig = matrix.eig()
            assert len(eig) == 19 and all(v.imag.contains(0) for v in eig)
            spectra[q] = sorted([v.real for v in eig], key=float)
            assert all(0 < v < 36 for v in spectra[q])
            assert all((matrix[rotation[i], rotation[j]] - matrix[i, j]).contains(0)
                       for i in range(19) for j in range(19))
            check_three_characters(matrix, identity, orbits)
    lower = max(values[2] for values in spectra.values())
    upper = min(values[3] for values in spectra.values())
    assert arb("8.6387") < lower < arb("8.6389") and arb("9.6728") < upper < arb("9.6730")
    assert all(values[0] > arb("1.09") for values in spectra.values())
    print("三个模式的共同截断区间两端", lower.str(13), upper.str(13))
    spectral_order(spectra, {1: -1, 2: 1, 3: -1, 4: 1}, ".025")
    spectral_order(spectra, {1: -7, 2: 1, 3: -1, 4: 4, 6: 3}, ".100")

    family, dimension = reps["family"], reps["dimension"]
    theta = shell["geo"]["theta"]
    angle = height * (theta(724) - theta(244))
    expected = {1: ("1.59882e-74", "1.59883e-74", "1.87606e-76", "1.87607e-76"),
                2: ("1.28675e-42", "1.28676e-42", "7.88481e-46", "7.88482e-46")}
    for steps in (1, 2):
        blocks = {q: b ** steps for q, b in walks.items()}
        faces = []
        for flat in (False, True):
            factors = []
            for rep in family:
                q = rep[0]
                block = blocks[0 if flat else abs(q)]
                value = (identity - block * (q * angle / 2).sin()**2).det()
                assert value.imag.contains(0) and value.real > 0
                factors.append(value.real ** dimension(rep))
            faces.append(prod(factors))
        l, h, lf, hf = map(arb, expected[steps])
        assert l < faces[0] < h and lf < faces[1] < hf
        print("步数、实际面权与平坦面权之比", steps, (faces[0] / faces[1]).str(14))
        # For s=1 the first holonomy term in -log W is from sin^6/3;
        # for s=2 it is from sin^4/2. Lower trace moments cancel exactly.
        power, divisor = (3, 192) if steps == 1 else (2, 32)
        change = sum(rep[0]**(2 * power) * dimension(rep) *
                     (trace(blocks[abs(rep[0])]**power) - trace(blocks[0]**power))
                     for rep in family) / divisor
        low, high = ("-4.66759", "-4.66758") if steps == 1 else ("-3.39187", "-3.39186")
        assert arb(low) < change < arb(high)
        print("负对数面权与平坦值的首个差项", 2 * power, change.str(14))


if __name__ == "__main__":
    verify()

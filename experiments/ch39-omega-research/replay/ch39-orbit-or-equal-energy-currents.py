"""前 256 枚零点的轨道等权或八站 Dirichlet 等能量检验。"""
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


SITES = (0, 2, 3, 5, 6, 9, 11, 13)
ORBIT = (0, 3, 4)


def check_support(shell):
    points = shell["geo"]["sites"]
    assert set(shell["soft"]["star"]["orbit"]) == set(ORBIT)
    stabilizers = [set(), set()]
    for permutation in permutations(range(3)):
        image = [points.index(tuple(point[j] for j in permutation)) for point in points]
        for group, support in zip(stabilizers, (SITES, ORBIT)):
            if {image[i] for i in support} == set(support):
                group.add(permutation)
        if permutation == (1, 2, 0):
            assert all({i, image[i], image[image[i]]} - set(SITES) for i in SITES)
    assert len(stabilizers[0]) == 2 and len(stabilizers[1]) == 3
    assert stabilizers[0] & stabilizers[1] == {(0, 1, 2)}


def averaged(green, left, right):
    return sum(green[i, j] for i in left for j in right)/(len(left)*len(right))


def check_equal_orbit(green):
    for tied in ((0, 3), ORBIT):
        vertices = [tied]+[(i,) for i in SITES if i not in tied]
        for left in vertices:
            for right in vertices:
                floor = ".0011"
                if (left, right) == ((5,), (6,)):
                    floor = "-.004"
                elif left == right and left in ((5,), (6,)):
                    floor = ".004"
                amplitude = averaged(green, left, right)
                assert amplitude.real+2*amplitude.imag > arb(floor)
    assert arb(".0011")/arb(5).sqrt() > arb(".00049")
    assert 2*(arb(".004")-arb(".0011")) > arb(".004")+arb(".0011")
    assert averaged(green, ORBIT, ORBIT).real > arb(".0419")
    # 在原八站加4、1。两端在 ORBIT 的各站均有相同的严格正权。
    values = [arb(".999")*averaged(green, (6,), v)+arb(".001")*averaged(green, ORBIT, v)
              for v in ((1,), (5,), ORBIT)]
    assert all((a.conjugate()*b).imag < arb("-1e-6") for a,b in zip(values, values[1:]+values[:1]))


def check_equal_energy(laplacian, green, energy_bounds, amplitude_bounds):
    for i in SITES:
        for j in SITES:
            special = {i,j} in ({6,3}, {6,5})
            assert green[i,j].real > (arb("-.00085") if special else arb(".006"))
            if i != j:
                assert abs(laplacian[i,j].real) <= 1
    assert laplacian[6,6].real == 18
    assert all(laplacian[i,i].real <= (10 if i in (3,5) else 18) for i in SITES)
    c = arb(110)/137
    assert c-(1-c)**2/(2*c) > arb(7)/9
    assert 19*(arb(7)/9)**2-1 > 10

    size = len(SITES)
    q = arb_mat([[laplacian[i,j].real for j in SITES] for i in SITES])
    assert all(q[i,i] > sum(abs(q[i,j]) for j in range(size) if i != j) for i in range(size))
    vector = q.solve(arb_mat([[1] for _ in SITES]))
    assert all(v > 0 for v in vector.entries())
    energy = 1/sum(vector.entries())
    weights = [energy*v for v in vector.entries()]
    amplitude = sum(weights[i]*weights[j]*green[k,l]
                    for i,k in enumerate(SITES) for j,l in enumerate(SITES))
    assert arb(energy_bounds[0]) < energy < arb(energy_bounds[1])
    assert arb(amplitude_bounds[0]) < amplitude.real < arb(amplitude_bounds[1])
    print("唯一极小能量与静态振幅", energy.str(13), amplitude.real.str(13))


def verify(replay):
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 192
    check_support(shell)
    rows = (
        (120, ("1.1358227698", "1.1358227700"), (".02406504479", ".02406504480")),
        (163, ("1.1357065862", "1.1357065863"), (".02409255315", ".02409255316")),
        (177, ("1.1357585435", "1.1357585436"), (".02403516866", ".02403516867")),
    )
    # 第四个可行编号独立用普适的支撑压缩判据核对等权约束。
    square = run_path(str(replay/"ch39-proton-next-79-heights.py"))
    laplacian = shell["laplacian"](5, acb.zeta_zero(205).imag)
    green = laplacian.inv()
    check_equal_energy(laplacian, green, ("1.1354354902", "1.1354354904"),
                       (".02418744109", ".02418744111"))
    for tied in ((0, 3), ORBIT):
        vertices = [tied]+[(i,) for i in SITES if i not in tied]
        for i in vertices:
            assert square["separated"]([averaged(green, i, j) for j in vertices])
        for (i, k), (j, m) in combinations_with_replacement(combinations(vertices, 2), 2):
            square["exclude_square"](*(averaged(green, x, y)
                                       for x, y in ((i, j), (i, m), (k, j), (k, m))))

    for index, energies, amplitudes in rows:
        laplacian = shell["laplacian"](5, acb.zeta_zero(index).imag)
        green = laplacian.inv()
        check_equal_orbit(green)
        check_equal_energy(laplacian, green, energies, amplitudes)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)

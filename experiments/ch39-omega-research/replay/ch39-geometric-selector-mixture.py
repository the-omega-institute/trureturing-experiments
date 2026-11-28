"""十九站磁子图多项式的正混合可实现指定的十二次谱选择器。"""
from functools import cache
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, arb_poly, ctx


# 位 0 至 6 依次按轨道最小站号排列。每对掩码共享同一边界约定。
SUPPORT = (
    (3, 116, "I"), (4, 4, "D"), (29, 60, "I"), (46, 108, "I"),
    (60, 35, "I"), (60, 87, "I"), (60, 124, "I"), (61, 124, "I"),
    (63, 127, "D"), (123, 127, "D"), (126, 116, "I"), (127, 35, "I"),
)


def squared_factors(energies):
    # x=t/9-1。所有十二次多项式在此坐标中仍取首一。
    return prod((arb_poly([1-e/9, 1])**2 for e in energies), start=arb_poly([1]))


def minimum(matrix):
    values = matrix.eig(multiple=True)
    assert len(values) == matrix.nrows()
    assert all(z.is_finite() and z.imag.contains(0) for z in values)
    return min(z.real.lower() for z in values).union(min(z.real.upper() for z in values))


def geometry(shell, pump):
    geo = shell["geo"]
    # 从整数坐标重新构造精确边相位的高精度包围。
    geo["edges"] = [(i, j, d, geo["unit"](i, j)*(-1 if d == 200 else 1))
                    for i, j, d, _ in geo["edges"]]
    sites = geo["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    orbits = [sorted({i, rotation[i], rotation[rotation[i]]}) for i in range(len(sites))
              if i == min(i, rotation[i], rotation[rotation[i]])]
    assert list(map(len, orbits)) == [3]*6+[1]
    print("原壳轨道（站号从零起）", orbits)
    frames = list(pump["rotation_bases"](shell))
    adjoint = pump["adjoint"]
    matrices = {}
    for charge in (-4, 2):
        base, links = pump["clock_matrices"](shell, charge)
        off = links[0]+links[1]
        matrices[charge] = base+off+adjoint(off)

    @cache
    def triple(charge, mask, mode):
        keep = sorted(i for k, orbit in enumerate(orbits) if mask & (1 << k) for i in orbit)
        matrix = matrices[charge]
        block = acb_mat([[matrix[i,j] for j in keep] for i in keep])
        if mode == "I":
            degrees = dict.fromkeys(keep, 0)
            for i, j, _, _ in geo["edges"]:
                if i in degrees and j in degrees:
                    degrees[i] += 1
                    degrees[j] += 1
            for row, i in enumerate(keep):
                block[row,row] = degrees[i]
        else:
            assert mode == "D"
        roots = []
        for frame in frames:
            columns = [j for j in range(frame.ncols()) if any(abs(frame[i,j]) > 0 for i in keep)]
            basis = acb_mat([[frame[i,j] for j in columns] for i in keep])
            roots.append(minimum(adjoint(basis)*block*basis))
        return roots

    return triple


def simplex_inverse(recipes):
    # 首项相同，零次项只平移作用。余下 11 个系数与概率归一化构成方阵。
    matrix = arb_mat([[p[k] for p in recipes] for k in range(1,12)]+[[1]*12])
    assert abs(matrix.det()) > 0
    return matrix.inv()


def barycentric(inverse, target):
    weights = inverse*arb_mat([[target[k]] for k in range(1,12)]+[[1]])
    smallest = min(weights[j,0].lower() for j in range(12))
    assert smallest > 0
    # 系数向量在无穷范数下移动时，到每个重心坐标为零的超平面的距离。
    radius = min((weights[j,0]/sum(abs(inverse[j,k]) for k in range(11))).lower()
                 for j in range(12))
    return weights, smallest, radius


def rationalize(weights, recipes, target):
    # 前十一项向下取整，末项补齐，使分母相同的有理权重严格归一。
    denominator = 10**110
    numbers = []
    for j in range(11):
        numerator = (denominator*weights[j,0]).floor().unique_fmpz()
        assert numerator is not None
        numbers.append(int(numerator))
    numbers.append(denominator-sum(numbers))
    assert min(numbers) > 0
    mixture = sum(((arb(n)/denominator)*p for n,p in zip(numbers, recipes, strict=True)),
                  start=arb_poly([0]))
    # 原变量作用差最多比较六个能量。磁拉普拉斯的加权度 <=18，
    # 故 t 在 [0,36]，|x|<=3。首项按有理权重和精确相消。
    error = 6*9**12*sum(3**k*abs(mixture[k]-target[k]) for k in range(1,12))
    return numbers, error.upper()


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    selector = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    ctx.prec = 640
    triple = geometry(shell, pump)
    recipes = [squared_factors(triple(-4, up, mode)+triple(2, down, mode))
               for up, down, mode in SUPPORT]
    target = squared_factors([arb(r) for r in selector["ROOTS"]])
    weights, smallest, radius = barycentric(simplex_inverse(recipes), target)
    assert smallest > arb("3.80e-9") and radius > arb("6.59e-8")
    print("最小概率与系数单形内半径", smallest.str(10), radius.str(10))
    numbers, error = rationalize(weights, recipes, target)
    assert error < arb("3.63e-91") < arb("5.47e-76")
    print("有理权重分子（公分母 10^110）", numbers)
    print("全编号作用差扰动上界", error.upper().str(10))


if __name__ == "__main__":
    verify()

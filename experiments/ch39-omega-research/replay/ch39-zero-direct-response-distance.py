"""从原始十九点响应的复区间矩阵验证 11431 的最近编号与误差半径。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
memory = run_path(str(root/"ch39-three-clock-minimal-memory.py"))
readout = run_path(str(root/"ch39-zero-inverse-response-distance.py"))
shell = readout["shell"]


def positive_pivots(matrix):
    # 厄米消元只用下三角；实对角取区间实部以利用自伴性。
    n = matrix.nrows()
    for k in range(n):
        pivot = matrix[k, k].real
        if pivot < 0:
            return False
        assert pivot > 0
        for i in range(k + 1, n):
            multiplier = matrix[i, k]/pivot
            for j in range(k + 1, i + 1):
                matrix[i, j] -= multiplier*matrix[j, k].conjugate()
    return True


def radius(gram, identity):
    # 半径使用 10^-12 作单位。上端为正定，下端已越过最大奇异值。
    low, high = arb(0), arb(3)
    assert positive_pivots(high**2*identity-gram)
    for _ in range(55):
        mid = (low+high)/2
        if positive_pivots(mid**2*identity-gram):
            high = mid
        else:
            low = mid
    return low.union(high)/10**12


def verify():
    ctx.prec = 320
    _, _, _, g0, g1, g2 = memory["coefficients"]()
    _, edges = shell["geometry"]()
    old_edges = [edge for edge in edges if edge[1] < 19]
    identity = acb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    phases = {n: 1-1/acb.zeta_zero(n) for n in (11430, 11431, 11432, 11433)}
    old = [shell["matrices"](range(19), old_edges, phases[n]) for n in (11431, 11432)]
    bounds = []
    for z in (acb(1), acb(-1, arb(3).sqrt())/2):
        offset = z*identity - z**2*acb_mat(g0) - z*acb_mat(g1) - acb_mat(g2)
        left, right = [(offset-matrix).inv() for matrix in old]
        difference = (left-right)*(10**12)/2
        norm = radius(difference.conjugate().transpose()*difference, identity)
        bounds.append(norm)
        print("Direct response half distance", norm.str(17))
    assert bounds[0] > bounds[1]
    for n in (11430, 11433):
        inverse_distance = arb(2).sqrt()*27/52136*abs(phases[n]-phases[11431])
        assert arb(".902")**2*inverse_distance/2 > max(bounds)
    return bounds


if __name__ == "__main__":
    verify()

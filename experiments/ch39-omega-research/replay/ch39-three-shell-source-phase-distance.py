"""截止 52 的同号门与已校准负号门所产生的第五壳源算子距离。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb


def frobenius_upper(matrix):
    return sum(z.abs_upper()**2 for z in matrix.entries()).sqrt()


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    basis, columns, levels = source["verify"](root, return_operator="spectral")
    station = basis * columns
    order = run_path(str(root / "ch39-two-shell-ordered-bond-gates.py"))
    negative = run_path(str(Path(__file__).with_name("ch39-three-shell-final-sign.py")))
    angles = list(map(arb, negative["ANGLES"]))
    for j in negative["ACTIVE"]:
        angles[j] += arb(0, arb(".00008").upper())

    head, first_leaf = zip(*order["LAYERS"][0][1])
    head_b, early = zip(*order["LAYERS"][2][1])
    head_c, late = zip(*order["LAYERS"][3][1])
    assert head == head_b == head_c
    second_round = {i for edge in order["LAYERS"][1][1] for i in edge}
    assert second_round.isdisjoint(head + first_leaf + early + late)
    assert len(set(head + first_leaf + early + late)) == 12
    assert 14 not in set(head + first_leaf + early + late) | second_round
    a, b, _, _ = order["block"](angles[:4])[0]
    diagonal = order["block"](angles[8:])[0][3]
    # Rows are indexed by the two heavy eigenvectors and the three rays.
    matrix = acb_mat([[a*station[h,k] + b*station[p,k], station[u,k], station[v,k]]
                      for k in (0, 1)
                      for h,p,u,v in zip(head, first_leaf, early, late, strict=True)])
    mid = matrix.mid()
    adj = mid.conjugate().transpose()
    left = (adj * mid).inv() * adj
    eye = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    singular_floor = (1 - frobenius_upper(eye - left*matrix)) / frobenius_upper(left)
    low, heavy, top = levels
    gap = min((heavy-low).lower(), (top-heavy).lower())
    phase_cost = sum(((1+1/abs(station[14,k]))/gap)**2 for k in (0, 1)).sqrt()
    bound = singular_floor*abs(diagonal)/phase_cost
    assert singular_floor > arb(".0061") and gap > arb(".3259")
    assert phase_cost < arb("272.20") and bound > arb(".000012")
    print("两条重线的三列最小奇异值下界、站点相位误差因子及同号算子距离下界",
          singular_floor.str(12), phase_cost.str(12), bound.str(12))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())

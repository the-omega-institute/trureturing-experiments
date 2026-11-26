"""105圈锐阶上完整维数的正权族，共用同一相位和全部正零点高度的全局真空。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    fixed = run_path(str(root / "ch39-height-independent-flavor-vacuum.py"))
    local = run_path(str(root / "ch39-mincut-reference-ckm-action.py"))
    global_model = run_path(str(root / "ch39-reference-flavor-global-winding.py"))
    lib = run_path(str(root / "ch39-mincut-global-winding.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root / "ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    tri = lib["terms"](shell)
    global_model["conjugation_symmetry"](root)
    original = fixed["minimum_order_weights"](lib, tri, x, y)
    solved = (30,35,38,40,105)
    free = [q for q in range(1,106) if q not in solved]
    assert len(free) == 100
    coefficients = {q:[arb(1), *fixed["gradient_columns"](lib,tri,q,x,y)]
                    for q in range(1,106)}
    matrix = arb_mat([coefficients[q] for q in solved]).transpose()
    assert not matrix.det().contains(0)
    # 原正权中心小幅放开，原来空缺的绕数先赋予相同的严格正权。
    centers = {q:original.get(q, arb("1e-8")) for q in free}
    radius = arb("1e-9")
    weights = {q:v+arb(0,radius) for q,v in centers.items()}
    target = [arb(int(k==0))-sum(weights[q]*coefficients[q][k] for q in free)
              for k in range(5)]
    dependent = matrix.solve(arb_mat([[v] for v in target])).entries()
    weights.update(zip(solved,dependent,strict=True))
    assert all(v>arb("8.9e-9") for v in weights.values())
    cap = arb(23)/25*arb(2)**-20/14**2
    delta = arb(0).union(cap)
    modes = [(q,[w*c for c in lib["coeff"](tri, q, delta)]) for q,w in weights.items()]
    drift = global_model["global_cover"](lib, tri, weights, modes, x, y)
    assert drift < arb("-3")
    print("独立系数数目、同时变化半宽、五个代回系数", len(free), radius.str(9),
          *[w.str(12) for w in dependent])


if __name__ == "__main__":
    verify()

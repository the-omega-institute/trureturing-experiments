"""从六质量、相对奇异谱及一行混合概率恢复三代 CKM 的模方和 CP 面积。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat, ctx


def verify():
    ctx.prec = 256
    root = Path(__file__).resolve().parent
    standard = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))["standard_matrix"]
    spectrum = run_path(str(root / "ch39-first-zero-relative-quark-spectrum.py"))["singular_values"]
    up = [arb(x) for x in (".00000615", ".00311", ".8616")]
    down = [arb(x) for x in (".0000135", ".000268", ".01401")]
    v = standard(*(arb(x).sin() for x in (".22704", ".04275", ".00377")), arb("1.139"))
    target = acb_mat([[up[i]*v[i][j]/down[j] for j in range(3)] for i in range(3)])
    squared = [s*s for s in spectrum(target)]
    first_moment = sum(squared)
    inverse_moment = sum(1/s for s in squared)
    a, b = ([x*x for x in masses] for masses in (up, down))
    def mass_boxes(centers, radii):
        return [center+arb(0, arb(radius).upper()) for center, radius in zip(centers, radii)]

    up_box = mass_boxes(up, (".00000014", ".00005", ".0043"))
    down_box = mass_boxes(down, (".0000002", ".000003", ".00009"))
    products = [u*d for u in up_box for d in down_box]
    assert all(not p.overlaps(q) for k, p in enumerate(products) for q in products[:k])
    first = [abs(z)**2 for z in v[0]]
    # 将最后一行依列和消去，剩下第二行的三个线性方程。
    constant = sum((a[0]*p+a[2]*(1-p))/mass for p, mass in zip(first, b))
    inverse_constant = sum((p/a[0]+(1-p)/a[2])*mass for p, mass in zip(first, b))
    rhs = arb_mat([[1], [(first_moment-constant)/(a[1]-a[2])],
                   [(inverse_moment-inverse_constant)/(1/a[1]-1/a[2])]])
    coefficient = arb_mat([[1, 1, 1], [1/x for x in b], b])
    assert not coefficient.det().contains(0)
    second = coefficient.solve(rhs).entries()
    recovered = [first, second, [1-p-q for p, q in zip(first, second)]]
    for i in range(3):
        for j in range(3):
            entry = recovered[i][j]
            assert entry > 0 and entry < 1
            assert entry.overlaps(abs(v[i][j])**2) and entry.rad() < arb("1e-40")
    lengths_sq = [p*q for p, q in zip(first, second)]
    x, y, z = lengths_sq
    j_squared = x*y - (z-x-y)**2/4
    j_abs = j_squared.sqrt()
    assert arb(".0000319") < j_abs < arb(".0000321")
    print("J 的绝对值", j_abs.str(18))

    def perpendicular(weights):
        return [weights[(i-1) % 3]-weights[(i+1) % 3] for i in range(3)]

    def moment_direction(rows, columns):
        return arb_mat([[u*v for v in perpendicular(columns)] for u in perpendicular(rows)])

    m = moment_direction(a, b)
    n = moment_direction([1/x for x in a], [1/x for x in b])
    # 由两个秩一方向消去 (0,2) 元素，保持 |V_ub|。
    direction = m*n[0, 2] - n*m[0, 2]
    direction = direction / abs(direction[2, 1])
    for t_sign, near, far in ((1, ".00113667", ".00113670"), (-1, ".00128453", ".00128456")):
        for t, inside in ((near, True), (far, False)):
            candidate = arb_mat(recovered) + direction*(t_sign*arb(t))
            assert all(entry > 0 for entry in candidate.entries())
            x, y, z = [candidate[0,j]*candidate[1,j] for j in range(3)]
            area = x*y - (z-x-y)**2/4
            signed_area = area if inside else -area
            assert signed_area > 0
        print("保持 |V_ub| 的实正交解，参数符号及绝对值区间", t_sign, near, far)


if __name__ == "__main__":
    verify()

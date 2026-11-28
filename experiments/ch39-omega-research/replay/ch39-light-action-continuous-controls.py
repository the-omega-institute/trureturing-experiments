"""在指定零点，核验沿真实高度曲线的局部谱作用及三个独立壳参数的能量坐标。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


def verify(replay):
    shell = run_path(str(replay / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay / "ch39-native-light-projector-pump.py"))
    ctx.prec = 256
    adj = pump["adjoint"]
    blocks = list(pump["rotation_bases"](shell))
    soft = shell["soft"]
    radial = acb_mat(19, 19)
    for site in soft["star"]["orbit"]:
        pole = soft["pole"]
        radial[pole, pole] += 1
        radial[site, site] += 1
        radial[pole, site] -= 1
        radial[site, pole] -= 1

    levels = []
    slopes = []
    curvatures = []
    for charge, index in ((-4, 2345666), (2, 964890)):
        height = acb.zeta_zero(index).imag
        lap = shell["laplacian"](charge, height)
        phases = shell["geo"]["phases"](height)
        theta = shell["geo"]["theta"]
        _, edges = pump["clock_matrices"](shell, charge)
        directions = []
        first = acb_mat(19, 19)
        second = acb_mat(19, 19)
        for edge, length in zip(edges, (724, 804), strict=True):
            amplitude = edge * phases[length]**charge
            derivative = acb(0, 1) * (amplitude - adj(amplitude))
            speed = charge * (theta(length) - theta(244))
            directions.append(derivative)
            first += speed * derivative
            second -= speed**2 * (amplitude + adj(amplitude))
        denom = height**2 + soft["point"]["study"]["offset"]
        deficit = arb(23) / 25 * arb(2)**-20 / denom
        first += (2*height*deficit/denom) * radial
        second += (2*deficit/denom - 8*height**2*deficit/denom**2) * radial
        directions.append(-radial)

        jacobian = []
        slope = [arb(0)] * 4
        curvature = [arb(0)] * 4
        for block in blocks:
            H = adj(block) * lap * block
            vals, vec = H.eig(right=True, algorithm="rump")
            order = sorted(range(len(vals)), key=lambda j: float(vals[j].real.mid()))
            assert all(z.is_finite() and z.imag.contains(0) for z in vals)
            lam = [vals[j].real for j in order]
            assert all(a < b for a, b in zip(lam, lam[1:]))
            vectors = []
            for j in order:
                v = acb_mat([[vec[i, j]] for i in range(H.nrows())])
                vectors.append(block*v/sum(abs(x)**2 for x in v.entries()).sqrt())
            ground = vectors[0]
            jacobian.append([(adj(ground)*op*ground)[0, 0].real for op in directions])
            e = lam[0]
            levels.append(e)
            tangent = first * ground
            de = (adj(ground)*tangent)[0, 0].real
            dde = (adj(ground)*second*ground)[0, 0].real + sum(
                2 * abs((adj(v)*tangent)[0, 0])**2 / (e-other)
                for other, v in zip(lam[1:], vectors[1:], strict=True)
            )
            for k in range(1, 5):
                slope[k-1] += k*e**(k-1)*de
                curvature[k-1] += k*e**(k-1)*dde
                if k > 1:
                    curvature[k-1] += k*(k-1)*e**(k-2)*de**2
        slopes.append(slope)
        curvatures.append(curvature)
        det = arb_mat(jacobian).det()
        print("荷与三个独立控制的雅可比行列式", charge, det.str(17))
        assert (det > arb(".0080958") if charge == -4 else det < arb("-.0036311"))
    assert all(abs(a-b) > arb(".14958") for i, a in enumerate(levels) for b in levels[i+1:])

    degree_two = arb_mat([row[:2] for row in slopes])
    print("高度方向驻点的二次子式", degree_two.det().str(15))
    assert degree_two.det() > arb(".022798")
    for degree in (3, 4):
        # 首一且除一次、二次及最高次外，其余系数取零。
        free = degree_two.solve(arb_mat([[-row[degree-1]] for row in slopes]))
        coefficient = [*free.entries(), arb(degree == 3), arb(degree == 4)]
        second_at_targets = [
            sum(c*v for c, v in zip(coefficient, row, strict=True))
            for row in curvatures
        ]
        print("首一次数、系数与两高度二阶导数", degree,
              [c.str(15) for c in coefficient], [d.str(15) for d in second_at_targets])
        if degree == 3:
            assert second_at_targets[0] < -6 and second_at_targets[1] > 2
        else:
            assert second_at_targets[0] < -452 and second_at_targets[1] < -74


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)

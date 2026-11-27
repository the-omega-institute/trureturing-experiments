"""第七至第九壳的引力六式，最后一个 52 弦轮使用相反的门符号。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


ANGLES = ("-2.29931640625", "-0.3143824576973494", "-1.8145357643647222", "0.62744140625",
          "0.484375", "2.5890145395207624", "-2.049072265625", "2.390612743505245",
          "-2.133292208417787", "0.6018105380446326", "0.277587890625", "-0.345947265625")
ACTIVE = (1, 2, 5, 7, 8, 9)


def centered_residual(projected, levels, weights, heat, field):
    low, heavy, top = levels
    newton, vacuum = -arb(1), -arb("1e-122")
    for i, weight in enumerate(weights):
        amplitude = low + (heavy-low)*abs(projected[i,0])**2 + (top-low)*abs(projected[i,1])**2
        lower = max(low.lower(), amplitude.lower())
        upper = min(top.upper(), amplitude.upper())
        assert 0 < lower <= upper
        domain = lower.union(upper)
        center = domain.mid()
        squared = center**2
        k2, k4, slope2, slope4 = (arb(0) for _ in range(4))
        for multiplicity, alpha, mass in field:
            k2 += alpha * heat(2, mass/squared)
            k4 += multiplicity * heat(3, mass/squared)
            exposure = (-mass/domain**2).exp()
            slope2 += alpha * exposure
            slope4 += multiplicity * exposure
        # Sum the signed field derivatives before bounding the input uncertainty.
        shift = domain - center
        newton += weight*(squared*k2/2 + shift*domain*slope2)/arb.pi()
        vacuum -= weight*(squared**2*k4/4 + shift*domain**3*slope4/2)/arb.pi()
    return arb_mat([[newton], [vacuum]])


def verify(root):
    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    source = load("fifth-full-positive-joint")
    fifth, columns, levels = source["verify"](root, return_operator="spectral")
    order = load("two-shell-ordered-bond-gates")
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = load("shortest-common-complex-gate")["gravity_system"]
    station = fifth * columns
    last = order["LAYERS"][-1]
    assert [j for j, _ in order["LAYERS"]] == [0, 1, 2, 2]
    flip = [i for pair in last[1] for i in pair]
    assert len(set(flip)) == len(flip) == 6
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_096):
        _, _, _, basis, _, weights, heat, field = geometry(root, fifth, columns, levels, number=number)
        shells.append((basis.conjugate().transpose(), weights, heat, field))

    def system(point, selected, values):
        mats = order["transport"](station, point, parameters=ANGLES, active=ACTIVE)
        # Each disjoint final gate is -G. This is exactly the original circuit
        # followed by a sign on its six final endpoints, for values and derivatives.
        for matrix in mats:
            for i in flip:
                for j in range(matrix.ncols()):
                    matrix[i,j] = -matrix[i,j]
        residuals, derivatives = [], []
        for left, weights, heat, field in selected:
            moved, *slopes = [left*m for m in mats]
            _, jac = gravity(moved, slopes, levels, weights, heat, field, False)
            if values:
                residuals.extend(centered_residual(moved, levels, weights, heat, field).tolist())
            derivatives.extend(jac.tolist())
        return arb_mat(residuals) if values else None, arb_mat(derivatives)

    center = [arb(ANGLES[i]) for i in ACTIVE]
    fixed, jacobian = system(center, shells[:3], True)
    pre = jacobian.mid().inv()
    assert not pre.det().contains(0)
    identity = arb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    radius = arb("0.00008")
    box = [v + arb(0, radius.upper()) for v in center]
    _, enclosure = system(box, shells[:3], False)
    move = pre * fixed
    defect = identity - pre * enclosure
    ratio = max(abs(move[i,0]).upper()/radius + sum(z.abs_upper() for z in row)
                for i, row in enumerate(defect.tolist()))
    print("六式自映比", ratio.str(12), flush=True)
    assert ratio < arb(".74")
    offset, _ = system(center, shells[3:], True)
    _, slopes = system(box, shells[3:], False)
    held = [offset[i,0] + arb(0, (radius*sum(z.abs_upper() for z in row)).upper())
            for i, row in enumerate(slopes.tolist())]
    newton, vacuum = 1 + held[0], arb("1e-122") + held[1]
    assert arb("1.1134") < newton < arb("1.1152")
    assert arb("-.00955") < vacuum < arb("-.00937")
    print("第十壳牛顿读数、真空读数", newton.str(14), vacuum.str(14), flush=True)


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())

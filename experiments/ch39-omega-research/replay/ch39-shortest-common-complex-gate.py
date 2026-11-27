"""三条有向极小边施同一个复旋转，复核第七壳的两个引力矩。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def gravity_system(projection, derivatives, levels, weights, heat, field, compute_values):
    low, heavy, top = levels
    gains = (heavy-low, top-low)
    target = (arb(1), arb("1e-122"))
    values = [-t for t in target]
    jacobian = [[arb(0) for _ in derivatives] for _ in target]
    for i, weight in enumerate(weights):
        amplitude = low + sum(gain*abs(projection[i, j])**2 for j, gain in enumerate(gains))
        left = max(low.lower(), amplitude.lower())**2
        right = min(top.upper(), amplitude.upper())**2
        assert 0 < left <= right
        dk, dw = arb(0), arb(0)
        for g, alpha, mass in field:
            exposure = (-mass/left).exp().union((-mass/right).exp())
            dk += alpha * exposure
            dw += g * exposure
            if compute_values:
                m2 = (left*heat(2, mass/left)).union(right*heat(2, mass/right))
                m4 = (left**2*heat(3, mass/left)).union(right**2*heat(3, mass/right))
                values[0] += weight*alpha*m2/(2*arb.pi())
                values[1] -= weight*g*m4/(4*arb.pi())
        slopes = (amplitude*dk/arb.pi(), -amplitude**3*dw/(2*arb.pi()))
        for j, direction in enumerate(derivatives):
            da = 2 * sum(gain*(projection[i, h].conjugate()*direction[i, h]).real
                         for h, gain in enumerate(gains))
            for row, slope in zip(jacobian, slopes, strict=True):
                row[j] += weight*slope*da
    return arb_mat([[value] for value in values]), arb_mat(jacobian)


def certify(project, center, radius, levels, weights, heat, field):
    box = [value + arb(0, radius) for value in center]
    residual, mid_jac = gravity_system(*project(center), levels, weights, heat, field, True)
    _, jac = gravity_system(*project(box), levels, weights, heat, field, False)
    preconditioner = mid_jac.mid().inv()
    defect = arb_mat([[1, 0], [0, 1]]) - preconditioner*jac
    offset = max(abs(value).upper() for value in (preconditioner*residual).entries())
    norm = max(sum(abs(value).upper() for value in row) for row in defect.tolist())
    ratio = offset/radius + norm
    print("有理坐标盒中心、半径", *(v.str(16) for v in center), radius.str(8))
    print("自映比上界、导数行列式", ratio.upper().str(16), jac.det().str(16), flush=True)
    assert ratio < 1 and not jac.det().contains(0)
    return box, ratio, jac.det()


def experiment(root, fifth, columns, levels):
    real_gates = run_path(str(root / "ch39-fifth-seventh-shortest-site-gates.py"))
    base, cosine, _, basis, station, weights, heat, field = real_gates["geometry"](
        root, fifth, columns, levels
    )
    pairs = ((6, 11), (7, 12), (8, 18))
    diagonal = sum(cosine[1:], cosine[0])
    fixed = base - diagonal

    def transfer(pairs):
        return acb_mat([[sum(basis[b, i].conjugate() * station[a, j] for a, b in pairs)
                         for j in range(2)] for i in range(19)])

    forward = transfer(pairs)
    reverse = transfer([(b, a) for a, b in pairs])

    def project(point):
        # 同一块 [[c,-conj(z)],[z,c]]，有理球面坐标满足 c²+|z|²=1。
        x, y = point
        denominator = 1 + x*x + y*y
        c = 1 - 2/denominator
        z = 2 * acb(x, y)/denominator
        projection = fixed + c*diagonal + z*forward - z.conjugate()*reverse
        derivatives = []
        for p, unit in zip(point, (acb(1), acb(0, 1)), strict=True):
            dc = 4*p / denominator**2
            dz = (2*unit - 2*p*z) / denominator
            derivatives.append(dc*diagonal + dz*forward - dz.conjugate()*reverse)

        return projection, derivatives

    center = list(map(arb, ("-.01781613411898", ".00257827694626")))
    box, ratio, determinant = certify(project, center, arb(".00002"), levels, weights, heat, field)
    assert ratio < arb(".206") and determinant < -139
    return box


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    return experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

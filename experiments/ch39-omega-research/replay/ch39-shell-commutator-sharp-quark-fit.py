"""在固定实际零点的十项最佳角点上确定壳实协方差加 i[H,P] 的全部解。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb, arb_mat


def diagonal(values):
    return acb_mat([[x if i == j else 0 for j in range(3)] for i, x in enumerate(values)])


def optimal_values(model, traces, affine):
    signs = (-1, 1, 1, -1, 1, 1, 1, 1, 1)

    def corner(radius):
        return [v + s * radius * w for v, s, w in
                zip(model["ref"][:9], signs, model["tol"][:9], strict=True)]

    def residual(radius):
        values = corner(radius)
        base, slope = affine(values)
        inverse_base, inverse_slope = affine(values, inverse=True)
        return inverse_base + inverse_slope * ((traces[0]-base)/slope) - traces[1]

    lo, hi = arb(".0574339450674740"), arb(".0574339450674742")
    assert residual(lo) < 0 < residual(hi)
    for _ in range(24):
        mid = (lo+hi)/2
        f = residual(mid)
        if f < 0:
            lo = mid
        elif f > 0:
            hi = mid
        else:
            break
    radius = lo.union(hi)
    values = corner(radius)
    base, slope = affine(values)
    delta = ((traces[0]-base)/slope).acos()
    assert 0 < delta and abs(delta-model["ref"][9]) < lo*model["tol"][9]
    return values + [delta], radius


def canonical_frame(model, matrix):
    levels, frame = model["eig"](model["adj"](matrix)*matrix)
    for j in range(3):
        pivot = max(range(3), key=lambda i: float(abs(frame[i,j]).mid()))
        phase = frame[pivot,j]
        assert all(abs(phase) > abs(frame[i,j]) for i in range(3) if i != pivot)
        for i in range(3):
            frame[i,j] *= phase.conjugate()/abs(phase)
    return levels, frame


def real_pencil_minors(model, covariance, shell_frame, generators):
    p, h0, h2 = generators
    ratio = h0[0,1].real/h0[1,2].real
    shell_coordinates = model["adj"](shell_frame)*covariance*shell_frame
    real_part = acb_mat([[shell_coordinates[i,j].real for j in range(3)] for i in range(3)])
    levels, basis = model["eig"](real_part)
    assert 0 < levels[0] < levels[1] < levels[2]
    minors = []
    for i in (1, 2):
        assert abs(basis[0,i]) > 0
        y, z = [(basis[j,i]/basis[0,i]).real for j in (1, 2)]
        h = y*(1-ratio*z)
        c = z*(ratio+z)-y*y
        t = ratio+z-ratio*y*y
        minor = (ratio*c-t)*(-t)-h*h
        assert abs(minor) > 0
        minors.append(minor)
    commutators = [a*real_part-real_part*a for a in (p, h0, h2)]
    determinant = arb_mat([[c[i,j].real for c in commutators]
                           for i,j in ((0,1), (0,2), (1,2))]).det()
    assert determinant < 0
    return minors, determinant


def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))
    orbit = load("two-zero-orbit-condition")
    model, _, _ = orbit["calibration"]()
    adj = model["adj"]
    at_zero = load("quadratic-gram-source")["at_zero"]
    up, down = at_zero(model, -4, 2345666), at_zero(model, 2, 964890)
    # 固定下部门各行的本征向量相位，消去底层求谱对相位原点的选择。
    row_factors = []
    for i in range(3):
        entry = down[i,0]
        assert abs(entry) > 0
        row_factors.append(entry.conjugate()/abs(entry))
    down = diagonal(row_factors)*down
    traces = load("two-zero-sharp-quark-fit")["normalized_traces"](up, down)
    values, radius = optimal_values(model, traces, load("nine-window-cp-cosine")["affine_trace"])
    mixing = acb_mat(model["load"]("shell-clock-ckm-blocks")["standard_matrix"](
        *(v.sin() for v in values[6:9]), values[9]))
    relative = up*down.inv()
    down_mass = diagonal(values[3:6])
    target = diagonal(values[:3])*mixing*down_mass.inv()
    scale = (abs(target.det())/abs(relative.det()))**(arb(1)/3)
    j_levels, j_frame = canonical_frame(model, relative)
    t_levels, t_frame = canonical_frame(model, target)
    assert all((j*scale**2).overlaps(t) for j, t in zip(j_levels, t_levels, strict=True))
    k0 = down.inv()*j_frame*adj(t_frame)*down_mass
    assert abs(k0.det()) > 0
    parts = []
    for j in range(3):
        column = acb_mat([[j_frame[i,j]] for i in range(3)])
        row = adj(acb_mat([[t_frame[i,j]] for i in range(3)]))
        parts.append((down.inv()*column)*(row*down_mass))
    reduced_gram = adj(t_frame)*down_mass*down_mass*t_frame
    assert all(abs(reduced_gram[i,j]) > 0 for i,j in ((0,1),(0,2),(1,2)))
    shell = load("shell-commutator-next-zeros")
    axes, commutator = shell["shell_axes"](model)
    shell_frame, h = load("shell-pencil-six-mass-window")["shell_frame"](model)
    pole = model["shell"]["soft"]["pole"]
    p = diagonal([sum(abs(model["sources"][pole,j])**2 for j in range(3)), 0, 0])
    words = load("two-zero-quadratic-source")
    equations = words["orbit_coefficients"](parts, adj, axes)
    trig = words["trigonometric_columns"]
    centers = (("2.334939524658", "4.995259501083"), ("2.726739122764", "2.824565186904"),
               ("5.217003118501", "2.578471812730"), ("5.690059534560", "5.261178983661"))
    certified, leaves = shell["certify_centers"](equations, trig, centers)
    print("最优半宽倍数、正 CP 相位、上部门标尺", radius.str(17), values[9].str(17), scale.str(17))
    print("两相位环面覆盖终格数", leaves)
    intervals = (("-.046", "-.045"), (".025", ".026"), ("-.0052", "-.0048"), (".034", ".035"))
    for (x0, y0), (lower, upper) in zip(certified, intervals, strict=True):
        half = arb("1e-9")
        pad = arb(0, half.upper())
        x, y = x0+pad, y0+pad
        vector, jacobian = trig(x0, y0)
        pre = (equations*jacobian).mid().inv()
        step = pre*equations*vector
        error = arb_mat([[1, 0], [0, 1]]) - pre*equations*trig(x,y)[1]
        assert all(abs(step[i,0])/half+sum(abs(error[i,j]) for j in range(2)) < 1
                   for i in range(2))
        k = orbit["phase_sum"](parts, orbit["phases_at"](x,y))
        assert abs(k.det()) > 0
        covariance = k*adj(k)
        coefficient = (commutator*covariance).trace().real/(commutator*commutator).trace().real
        assert arb(lower) < coefficient < arb(upper)
        pair = [model["eig"](adj(yukawa)*yukawa) for yukawa in (scale*up*k, down*k)]
        observations = [v.sqrt() for levels, _ in pair for v in levels]
        observations.extend(model["mixing_angles"](pair[0][1], pair[1][1]))
        assert all(v.overlaps(w) for v,w in zip(observations, values, strict=True))
        print("两相位与交换子系数", x.str(14), y.str(14), coefficient.str(13))
        minors, determinant = real_pencil_minors(model, covariance, shell_frame, (p, h, h*h))
        print("中、大谱线所属实铅笔的二阶子式", [minor.str(8) for minor in minors])
        print("排除单实谱方向的交换子行列式左右界",
              determinant.lower().str(8), determinant.upper().str(8))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)

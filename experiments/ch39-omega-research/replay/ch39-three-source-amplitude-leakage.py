"""共同厄米汤川振幅的两种平方次序及非零泄漏。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def real_coordinates(matrix, common):
    diagonal_and_cross = common["coordinates"](matrix)
    return arb_mat([[v if i < 3 else arb(2).sqrt()*v] for i, v in enumerate(diagonal_and_cross)])


def target_amplitudes(ckm, common):
    y_u = list(map(arb, (".00000615", ".00311", ".8616")))
    y_d = list(map(arb, (".0000135", ".000268", ".01401")))

    def diagonal(values):
        return acb_mat([[v if i == j else 0 for j in range(3)] for i, v in enumerate(values)])

    adj = common["adj"]
    v = acb_mat(ckm["standard_matrix"](*ckm["centers"]))
    a, b = adj(v)*diagonal(y_u)*v, diagonal(y_d)
    # 共同味基的第二轴选最重下型本征矢，第三轴选上型算子的正交残差。
    residual_norm = (abs(a[0, 2])**2 + abs(a[1, 2])**2).sqrt()
    assert residual_norm > 0
    z = [acb(0, 1)*a[j, 2]/residual_norm for j in (0, 1)] + [acb(0)]
    x = [-z[1].conjugate(), z[0].conjugate(), acb(0)]
    unitary = common["columns"](x, [acb(0), acb(0), acb(1)], z)
    assert (adj(unitary)*unitary).contains(diagonal([1, 1, 1]))
    return [adj(unitary)*matrix*unitary for matrix in (a, b)], y_u, y_d, v


def moment_frame(matrix, frame, sites, common):
    adj = common["adj"]
    levels, vectors = matrix.eig(right=True, algorithm="rump")
    assert all(v.imag.contains(0) for v in levels)
    order = sorted(range(19), key=lambda i: levels[i].real.mid())
    spectrum = [levels[j].real for j in order]
    projections = []
    data = []
    for j in order:
        v = adj(frame)*acb_mat([[vectors[i,j]] for i in sites])
        norm_sq = sum(abs(vectors[i,j])**2 for i in range(19))
        assert norm_sq > 0
        weight = v*adj(v)*(1/norm_sq)
        h = weight.trace().real
        assert h > 0
        projections.append(weight)
        data.append((h, real_coordinates(weight, common)))
    assert sum(projections, acb_mat(3, 3)).contains(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
    metric = sum((d*d.transpose()*(1/h) for h, d in data), arb_mat(6, 6))
    return spectrum, projections, data, metric


def verify(root):
    ray = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ckm = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    common = run_path(str(root / "ch39-signed-three-source-quark-pair.py"))
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    ctx.prec = 320
    frame = common["adapted_frame"](native["gram_cholesky"](root))
    sites = [ray["soft"]["pole"], *ray["boundary_sites"]()]
    height = acb.zeta_zero(2_306_888_100_085).imag
    matrices, up, down, ckm_matrix = target_amplitudes(ckm, common)
    all_levels, left_projections = [], []
    window = [
        ((".00000601", ".00000629"), (".00306", ".00316"), (".8573", ".8659")),
        ((".0000133", ".0000137"), (".000265", ".000271"), (".01392", ".01410")),
    ]
    for q, target, masses, box in zip((-4,2), matrices, (up, down), window, strict=True):
        spectrum, projections, data, metric = moment_frame(ray["laplacian"](q,height), frame, sites, common)
        all_levels.extend(spectrum)
        eig = acb_mat(metric).eig(algorithm="rump")
        assert all(z.imag.contains(0) for z in eig)
        eig = sorted((z.real for z in eig), key=lambda x: x.mid())
        assert all(a < b for a, b in zip(eig,eig[1:]))
        assert 0 < eig[0] and eig[-2] < 1 and eig[-1].contains(1)
        factor = 1/eig[-2]-1
        spread = sum((arb(box[j][0])-arb(box[i][1]))**2 for j in range(3) for i in range(j))/3
        floor = factor*spread
        max_trace = sum(arb(high)**2 for low,high in box)
        width = max_trace - sum(arb(low)**2 for low,high in box)
        assert floor > width
        print(q, "最大无迹本征值、泄漏因子、整窗绝对与相对下界、平方质量迹窗宽",
              *(v.str(12) for v in (eig[-2],factor,floor,floor/max_trace,width)))

        target_vector = real_coordinates(target, common)
        dual = metric.solve(target_vector, algorithm="precond")
        weights = [(d.transpose()*dual)[0,0]/h for h,d in data]
        block = sum((w*p for w,p in zip(weights,projections)), acb_mat(3,3))
        assert (block-target).contains(acb_mat(3,3))
        full_square = sum((w**2*p for w,p in zip(weights,projections)), acb_mat(3,3))
        leakage = full_square-target*target
        trace = leakage.trace().real
        cost = (target_vector.transpose()*(dual-target_vector))[0,0]
        assert trace.overlaps(cost)
        leak_eig = sorted((v.real for v in leakage.eig(algorithm="rump")),key=lambda v:v.mid())
        pos = sum(w>0 for w in weights)
        neg = sum(w<0 for w in weights)
        assert (pos,neg) == (14,5)
        bounds = (".00888523", ".00888525") if q == 2 else ("10.06537", "10.06539")
        assert arb(bounds[0]) < trace < arb(bounds[1])
        print(q,"最小泄漏中心见证正负数、迹与泄漏谱",pos,neg,trace.str(12),*[v.str(12) for v in leak_eig])
        squared_levels, projectors = flavor["projections"](target*target)
        assert all(v > 0 and v.overlaps(y**2) for v,y in zip(squared_levels,masses))
        left_projections.append(projectors)
    union = sorted(all_levels,key=lambda v:v.mid())
    assert all(a < b for a,b in zip(union,union[1:]))
    pu, pd = left_projections
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(ckm_matrix[i,j])**2)
               for i in range(3) for j in range(3))
    j = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    assert j > 0 and j.overlaps((ckm_matrix[0,0]*ckm_matrix[1,1]*ckm_matrix[0,1].conjugate()*ckm_matrix[1,0].conjugate()).imag)
    print("有符号 J",j.str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

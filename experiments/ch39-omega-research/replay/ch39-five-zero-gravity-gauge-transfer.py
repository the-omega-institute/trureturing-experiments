"""在五个实际零点背景检验共同的引力矩与零动量一圈规范阈值。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


BASE = 2_306_888_100_085
CANDIDATES = (2_345_666, 964_890, 27_804_036, 57_113_204)
ORIGINAL_SITES = ("6.5e-13", ".03", ".5", "1")
SHARED_SITES = (
    "4.217e-13", "3.1623e-12", "5.6234e-12", "1e-11", "1.7783e-11",
    "3.1623e-11", "5.6234e-11", "1e-10", "2.3714e-10", "4.217e-10",
    "1.3335e-09", "3.1623e-09", "7.4989e-09", "1.7783e-08", "3.1623e-08",
    "5.6234e-08", "1e-07", "1.7783e-07", "3.1623e-07", "7.4989e-07",
    "2.3714e-06", "3.1623e-05", ".03", "0.04217", ".5", "1",
)


def dot(a, b):
    return sum(x*y for x, y in zip(a, b, strict=True))


def apply_weights(columns, weights):
    return [dot(row, weights) for row in zip(*columns, strict=True)]


def certified_weights(matrix, target):
    """预条件剩余量给出精确解每一分量的绝对误差上界。"""
    n = matrix.nrows()
    midpoint = arb_mat([[matrix[i, j].mid() for j in range(n)] for i in range(n)])
    with ctx.workprec(640):
        approx_inverse = midpoint.inv()
    inverse = arb_mat([[approx_inverse[i, j].mid() for j in range(n)] for i in range(n)])
    identity = arb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    error = identity - inverse*matrix
    norm = max(sum(abs(error[i, j]) for j in range(n)).upper() for i in range(n))
    assert norm < arb("1e-14"), norm
    approx = arb_mat([[x.mid()] for x in (inverse*target).entries()])
    residual = inverse*(target-matrix*approx)
    radius = max(abs(v).upper() for v in residual.entries())/(1-norm)
    assert radius < arb(".02")
    return [w+arb(0, radius.upper()) for w in approx.entries()], norm, radius


def verify():
    root = Path(__file__).resolve().parent
    local = run_path(str(root/"ch39-signed-gravity-unification.py"))
    spectra = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))
    indices = (BASE, *CANDIDATES, BASE+1)
    single, lib, rows = spectra["gravity_fields"](root, indices, tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 320
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    ratio = (beta[0]-beta[2])/(beta[0]-beta[1])
    d = (1-ratio, ratio, -1)
    charge_weights = local["gauge_weights"](model)
    epsilon = arb("1e-122")
    def columns(row, sites):
        return [local["threshold_column"](row, single["moments"], lib["expint"],
                                          charge_weights, arb(x)) for x in sites]

    # 在已发表编号上复原四节点精确解。
    original = columns(rows[0], ORIGINAL_SITES)
    matrix = arb_mat([c[:3]+[dot(d, c[3:])] for c in original]).transpose()
    assert not matrix.det().contains(0)
    weights = matrix.solve(arb_mat([[1], [1], [epsilon], [-dot(d, initial)]])).entries()
    shifts = apply_weights(original, weights)[3:]
    intercept = [a+b for a, b in zip(initial, shifts, strict=True)]
    logscale = (intercept[0]-intercept[1])/(beta[0]-beta[1])
    scale = model["mz"]*logscale.exp()
    common = intercept[0]-beta[0]*logscale
    assert arb("5.03052e16") < scale < arb("5.03053e16") and arb("37.8345") < common < arb("37.8346")
    assert min(intercept) > 0
    print("共同目标能标 (GeV)、统一逆耦合", scale.str(13), common.str(13))

    floors = (".0015051", ".0278667", ".0374642", ".0204877")
    for number, row, floor in zip(CANDIDATES, rows[1:-1], floors, strict=True):
        _, k, w, *shift = apply_weights(columns(row, ORIGINAL_SITES), weights)
        delta = [u-v for u, v in zip(shift, shifts, strict=True)]
        residual = dot(d, delta)
        best = abs(residual)/sum(abs(a) for a in d)
        assert k > 1 and w > epsilon and residual > 0
        assert arb(floor) < best < arb(floor)+arb("1e-7")
        print(number, "K−1, W−ε, d·A, 最小最大逆耦合误差",
              *[v.str(12) for v in (k-1, w-epsilon, residual, best)])

    # 归一化、基准背景五行，其余四个背景各用五条差行。
    joint = [columns(row, SHARED_SITES) for row in rows[:-1]]
    basis = []
    for base, *others in zip(*joint, strict=True):
        differences = [v-u for other in others for u, v in zip(base[1:], other[1:], strict=True)]
        basis.append(base+differences)
    matrix = arb_mat(basis).transpose()
    target = arb_mat([[v] for v in [1, 1, epsilon, *shifts, *([0]*20)]])
    assert matrix.nrows() == matrix.ncols() == 26
    fitted, norm, radius = certified_weights(matrix, target)
    assert all(w > 0 or w < 0 for w in fitted)
    cost = sum(-w for w in fitted if w < 0)
    assert arb("2.2806415e13") < cost < arb("2.2806416e13")
    print("26 节点解的预条件误差范数、分量误差上界、负权总量",
          norm.str(8), radius.str(8), cost.str(14))
    print("按节点递增的共同权", *[v.str(12) for v in fitted])
    # 对下一枚的差列计算留出残差，以保留共同大数项的消去。
    held = columns(rows[-1], SHARED_SITES)
    differences = [[v-u for u, v in zip(c[1:], e[1:], strict=True)]
                   for c, e in zip(joint[0], held, strict=True)]
    dk, dw, *dshift = apply_weights(differences, fitted)
    projection = dot(d, dshift)
    print(BASE+1, "第六背景 K−1, W−ε, d·A", *[v.str(12) for v in (dk, dw, projection)])
    assert arb("-.003361") < dk < arb("-.003359")
    assert arb("-2e-9") < dw < 0 and 179874 < projection < 179876


if __name__ == "__main__":
    verify()

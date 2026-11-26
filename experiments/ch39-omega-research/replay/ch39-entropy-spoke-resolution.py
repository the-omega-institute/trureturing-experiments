"""核验一条中性辐边对高度的反演以及固定二进读数的尾部。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx

root = Path(__file__).parent
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
flavor = run_path(str(root / "ch39-cp-ensemble-test.py"))
scale = 2_903_600_000


def radial_coefficients(offset, strength):
    # 抬升常数零模，在固定的十九维空间展开正谱预解式。
    size = len(shell["geo"]["sites"])
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    laplacian = arb_mat(size, size)
    spokes = arb_mat(size, size)
    for i, j, _, _ in shell["geo"]["edges"]:
        matrices = [laplacian]
        if pole in (i,j) and (i in tuned or j in tuned):
            matrices.append(spokes)
        for matrix in matrices:
            matrix[i,i] += 1
            matrix[j,j] += 1
            matrix[i,j] -= 1
            matrix[j,i] -= 1
    lift = arb_mat([[arb(1)/size]*size for _ in range(size)])
    identity = arb_mat([[int(i == j) for j in range(size)] for i in range(size)])
    inverse = [(laplacian+lift).inv()]
    matrix_coefficients = []
    eta = [arb(1)]
    for k in range(5):
        eta.append(inverse[k].trace()-(-1)**k)
        shift = strength/scale**2*(-offset/scale**2)**k
        matrix_coefficients.append(identity-spokes*shift if k == 0 else -spokes*shift)
        rest = arb_mat(size,size)
        for j in range(k+1):
            rest += matrix_coefficients[j]*inverse[k-j]
        inverse.append(-inverse[0]*rest)
    assert eta[1] > 0
    joint = run_path(str(root / "ch39-five-neutral-quark-probability.py"))
    ctx.prec, ctx.cap = 512, 3
    _, _, poles, slopes = joint["model"]["shared"]["flavor_data"]()
    gradient = [arb(0), *[sum(s*(-z)**k for s,z in zip(slopes,poles,strict=True)) for k in range(6)]]
    mixed = [sum(eta[j]*gradient[k-j] for j in range(min(k,5)+1)) for k in range(7)]
    polynomial = arb_poly([0,1])
    for j in range(1,6):
        polynomial *= arb_poly([j*j,1])
    assert all(polynomial(-j*j).contains(0) for j in range(1,6))
    transverse = sum(polynomial[k]*mixed[k] for k in range(1,7))
    assert arb("-.0009790") < transverse < arb("-.0009789")
    print("中性响应在零端的一阶系数与六阶线性组合", eta[1].str(14), transverse.str(14))


def verify():
    ctx.prec = 384
    tuned = shell["soft"]["star"]["orbit"]
    pole = shell["soft"]["pole"]
    assert len(tuned) == 3 and pole not in tuned
    offset = shell["soft"]["point"]["study"]["offset"]
    strength = arb(23)/25*arb(2)**-20
    radial_coefficients(offset, strength)
    first = acb.zeta_zero(8359).imag
    for h in (first, arb(15*10**12), arb(24027719528722), arb(32*10**12)):
        shift = strength/(h*h+offset)
        laplacian = shell["laplacian"](0, h)
        for vertex in tuned:
            actual = 1+laplacian[pole,vertex].real
            assert actual.overlaps(shift) and actual > 0
            denominator = strength-offset*actual
            assert denominator > 0
            x = scale**2*actual/denominator
            assert x.overlaps((scale/h)**2)
    def binary_bits(h):
        delta = strength/(h*h+offset)
        return int((-delta.log()/arb(2).log()).ceil().unique_fmpz())
    assert [binary_bits(h) for h in (first, arb(15*10**12), arb(24027719528722))] == [47,108,110]
    # 六面味坐标在高处参考邻域确有变化，U 对 c 的导数在参考点非零。
    model = flavor["flavor"]
    reference = model["singlet"](11661)[1]
    axes = [model["singlet_at"]({724:v,804:acb(1)})[1] for v in (acb(1),acb(0,1))]
    assert sum(v*v for v in axes) > reference**2
    up, _ = flavor["two"]["probabilities"](reference, flavor["two"]["spectral_data"]())
    assert up[1] < 0
    joint = run_path(str(root / "ch39-five-neutral-quark-probability.py"))
    rows, *_ = joint["columns_at"](scale, (8359,12336))
    difference = sum(abs(a-b) for a,b in zip(rows[0][0][1:],rows[1][0][1:],strict=True))
    multiplier = -arb(".006779").log()/difference
    assert multiplier > arb("308.52")
    print("按七个原读数单位的最大绝对乘子下界", multiplier.str(13))
    print("首枚及两处高度界，直接向下读取辐边亏额所需的二进小数位", 47,108,110)


if __name__ == "__main__":
    verify()

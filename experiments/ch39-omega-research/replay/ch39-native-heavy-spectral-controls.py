"""在原生十九点壳上，用一个背景的谱迹读出两钟角及径向亏额。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


NUMBER = 2_306_888_100_085
N = 19


def soft_edge(shell, i, j):
    return (shell["soft"]["pole"] in (i, j)
            and (i in shell["soft"]["star"]["orbit"] or j in shell["soft"]["star"]["orbit"]))


def first_two(shell):
    # 精确整数系数。从无权度 r_i 及逐点软边数 v_i 展开 sum_i(r_i-δv_i)^2+2 sum_e w_e^2。
    edges = shell["geo"]["edges"]
    degree = [0] * N
    affected = [0] * N
    for i, j, distance, phase in edges:
        for k in (i, j):
            degree[k] += 1
            affected[k] += int(soft_edge(shell, i, j))
    count = sum(affected) // 2
    assert (len(edges), count) == (135, 3)
    pairs = [(i, j) for i, j, _, _ in edges]
    assert len(set(pairs)) == len(pairs) and all(0 <= i < j < N for i, j in pairs)
    c0 = sum(r*r for r in degree) + 2*len(edges)
    c1 = -2*sum(r*v for r, v in zip(degree, affected)) - 4*count
    c2 = sum(v*v for v in affected) + 2*count
    assert (sum(degree), -sum(affected), c0, c1, c2) == (270, -6, 4290, -180, 18)
    print("一阶、二阶迹作为亏额多项式的精确系数", (270, -6), (c0, c1, c2))


def system(shell, x, y, delta, charge):
    matrix = acb_mat(N, N)
    dx, dy, dd = (acb_mat(N, N) for _ in range(3))
    clocks = {724: acb(0, x).exp(), 804: acb(0, y).exp()}
    for i, j, distance, direction in shell["geo"]["edges"]:
        spoke = shell["soft"]["pole"] in (i, j)
        soft = soft_edge(shell, i, j)
        phase = acb(1) if spoke else (direction*clocks.get(distance, 1))**charge
        weight = 1-delta if soft else arb(1)
        matrix[i, i] += weight
        matrix[j, j] += weight
        matrix[i, j] = -weight*phase
        matrix[j, i] = matrix[i, j].conjugate()
        if soft:
            dd[i, i] -= 1
            dd[j, j] -= 1
            dd[i, j] = phase
            dd[j, i] = phase.conjugate()
        if distance in clocks:
            assert not spoke
            axis = dx if distance == 724 else dy
            axis[i, j] = -acb(0, charge)*weight*phase
            axis[j, i] = axis[i, j].conjugate()
    return matrix, (dx, dy, dd)


def diagonal_trace(matrix):
    z = sum(matrix[i, i] for i in range(N))
    assert z.imag.contains(0)
    return z.real


def trace_product(a, b):
    z = sum(a[i, j]*b[j, i] for i in range(N) for j in range(N))
    assert z.imag.contains(0)
    return z.real


def jacobians(shell, x, y, delta):
    matrix, axes = system(shell, x, y, delta, 1)
    other, other_axes = system(shell, x, y, delta, 3)
    squared = matrix*matrix
    cubed = squared*matrix
    other_squared = other*other
    linear = [0, 0, -6]
    cubic = [3*trace_product(squared, axis) for axis in axes]
    quartic = [4*trace_product(cubed, axis) for axis in axes]
    other_cubic = [3*trace_product(other_squared, axis) for axis in other_axes]
    return arb_mat([linear, cubic, quartic]), arb_mat([linear, cubic, other_cubic])


def readings(shell, point, mixed):
    x, y, deficit = point
    native, _ = system(shell, x, y, deficit, 1)
    squared = native*native
    cubic = squared*native
    if mixed:
        other, _ = system(shell, x, y, deficit, 3)
        final_trace = diagonal_trace(other**3)
    else:
        final_trace = diagonal_trace(squared*squared)
    return arb_mat([[diagonal_trace(native)], [diagonal_trace(cubic)], [final_trace]])


def recover(shell, observation, inverse, bounds, mixed):
    # 每步取中点压住舍入半径。对最后所得盒内点另用原始迹残差给出严格误差界。
    point = [((lo+hi)/2).mid() for lo, hi in bounds]
    for _ in range(12):
        step = inverse*(observation-readings(shell, point, mixed))
        point = [(v+step[i, 0]).mid().max(lo).min(hi) for i, (v, (lo, hi)) in enumerate(zip(point, bounds))]
    assert all(lo < v < hi for v, (lo, hi) in zip(point, bounds))
    return point


def norm(matrix):
    return max(sum(abs(matrix[i, j]) for j in range(matrix.ncols())).upper()
               for i in range(matrix.nrows()))


def bound(center, boxes):
    identity = arb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    preconditioner = center.inv()
    contraction = max(norm(identity - preconditioner*jac) for jac in boxes)
    assert contraction < 1
    inverse_bound = norm(preconditioner)/(1-contraction)
    determinants = [arb_mat([[jac[i, j] for j in (0, 1)] for i in (1, 2)]).det() for jac in boxes]
    slopes = [jac[1, 1] for jac in boxes]
    print("全盒角雅可比行列式绝对值、荷一三阶迹的第二角斜率的下界",
          min(abs(z).lower() for z in determinants).str(10), min(z.lower() for z in slopes).str(10))
    print("预乘矩阵范数、整盒压缩率、绝对误差放大上界",
          *[a.str(14) for a in (norm(preconditioner), contraction, inverse_bound)])
    return contraction, inverse_bound


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    angles = run_path(str(root/"ch39-mincut-reference-ckm-action.py"))
    flavor = run_path(str(root/"ch39-chiral-neutral-flavor-phase.py"))
    ctx.prec = 256
    first_two(shell)
    height = acb.zeta_zero(NUMBER).imag
    clocks = shell["geo"]["phases"](height)
    x, y = (clocks[k].arg() for k in (724, 804))
    offset = shell["soft"]["point"]["study"]["offset"]
    deficit = arb(23)/25*arb(2)**-20/(height**2+offset)
    cap = arb(23)/25*arb(2)**-20/14**2
    radius = arb("1e-4")
    assert offset > 0 and acb.zeta_zero(1).imag > 14 and 0 < deficit < arb("1e-12") < cap
    reference_gap = angles["reference_angles"](flavor)[0]
    reference_x = arb("-.01046")
    reference_y = reference_x-reference_gap
    assert abs(x-reference_x) + arb("9.2e-6") < radius
    assert abs(y-reference_y) + arb("9e-6") < radius
    assert arb("1e-9") < radius
    print("零点编号、钟角、实际亏额、整盒亏额上界", NUMBER,
          *[v.str(19) for v in (x, y, deficit, cap)])
    for charge in (1, 3):
        native = shell["laplacian"](charge, height)
        test, _ = system(shell, x, y, deficit, charge)
        assert all(z.contains(0) for z in (native-test).entries())
        moments = [diagonal_trace(native**j) for j in (1, 2, 3, 4)]
        assert (moments[0]-(270-6*deficit)).contains(0)
        assert (moments[1]-(4290-180*deficit+18*deficit**2)).contains(0)
        print("荷与实际零点上的四个原始谱迹", charge, *[z.str(38) for z in moments])
    center = jacobians(shell, x, y, deficit)
    entire = jacobians(shell, x+arb(0, radius), y+arb(0, radius), arb(0).union(cap))
    print("同荷 +1，一、三、四阶")
    contraction, inverse_bound = bound(center[0], [entire[0]])
    assert contraction < arb(".111") and inverse_bound < 33
    assert entire[0][1, 1] > 24
    print("双荷 +1、+3，一、三、三阶")
    # 双荷的两个三阶导数较接近平行。细分仍共用同一预乘矩阵，整条连接线适用同一范数界。
    count = 8
    subradius = radius/count
    pieces = [(2*i+1-count)*subradius+arb(0, subradius) for i in range(count)]
    mixed = [jacobians(shell, x+a, y+b, arb(0).union(cap))[1] for a in pieces for b in pieces]
    contraction, inverse_bound = bound(center[1], mixed)
    assert contraction < arb(".342") and inverse_bound < 3274
    # 同一组预乘矩阵与参数盒，移开参考零点后合成读数，复核两个逆均能仅凭各自三项数据还原。
    bounds = [(x-radius, x+radius), (y-radius, y+radius), (arb(0), cap)]
    hidden = [x+radius*5/8, y-radius/3, 2*cap/3]
    for mixed, jac, error_bound in zip((False, True), center, (33, 3274)):
        observed = readings(shell, hidden, mixed)
        point = recover(shell, observed, jac.inv(), bounds, mixed)
        residual = norm(observed-readings(shell, point, mixed))
        assert max(abs(a-b) for a, b in zip(point, hidden)) < arb("1e-30")
        assert error_bound*residual < arb("1e-30")
        print("盒内另取背景的反演，最大原始迹残差、严格坐标误差界",
              residual.str(10), (error_bound*residual).str(10))
    print("同荷三阶及以下至多读出一项角标量。两荷二阶及以下均与两钟角无关。")


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent
    verify(root)

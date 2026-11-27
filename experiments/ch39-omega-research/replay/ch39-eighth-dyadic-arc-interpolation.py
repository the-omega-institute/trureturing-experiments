"""以三次平分生成八等分的九个同步点，检验半角坐标的双二次插值。"""
from flint import acb, arb, arb_mat, ctx


def eighth_step(endpoint):
    # 两个终点都在闭第一象限，三次取主平方根即沿短弧平分。
    for _ in range(3):
        endpoint = endpoint.sqrt()
    return endpoint


def half_tangent(unit):
    return unit.imag / (1 + unit.real)


def verify():
    ctx.prec = 384
    sphere_end = acb(287, 24 * arb(2).sqrt()) / 289
    state_end = acb(0, 1)
    sphere_step, state_step = map(eighth_step, (sphere_end, state_end))
    normalisation = half_tangent(sphere_end)
    rows = []
    for k in range(9):
        r = half_tangent(sphere_step**k) / normalisation
        t = half_tangent(state_step**k)
        rows.append([r**a * t**b for a in range(3) for b in range(3)])
    determinant = arb_mat(rows).det()
    assert arb('3.4547248203624623035e-26') < determinant < arb('3.4547248203624623037e-26')
    print('九点半角坐标的双二次取值行列式', determinant.str(22))


if __name__ == '__main__':
    verify()

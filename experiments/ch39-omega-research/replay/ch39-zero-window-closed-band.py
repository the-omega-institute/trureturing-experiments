"""从相邻三角矩与有序四元组枚举分别计算实际零点之和的闭碰撞带。"""
from itertools import combinations_with_replacement

from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    heights = [z.imag for z in acb.zeta_zeros(11422, 19)]
    pairs = [(x + y, 1 if i == j else 2)
             for (i, x), (j, y) in combinations_with_replacement(enumerate(heights), 2)]
    radius, step = arb(1) / 100, arb(1) / 100000
    outer_edge = radius + step
    mass = 0
    ramp, outer_ramp = arb(0), arb(0)
    outside = []
    for x, a in pairs:
        for y, b in pairs:
            distance = abs(x - y)
            if distance < radius:
                mass += a * b
            else:
                assert distance > outer_edge
                outside.append(distance)
            ramp += a * b * (radius - distance).max(0)
            outer_ramp += a * b * (outer_edge - distance).max(0)

    recovered = ((outer_ramp - ramp) / step).unique_fmpz()
    trivial = sum(weight**2 for _, weight in pairs)
    assert sum(weight for _, weight in pairs)**2 == 19**4
    assert recovered == mass == 775 and trivial == 703
    nearest = min(d.lower() for d in outside).union(min(d.upper() for d in outside))
    assert arb(".01095669423") < nearest < arb(".01095669425")
    print("总有序四元组、闭带数、不同无序对贡献、带外最近距离",
          19**4, mass, mass - trivial, nearest.str(14))


if __name__ == "__main__":
    verify()

"""核对书中正 Fejér 核在第 11431 枚附近九个实际零点上的二阶和。"""
from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    heights = [acb.zeta_zero(index).imag for index in range(11427, 11436)]
    # Keep ordered-pair multiplicities after combining the two permutations.
    pair_sums = [(x+y, 1 if i == j else 2)
                 for i, x in enumerate(heights) for j, y in enumerate(heights) if j >= i]
    assert sum(weight for _, weight in pair_sums) == 81
    width = arb(1)/4
    hard_count, triangular = 0, arb(0)
    for x, wx in pair_sums:
        for y, wy in pair_sums:
            difference = abs(x-y)
            if difference < width:
                hard_count += wx*wy
                triangular += wx*wy*(1-difference/width)
            else:
                assert difference > width
    assert hard_count == 345
    assert arb("245.55588760255") < triangular < arb("245.55588760257")
    print("有序四元组总数、带内数、Fejér 三角权重", 9**4, hard_count, triangular)


if __name__ == "__main__":
    verify()

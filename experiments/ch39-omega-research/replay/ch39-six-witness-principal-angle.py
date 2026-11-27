"""包围实现指定六壳十八式见证所需的最小最大主角。"""

from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, arb_mat


def positive_sqrt(matrix, determinant):
    """用 Cayley–Hamilton 恒等式计算二阶正定矩阵的正平方根。"""
    identity = arb_mat([[1, 0], [0, 1]])
    root = determinant.sqrt()
    return (matrix + root * identity) / (matrix[0, 0] + matrix[1, 1] + 2 * root).sqrt()


def certificate(source, image_weight, levels):
    floor, first, top = levels
    assert floor < first < top
    beta = image_weight / (top - floor)
    assert 0 < beta < 1
    a = sum(abs(source[i, 0]) ** 2 for i in range(7))
    c = sum(abs(source[i, 1]) ** 2 for i in range(7))
    off = abs(sum(source[i, 0].conjugate() * source[i, 1] for i in range(7)))
    compression = arb_mat([[a, off], [off, c]])
    identity = arb_mat([[1, 0], [0, 1]])
    # 用 Cauchy–Binet 求行列式，保留两条投影近乎平行时的小正量。
    determinant = sum(
        abs(source[i, 0] * source[j, 1] - source[i, 1] * source[j, 0]) ** 2
        for i in range(7) for j in range(i)
    )
    complement_det = 1 - a - c + determinant
    assert 0 < a < 1 and determinant > 0 and complement_det > 0

    angle = arb(".417425")
    x, y = angle.cos(), angle.sin()
    h = x * x * a + 2 * x * y * off + y * y * c
    k = y * y * beta
    assert 0 < k < h < 1
    lower = (h * (1 - k)).sqrt() - (k * (1 - h)).sqrt()

    source_in = positive_sqrt(compression, determinant)
    source_out = positive_sqrt(identity - compression, complement_det)
    image_in = arb_mat([[0, 0], [0, beta.sqrt()]])
    image_out = arb_mat([[1, 0], [0, (1 - beta).sqrt()]])
    turn = arb("-1.111250")
    rotation = arb_mat([[turn.cos(), -turn.sin()], [turn.sin(), turn.cos()]])
    overlap = source_in * image_in + source_out * rotation * image_out
    trace = sum(z * z for z in overlap.entries())
    product = overlap.det() ** 2
    small = 2 * product / (trace + (trace**2 - 4 * product).sqrt())
    upper = (1 - small).sqrt()

    degrees = 180 / arb.pi()
    angle_low, angle_high = lower.asin() * degrees, upper.asin() * degrees
    generator = arb(3).sqrt() * lower
    distance = (2 - 2 * (1 - lower**2).sqrt()).sqrt()
    assert lower > arb(".941781") and upper < arb(".941783")
    assert angle_low > arb("70.3529") and angle_high < arb("70.3532")
    assert generator > arb("1.63121") and distance > arb("1.15219")
    print("七维字符的源 Gram 对角、非对角模、行列式", *(v.str(14) for v in (a, c, off, determinant)))
    print("像的高重线投影范数平方", beta.str(14))
    print("最小投影交换子范数的下、上界", lower.lower().str(16), upper.upper().str(16))
    print("最小最大主角（度）的下、上界", angle_low.lower().str(16), angle_high.upper().str(16))
    print("三循环交换子、至保持七维字符的酉群距离的下界", generator.lower().str(16), distance.lower().str(16))
    return lower, upper


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    sites = lib("three-ray-chiral-boundary")["geo"]["sites"]
    forms, _ = lib("cyclic-two-characters-nine")["characters"](sites)
    restrict = forms[0].conjugate().transpose()
    source = restrict * (fifth * columns)
    weighted = lib("six-unrestricted-eighteen")["experiment"](root, fifth, columns, levels)
    # 六壳见证的首列是谱隙较大的重列。另一条重列由构造严格位于十二维正交补。
    top_column = acb_mat([[weighted[i, 0]] for i in range(weighted.nrows())])
    image_weight = sum(abs(z) ** 2 for z in (restrict * top_column).entries())
    return certificate(source, image_weight, levels)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))

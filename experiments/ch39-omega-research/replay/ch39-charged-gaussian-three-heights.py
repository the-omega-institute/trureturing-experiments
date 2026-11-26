"""三处实际高度上，任意磁边相位的同号高斯行列式选权障碍。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))


def verify():
    ctx.prec, ctx.cap = 320, 3
    h1, h2 = [acb.zeta_zero(n).imag for n in (8359, 12336)]
    ends = [arb(23672421776), arb(237077947240)]
    fractions = [(h2/h1).log()/(h/h1).log() for h in ends]
    angle = fractions[0].union(fractions[1])
    ratios = [(h/h2)**2 for h in ends]
    first_ratio = (h2/h1)**2
    assert arb(".019") < angle < arb(".024")
    assert arb("1.9") < first_ratio < 2 and arb(10)**12 < ratios[0] and ratios[1] < arb(10)**15
    assert (1000000/h2)**2 > 7000

    edges = shell["geo"]["edges"]
    degrees = [sum(k in (i, j) for i,j,_,_ in edges) for k in range(19)]
    assert len(edges) == 135 and max(degrees) == 18
    pole = shell["soft"]["pole"]
    assert degrees[pole] == 18
    assert shell["soft"]["point"]["study"]["offset"] > 1
    # 三种质量范围中，对行列式对数严格弦差的有理下界。
    trace_lower = arb(260)
    small = arb(".019")*trace_lower/(arb("1e-12")+arb(".001")+arb(36)/7000)-270
    middle = 19*(arb(".976")*arb(".9")/arb("2.1")-arb(".024")*arb(1001).log())-arb(270)/7000
    large = 19*(arb(".976")*arb(".9")/21-arb(".024"))-arb(270)/7000
    assert small > 0 and middle > 0 and large > 0
    # 在第三枚相对第二枚的质量零归一行列式上界，先消掉共同的尺度幂。
    gain = arb(10)**15*(arb("1e-12")+arb(36)/7000)**19
    assert gain < 1
    print("首三处高度的对数插值系数范围", angle.str(10))
    print("三段弦差对应的正裕量", small.str(12), middle.str(12), large.str(12))
    print("质量零的 D/x 在高点相对第二枚的上界", gain.str(12))


if __name__ == "__main__":
    verify()

"""相邻变化受控的双钟，可在周期四维格上有整体绕数。"""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path


def verify():
    root = Path(__file__).parent
    radius = 2047
    run_path(str(root/"ch39-plaquette-overlap-gap.py"))["tent"](radius)
    square = F(6,2*(radius+1)**2+1)
    slope = F(1, 2**22)
    # 一维归一帐篷移动一格的范数严格小于这个有理数。
    shift = F(7,4*(radius+1))
    assert square < shift**2
    gap = F(61,640)
    assert (1-gap+F(17,1000))**2 < 1-F(30,200)
    loss = 8*shift+4*radius*slope
    assert F(17,1000)-loss > F(8,1000)

    size = 2**26
    assert size >= 2*radius+3
    # 基本钟在一条边上最多一个向前绕一步，逐顶点可有独立复弦距误差。
    chord_error = F(1,2**32)
    step = 2*F(22,7)/size
    scale = F(24,5)
    root_two_upper = F(99,70)
    assert root_two_upper**2 > 2
    normalized_change = (root_two_upper*6*(step+4*chord_error)+6*chord_error)/scale
    assert normalized_change < slope
    print("帐篷半径、每格质量变化、谱隙损失有理上界、整体绕行的每格质量变化上界",
          radius, slope, loss, float(normalized_change))


if __name__ == "__main__":
    verify()

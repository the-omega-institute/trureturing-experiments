"""三条软辐边给出单调一阶迹，计算识别给定编号所需的绝对精度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    geo = shell["geo"]
    assert len(geo["edges"]) == 135 and len(shell["soft"]["star"]["orbit"]) == 3
    offset = shell["soft"]["point"]["study"]["offset"]
    assert offset > 0
    def missing(index):
        height = acb.zeta_zero(index).imag
        return 6*arb(23)/25*arb(2)**-20/(height**2+offset)

    def bits(bound):
        exponent = (-bound.log()/arb(2).log()).ceil().unique_fmpz()
        assert exponent is not None
        n = int(exponent)
        assert arb(2)**-n < bound < arb(2)**(1-n)
        return n

    for n, expected in ((2345666,(80,74,60)),(964890,(76,70,57))):
        center = missing(n)
        cutoffs = []
        for width in (1,65):
            left = (missing(n-width)-center)/2
            right = (center-missing(n+width))/2
            assert left > 0 and right > 0
            cutoffs.append(min(left.lower(),right.lower()))
        cutoffs.append(center/2)
        assert tuple(map(bits,cutoffs)) == expected
        print("实际编号、单号/±64/尾段阈值、绝对二进精度",
              n,[x.str(13) for x in cutoffs],expected)


if __name__ == "__main__":
    verify()

"""前两枚最小割三窗零点的二分读数与最短二进小数精度。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 256
    soft = shell["soft"]
    pole, outer = soft["pole"], set(soft["star"]["orbit"])
    edges = [(i, j) for i, j, _, _ in shell["geo"]["edges"]]
    star = [edge for edge in edges if pole in edge and outer.intersection(edge)]
    assert len(edges) == 135 and len(star) == len(outer) == 3
    assert all(sum(v in edge for edge in edges) == 10 for v in outer)
    subset = set(sorted(outer)[:2])
    assert sum(len(subset.intersection(edge)) == 1 for edge in star) == 2
    defects = []
    for index in (8236, 8359):
        h = acb.zeta_zero(index).imag
        delta = (arb(23)/25)*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])
        defects.append(delta*(2-delta))

    for count, first_bits in ((1, 48), (2, 47), (3, 48)):
        for bits in range(first_bits+1):
            # 对每个切割，整边数扣去 count*(2δ−δ²)。取整后只需比较扣除的整数。
            correction = [(count*d*2**bits).ceil().unique_fmpz() for d in defects]
            assert all(v is not None for v in correction)
            if bits < first_bits:
                assert correction[0] == correction[1]
            else:
                assert correction[0] != correction[1]
                print("切过软边数、首次分辨位数、两个整数修正", count, bits, *correction)


if __name__ == "__main__":
    verify()

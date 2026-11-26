"""认证第一、第二枚实际零点之间的未改权壳在半径 9 处的谱计数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 112
    heights = [acb.zeta_zero(index).imag for index in (1, 2)]
    intervals = 128
    for charge in (1, -4, 2, -3, 6):
        initial, final = [shell["laplacian"](charge, height) for height in heights]
        direction = final - initial
        bound = max(sum(abs(direction[i, j]) for j in range(19)).upper()
                    for i in range(19))
        margin = arb(100)
        for k in range(intervals):
            midpoint = arb(2*k + 1)/(2*intervals)
            values = sorted((v.real for v in (initial + midpoint*direction).eig(multiple=True)),
                            key=float)
            radius = bound/(2*intervals)
            assert all(v + radius < 9 for v in values[:3])
            assert all(v - radius > 9 for v in values[3:])
            local = min((abs(v - 9) - radius).lower() for v in values)
            margin = min(margin, local)
        assert margin > arb("0.54")
        print("荷、低于整数九的态数、离九的严格距离下界", charge, 3,
              margin.str(12), flush=True)


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

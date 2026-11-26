"""认证原壳可单独作为第四带的循环特征线在一块去点钟角域内延拓。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    sectors = run_path(str(root / "ch39-native-light-projector-pump.py"))
    ctx.prec = 112
    bases = list(sectors["rotation_bases"](shell))
    adj = sectors["adjoint"]
    deficit = 6*arb("1e-8")
    margin = arb("0.0001")

    @cache
    def phase(n, i):
        return acb(0, arb.pi()*(2*i+1)/n).exp()

    @cache
    def variation(n):
        return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin() + deficit

    for charge in (1, -4, 2, -3, 6):
        base, clocks = sectors["clock_matrices"](shell, charge)
        blocks = [[adj(q)*x*q for x in (base, *clocks)] for q in bases]
        misses = [None]*3
        leaves = 0
        finest = 8
        todo = [(8, i, j) for i in range(8) for j in range(8)]
        while todo:
            n, i, j = todo.pop()
            levels = []
            for d, link1, link2 in blocks:
                off = link1*phase(n, i) + link2*phase(n, j)
                ev = [v.real for v in (d+off+adj(off)).eig(multiple=True)]
                # 使用各区间端点的次序统计，允许同一特征标内有更高的重根。
                lo = sorted(v.lower() for v in ev)
                hi = sorted(v.upper() for v in ev)
                gap = min((lo[1]-hi[0]).lower(), (lo[2]-hi[1]).lower())
                levels.append((lo[1], hi[1], gap))
            fine = True
            for k, (low, _, gap) in enumerate(levels):
                other_high = min(levels[m][1] for m in range(3) if m != k)
                excess = low - other_high
                if misses[k] is None and excess > 2*deficit + arb("0.2"):
                    misses[k] = (n, i, j)
                if not (gap > 2*variation(n) + margin or excess > 2*variation(n) + margin):
                    fine = False
            if fine:
                leaves += 1
                finest = max(finest, n)
            else:
                assert n < 65536, (charge, n, i, j)
                todo.extend((2*n, 2*i+a, 2*j+b) for a in range(2) for b in range(2))
        assert all(m is not None for m in misses)
        print("荷、覆盖叶数、最细每圈分格数、三个特征标各自排除的角中点 (n,i,j)",
              charge, leaves, finest, *misses, flush=True)


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

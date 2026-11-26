"""对五荷的原壳前五态给出各循环特征标秩固定时的全亏额延拓。"""
from functools import cache
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


RANKS = ((3, 1, 1), (1, 3, 1), (1, 1, 3), (2, 2, 1), (2, 1, 2), (1, 2, 2))


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
    def twice_variation(n):
        return 2*(4*arb(2).sqrt()*(arb.pi()/(2*n)).sin() + deficit)

    for charge in (1, -4, 2, -3, 6):
        base, clocks = sectors["clock_matrices"](shell, charge)
        blocks = [[adj(q)*x*q for x in (base, *clocks)] for q in bases]
        misses = {}
        leaves = 0
        finest = 8
        todo = [(8, i, j) for i in range(8) for j in range(8)]
        while todo:
            n, i, j = todo.pop()
            levels = []
            for d, c1, c2 in blocks:
                off = c1*phase(n, i) + c2*phase(n, j)
                ev = [z.real for z in (d+off+adj(off)).eig(multiple=True)]
                levels.append((sorted(v.lower() for v in ev), sorted(v.upper() for v in ev)))
            fine = True
            for ranks in RANKS:
                cuts = [(lo[m-1], hi[m-1], lo[m], hi[m])
                        for m, (lo, hi) in zip(ranks, levels)]
                local_gap = min((u-l).lower() for _, l, u, _ in cuts)
                incompatibility = max(l for l, _, _, _ in cuts) - min(u for _, _, _, u in cuts)
                if ranks not in misses and incompatibility > 2*deficit + arb("0.1"):
                    misses[ranks] = (n, i, j)
                if not (local_gap > twice_variation(n) + margin
                        or incompatibility > twice_variation(n) + margin):
                    fine = False
            if fine:
                leaves += 1
                finest = max(finest, n)
            else:
                assert n < 65536, (charge, n, i, j)
                todo.extend((2*n, 2*i+a, 2*j+b) for a in range(2) for b in range(2))
        assert all(ranks in misses for ranks in RANKS)
        print("荷、覆盖叶数、最细每圈分格数、依序六种秩的排除角中点 (n,i,j)",
              charge, leaves, finest, *(misses[r] for r in RANKS), flush=True)


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)

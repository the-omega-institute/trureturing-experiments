"""认证五个荷在有效双钟全环面与对称软辐边区间上的三个最低本征值互异。"""
from functools import cache
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def certify(shell,pump,bases,charge,required):
    adj = pump["adjoint"]
    base,links = pump["clock_matrices"](shell,charge)
    blocks = [tuple(adj(frame)*matrix*frame for matrix in (base,*links)) for frame in bases]

    @cache
    def phase(n,i):
        return acb(0,arb.pi()*(2*i+1)/n).exp()

    @cache
    def variation(n):
        return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin()+6*arb("1e-8")

    cells = [(8,i,j) for i in range(8) for j in range(8)]
    leaves,finest = 0,8
    margin = arb(10)
    order = None
    while cells:
        n,i,j = cells.pop()
        z,w = phase(n,i),phase(n,j)
        lowest,second = [],[]
        for center,first_link,second_link in blocks:
            off = z*first_link + w*second_link
            roots = (center+off+adj(off)).eig(multiple=True)
            assert all(r.imag.contains(0) for r in roots)
            ordered = sorted((r.real for r in roots), key=lambda r:r.mid())
            assert all(ordered[0] < r for r in ordered[1:])
            lowest.append(ordered[0])
            second.append(min(r.lower() for r in ordered[1:]))
        cost = 2*variation(n)
        gap = min(abs(x-y).lower() for x,y in combinations(lowest,2))-cost
        heavy_gap = min(second)-max(x.upper() for x in lowest)-cost
        if gap > required and heavy_gap > 0:
            here = tuple(sorted(range(3),key=lambda k:lowest[k].mid()))
            if order is None:
                order = here
            assert order == here
            margin = min(margin,gap.lower())
            leaves += 1
            finest = max(finest,n)
        else:
            assert n < 4096, (charge,n,i,j)
            cells.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("荷、三字符按能量递增的编号、终格数、最细分格、全域间隔下界",
          charge,order,leaves,finest,margin.str(12),flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec = 160
    bases = list(pump["rotation_bases"](shell))
    for q,bound in ((1,".0000047671"),(-4,".0015440"),(2,".00077287"),(-3,".00033486"),(6,".0024642")):
        certify(shell,pump,bases,q,arb(bound))

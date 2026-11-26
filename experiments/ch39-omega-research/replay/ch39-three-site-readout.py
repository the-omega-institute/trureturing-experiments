"""比较三内点与两组三外点对整张双钟环面上轻模的读取。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


@lru_cache(maxsize=None)
def phase(n, i):
    # 区间 [(i/n),(i+1)/n] 圈的中点。
    return acb(0,arb.pi()*(2*i+1)/n).exp()


def lowest(matrix):
    values, columns = matrix.eig(right=True)
    indices = sorted(range(len(values)),key=lambda i: float(values[i].real))
    low = indices[0]
    assert all(values[low].real < values[j].real for j in indices[1:])
    second = min(values[j].real.lower() for j in indices[1:])
    v = [columns[j,low] for j in range(matrix.nrows())]
    norm = sum(abs(a)**2 for a in v).sqrt()
    return [a/norm for a in v], values[low].real, second-values[low].real


def inner_dark(block, clocks, inner):
    outer = [j for j in range(block.nrows()) if j != inner]
    principal = acb_mat([[block[i,j] for j in outer] for i in outer])
    v,energy,_ = lowest(principal)
    assert energy < arb("8.685")
    amplitudes = [abs(sum(a.conjugate()*m[i,inner] for a,i in zip(v,outer)))
                  for m in [block,*clocks]]
    assert amplitudes[0] > max(amplitudes[1:])
    return amplitudes[0]-amplitudes[1]-amplitudes[2]


def outer_readout(block, clocks, n, i, j):
    delta = 2*arb(2).sqrt()*arb.pi()/n + 6*arb("1e-8")
    edge = phase(n,i)*clocks[0]+phase(n,j)*clocks[1]
    matrix = block+edge+edge.conjugate().transpose()
    v,_,gap = lowest(matrix)
    if gap > delta:
        leakage = delta/(gap-delta)
        if min(abs(v[0]),abs(v[1]))-leakage > arb(1)/4:
            return 1
    assert n < 96
    return sum(outer_readout(block,clocks,2*n,2*i+di,2*j+dj)
               for di in range(2) for dj in range(2))


def verify():
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    common = run_path(str(root/"ch39-all-common-shell-ranks.py"))
    ctx.prec = 160
    bases = list(common["character_bases"](shell))
    sites = shell["geo"]["sites"]
    for q in (1,2,3,4,6):
        base,clock = common["build_family"](shell,q)
        stats = []
        for k,basis in enumerate(bases):
            orbits = [[i for i in range(19) if not basis[i,j].is_zero()]
                      for j in range(basis.ncols())]
            assert set(orbits[0]) == {0,3,4} and set(orbits[1]) == {1,2,5}
            assert sorted(sites[orbits[0][0]]) == sorted(sites[orbits[1][0]]) == [1,5,29]
            inner, = [j for j,o in enumerate(orbits) if set(o) == {6,7,8}]
            adjoint = basis.conjugate().transpose()
            block = adjoint*base*basis
            clocks = [adjoint*b*basis for b in clock]
            margin = inner_dark(block,clocks,inner)
            if k == 1:
                assert margin < arb("-.14")
            else:
                assert margin > arb(".14")
            cells = sum(outer_readout(block,clocks,8,i,j) for i in range(8) for j in range(8))
            stats.append((k,margin.str(5),cells))
        print("荷、字符内点三角余量与外点覆盖格数",q,stats)


if __name__ == "__main__":
    verify()

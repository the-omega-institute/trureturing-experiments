"""多枚零点的轻投影同时固定一个局域作用的最少背景数。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


local = run_path(str(Path(__file__).with_name("ch39-local-light-generator.py")))
mix, adjoint, sites = local["mix"], local["adjoint"], local["sites"]


def character_basis(k):
    rotated = [sites.index((y,z,x)) for x,y,z in sites]
    orbits = []
    for a in range(19):
        group = tuple(dict.fromkeys((a,rotated[a],rotated[rotated[a]])))
        if a == min(group) and (k == 0 or len(group) == 3):
            orbits.append(group)
    omega = acb(-1, arb(3).sqrt())/2
    return acb_mat([[omega**(k*o.index(s)) if s in o else 0 for o in orbits]
                    for s in range(19)])


def dark_columns(frames, k):
    basis = character_basis(k)
    n = min(len(frames), basis.ncols())
    vectors = acb_mat([[frame[i,k] for frame in frames[:n]] for i in range(19)])
    equations = adjoint(vectors)*basis
    def columns(indices):
        return acb_mat([[equations[i,j] for j in indices] for i in range(n)])
    pivots = max(combinations(range(basis.ncols()), n),
                 key=lambda ids: float(abs(columns(ids).det())))
    square = columns(pivots)
    assert abs(square.det()) > 0
    free = [j for j in range(basis.ncols()) if j not in pivots]
    if not free:
        return []
    solved = -square.solve(columns(free))
    coeff = acb_mat(basis.ncols(),len(free))
    for j,idx in enumerate(free):
        coeff[idx,j] = 1
    for i,idx in enumerate(pivots):
        for j in range(len(free)):
            coeff[idx,j] = solved[i,j]
    dark = basis*coeff
    return [[dark[i,j] for i in range(19)] for j in range(len(free))]


def verify():
    ctx.prec = 512
    numbers = (28325,48067,58214,22027,11661,11934,12207)
    frames = [mix["light_modes"](n) for n in numbers]
    for p,frame in enumerate(frames[1:]):
        overlap = adjoint(frames[0])*frame
        diagonal = [abs(overlap[k,k])**2 for k in range(3)]
        assert all(v > 0 for v in diagonal)
        if p == 0:
            assert all(not v.overlaps(w) for i,v in enumerate(diagonal) for w in diagonal[:i])

    # 四个秩三投影总可在十三点团内留下一个共同盲向。
    edges = local["forbidden"](804)
    assert not any(i>=6 and j>=6 for i,j in edges)
    assert len(sites[6:]) == 13
    dark = sum((dark_columns(frames[:5],k) for k in range(3)), [])
    assert len(dark) == 4
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    matrix = acb_mat([[int(rot[i]==j), int(rot[j]==i)]
                      + [a[i]*b[j].conjugate() for a in dark for b in dark]
                      for i,j in edges])
    assert matrix.nrows() == 72 and matrix.ncols() == 18
    local["positive"](adjoint(matrix)*matrix)
    print("五背景的公共正交补维数",len(dark),"；缺边约束核对的独立非标量方向",matrix.ncols())
    assert all(len(dark_columns(frames,k)) == 0 for k in range(3))
    print("七背景分别张成轮换特征空间的维数", *(character_basis(k).ncols() for k in range(3)))


if __name__ == "__main__":
    verify()

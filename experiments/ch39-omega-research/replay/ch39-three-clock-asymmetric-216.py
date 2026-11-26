"""截止 216 的三组固定列关系联合迫使外点平方范数为负。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

plane = run_path(str(Path(__file__).with_name("ch39-three-clock-plane-completion.py")))


def verify():
    points,distance,_,delta,_,_ = plane["shell_data"]()
    rotation = [points.index((p[1],p[2],p[0])) for p in points]
    rows = [0,2,4,6,8,10]
    cols = [1,3,7,12,15,16,17,18]
    relations = []
    for _ in range(3):
        assert all(distance[i][j]>=224 for i in rows for j in cols)
        block = plane["pairs"]["submatrix"](delta,rows,cols)
        selected = next(c for c in combinations(range(8),6)
                        if plane["pairs"]["submatrix"](block,range(6),c).det())
        pivot = plane["pairs"]["submatrix"](block,range(6),selected)
        for j in sorted(set(range(8))-set(selected)):
            coefficients = pivot.solve(plane["pairs"]["submatrix"](block,range(6),[j]))
            relation = [fmpq(0)]*19
            relation[cols[j]] = 1
            for k,l in enumerate(selected):
                relation[cols[l]] = -coefficients[k,0]
            relations.append(relation)
        rows,cols = [rotation[i] for i in rows],[rotation[j] for j in cols]
    free = [j for j in range(19) if distance[0][j]<224]
    assert free==[0,2,4,6]
    fixed = sorted(set(range(19))-set(free))
    system = fmpq_mat([[r[j] for j in free] for r in relations])
    rhs = fmpq_mat([[-sum(delta[0,j]*r[j] for j in fixed)] for r in relations])
    chosen = next(c for c in combinations(range(6),4) if plane["pairs"]["submatrix"](system,c,range(4)).det())
    values = plane["pairs"]["submatrix"](system,chosen,range(4)).solve(plane["pairs"]["submatrix"](rhs,chosen,[0]))
    assert system*values==rhs and values[0,0] < -fmpq(3733,10000)
    print("For cutoff below 224, rank six forces scaled outer squared norm",float(values[0,0]))


if __name__ == "__main__":
    verify()

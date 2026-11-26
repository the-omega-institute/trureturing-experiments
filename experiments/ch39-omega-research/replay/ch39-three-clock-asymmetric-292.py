"""在不要求补全对称时，用一对外点的内积排除截止 292 的秩五 Gram 补全。"""
from pathlib import Path
from runpy import run_path

from flint import fmpq

pairs = run_path(str(Path(__file__).with_name("ch39-three-clock-outer-pairs.py")))
previous = pairs["previous"]


def verify():
    _,_,_,g0,g1,g2 = previous["MEM"]["coefficients"]()
    delta = 10**10*(g2-g1*g0.solve(g1))
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17,17,17)})
    rows,cols = (1,3,5,7,9),(0,10,13,17,18)
    # 同一组已指定系数分别适用于本点对的两个方向。
    swap = [points.index((p[1],p[0],p[2])) for p in points]
    assert swap[0]==2 and swap[2]==0
    def coefficients(i,j,R,C):
        def far(r,c):
            return local["squared"](points[r],points[c])>=296
        assert C[0]==i and all(far(r,c) for r in R for c in (*C,j))
        assert all(far(i,c) for c in C[1:])
        block = pairs["submatrix"](delta,R,C)
        vector = block.solve(pairs["submatrix"](delta,R,[j]))
        return vector[0,0],sum(delta[i,c]*vector[k,0] for k,c in enumerate(C) if k)
    a,b = coefficients(0,2,rows,cols)
    assert (a,b)==coefficients(2,0,[swap[i] for i in rows],[swap[j] for j in cols])
    assert a>fmpq(11556,10000) and b>fmpq(579,10000)
    print("Rank five would require W02 = a W00 + b and W20 = a W22 + b, scaled by 1e10.")
    print("a =",float(a),"b =",float(b),"next permitted cutoff at least 296.")


if __name__ == "__main__":
    verify()

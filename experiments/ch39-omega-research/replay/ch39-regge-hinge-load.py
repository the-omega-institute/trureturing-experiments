"""在四维坐标链的平直胞上检查铰链负载的第一变分。"""
from collections import Counter, defaultdict
from itertools import permutations
from pathlib import Path
from runpy import run_path

import sympy as sp


def verify():
    flat = run_path(str(Path(__file__).with_name("ch39-regge-flat-cell.py")))
    simplex_edges = flat["EDGES"]
    q = sp.Matrix([sp.Rational((j-i)*(5-j+i),5) for i,j in simplex_edges])
    shape_change = sp.Matrix([int(j-i in (1,4)) for i,j in simplex_edges])
    local = []
    for face,cos,area_angle in flat["canonical"]():
        # 面积对平方边长一次齐次。先沿 q 收缩 dA，再除以该面面积。
        area = sp.sqrt(sp.Rational(1,5) if cos == 0 else sp.Rational(3,20))
        angle_change = (q.T*area_angle*shape_change)[0]*sp.sqrt(5)/area
        local.append((face,angle_change))

    per_hinge = defaultdict(list)
    for order in permutations(range(4)):
        vertices = [sum(1<<i for i in order[:j]) for j in range(5)]
        for (i,j,k),change in local:
            face = (vertices[j]-vertices[i],vertices[k]-vertices[j])
            per_hinge[face].append(change)
    coefficients = Counter((len(changes),sp.simplify(sum(changes))) for changes in per_hinge.values())
    assert coefficients == {(4,3):30,(6,-3*sp.sqrt(3)):20}
    assert sp.simplify(30*3/sp.sqrt(5)-20*3*sp.sqrt(3)*sp.sqrt(sp.Rational(3,20))) == 0
    # 有理局部负载沿 Pell 递推逼近两类几何面积比。
    a,b = 2,1
    for _ in range(100):
        assert a*a-3*b*b == 1
        a,b = 2*a+3*b,a+2*b
    print("每平移胞的铰链价数、角亏导数及对应面数",coefficients)


if __name__ == "__main__":
    verify()

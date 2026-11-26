"""从二十五个线性角亏读数同时恢复四维曲率与五个非度量分量。"""
from collections import defaultdict
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
from runpy import run_path

import sympy as sp


PLANES = list(combinations(range(4), 2))
PAIRS = list(combinations_with_replacement(range(6), 2))


def curvature_basis():
    # 第一 Bianchi 恒等式 R01,23 = R02,13 - R03,12。
    dependent = PAIRS.index((0,5))
    free = [i for i in range(21) if i != dependent]
    lift = sp.zeros(21,20)
    for j,i in enumerate(free):
        lift[i,j] = 1
    lift[dependent,free.index(PAIRS.index((1,4)))] = 1
    lift[dependent,free.index(PAIRS.index((2,3)))] = -1
    return lift, free


def readouts():
    flat = run_path(str(Path(__file__).with_name("ch39-regge-flat-cell.py")))
    edges = flat["EDGES"]
    q = sp.Matrix([sp.Rational((j-i)*(5-j+i),5) for i,j in edges])
    metric_indices = list(combinations_with_replacement(range(4),2))
    metric = sp.Matrix(15,10,lambda a,k: (1+int(metric_indices[k][0]!=metric_indices[k][1]))
                       *int(all((a+1)&(1<<d) for d in metric_indices[k])))
    heavy = sp.Matrix.hstack(*metric.T.nullspace())
    lift, _ = curvature_basis()
    steps = [[(mask>>i)&1 for i in range(4)] for mask in range(16)]
    rows = defaultdict(lambda: sp.zeros(1,25))
    sine = {}
    areas = {}
    for face, cosine, mixed in flat["canonical"]():
        area = sp.sqrt(sp.Rational(1,5) if cosine==0 else sp.Rational(3,20))
        sin_angle = sp.sqrt(1-cosine**2)
        dc = (q.T*mixed*sp.sqrt(5)*sin_angle/area).applyfunc(sp.simplify)
        for order in permutations(range(4)):
            v = [sum(1<<d for d in order[:j]) for j in range(5)]
            i,j,k = face
            key = v[j]-v[i], v[k]-v[j]
            sine[key] = sin_angle
            # 面定位于首顶点；平方边长采用 Riemann 正规坐标二阶面积度量。
            for col,(s,t) in enumerate(edges):
                if dc[col] == 0:
                    continue
                start = [h-g for h,g in zip(steps[v[s]],steps[v[i]])]
                edge = steps[v[t]-v[s]]
                biv = [edge[b]*start[d]-edge[d]*start[b] for b,d in PLANES]
                curvature = sp.Matrix([[-sp.Rational(1,3)*(1+int(b!=d))*biv[b]*biv[d]
                                        for b,d in PAIRS]])*lift
                rows[key] += dc[col]*curvature.row_join(heavy[v[t]-v[s]-1,:])
            sides = [edges.index(e) for e in combinations(face,2)]
            total = sum(q[col] for col in sides)
            da = sp.zeros(1,15)
            for col in sides:
                s,t = edges[col]
                da[v[t]-v[s]-1] = (total/2-q[col])/(8*area)
            if key in areas:
                assert areas[key]==da
            areas[key]=da
    assert len(rows)==50 and all(rows[a,b]==rows[b,a] for a,b in rows)
    keys = sorted((a,b) for a,b in rows if a<b)
    matrix = sp.Matrix.vstack(*(rows[key] for key in keys))
    assert matrix.det(method="domain-ge")==-2660511744
    return matrix, [sine[key] for key in keys], sp.Matrix.vstack(*(areas[key] for key in keys))


def operator_coefficients(matrix):
    # 每个角亏的曲率算子系数，使用二形式内积升指标。
    lift,_ = curvature_basis()
    measured = lift*matrix.inv()[:20,:]
    inv = sp.eye(4)+sp.ones(4)
    wedge = sp.Matrix(6,6,lambda i,j: inv[PLANES[i][0],PLANES[j][0]]*inv[PLANES[i][1],PLANES[j][1]]
                       - inv[PLANES[i][0],PLANES[j][1]]*inv[PLANES[i][1],PLANES[j][0]])
    return [wedge*sp.Matrix(6,6,lambda i,j: measured[PAIRS.index(tuple(sorted((i,j)))),k])
            for k in range(25)]


def verify():
    matrix, sine, _ = readouts()
    assert matrix.inv()*matrix==sp.eye(25)
    operators = operator_coefficients(matrix)
    bound = sp.zeros(6)
    for op in operators:
        bound += op.applyfunc(abs)
    print("二十五行行列式、四价和六价行数、迹三次多项式系数绝对值和上界",
          matrix.det(method="domain-ge"),sine.count(1),sine.count(sp.sqrt(3)/2),(bound**3).trace())


if __name__=="__main__":
    verify()

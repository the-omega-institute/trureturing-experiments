"""Exact zero-momentum Hessian of the 4D Regge area-deficit action on the rational metric."""
from collections import defaultdict
from itertools import combinations, permutations
import sympy as sp

VERTICES = range(5)
EDGES = list(combinations(VERTICES, 2))
R = sp.Rational


def canonical():
    coords = sp.Matrix(4, 4, lambda i,j: int(int(j) <= int(i)))
    metric = sp.eye(4)-sp.ones(4)/5
    inv = (coords*metric*coords.T).inv()
    grad = (-sp.ones(1,4)).col_join(sp.eye(4))
    normals = grad*inv*grad.T
    assert all(normals[i,i] == 2 for i in VERTICES)
    dq_normals = []
    for edge in EDGES:
        def dq(i, j):
            return int(tuple(sorted((int(i),int(j)))) == edge)
        dg = sp.Matrix(4, 4, lambda i,j: R(dq(0,i+1)+dq(0,j+1)-dq(i+1,j+1), 2))
        dq_normals.append(-grad*inv*dg*inv*grad.T)
    contributions = []
    for u, v in EDGES:
        triangle = tuple(i for i in VERTICES if i not in (u,v))
        face_edges = list(combinations(triangle, 2))
        def length(edge):
            i, j = edge
            return R((j-i)*(5-j+i),5)
        lengths = [length(edge) for edge in face_edges]
        cos = -normals[u,v]/2
        area2 = (2*sum(x*y for x,y in combinations(lengths,2))-sum(x*x for x in lengths))/16
        assert (cos, area2) in [(0, R(1,5)), (R(1,2), R(3,20))]
        divisor = 16*sp.sqrt(5*area2*(1-cos**2))
        a = sp.zeros(10, 1)
        for edge, value in zip(face_edges, lengths):
            a[EDGES.index(edge)] = (sum(lengths)-2*value)/divisor
        dc = sp.Matrix([-d[u,v]/2-cos*(d[u,u]+d[v,v])/4 for d in dq_normals])
        contributions.append((triangle, cos, a*dc.T))
    return contributions


def verify():
    parts = canonical()
    local = sum((m for _, _, m in parts), sp.zeros(10))
    H = sp.zeros(15)
    stars = defaultdict(list)
    for order in permutations(range(4)):
        vertices = [sum(1 << i for i in order[:k]) for k in VERTICES]
        indices = [vertices[j]-vertices[i]-1 for i, j in EDGES]
        for i, a in enumerate(indices):
            for j, b in enumerate(indices):
                H[a,b] += local[i,j]
        for (i,j,k), cos, _ in parts:
            shape = (vertices[j]-vertices[i], vertices[k]-vertices[j])
            stars[shape].append(cos)
    assert len(stars) == 50 and all(sorted(v) in ([0]*4, [R(1,2)]*6) for v in stars.values())
    metric_columns = [(i,j) for i in range(4) for j in range(i,4)]
    metric_modes = sp.Matrix(15, 10, lambda m, k: (1+int(metric_columns[k][0]!=metric_columns[k][1]))*
                             int(all((m+1) & (1<<v) for v in metric_columns[k])))
    assert H == H.T and metric_modes.rank() == 10 and H*metric_modes == sp.zeros(15,10)
    assert H.eigenvals() == {0: 10, R(-6,5): 4, R(-11,5): 1}
    print('Zero momentum squared-edge Hessian divided by sqrt(5):', H.eigenvals())
    print('Triangle orbits with four or six incident four-simplices:',
          sum(len(s)==4 for s in stars.values()), sum(len(s)==6 for s in stars.values()))


if __name__ == '__main__':
    verify()

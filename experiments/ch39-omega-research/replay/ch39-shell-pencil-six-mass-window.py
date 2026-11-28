"""实际夸克制备盒中，壳三循环与极点的全部实铅笔谱源的六质量迹比。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb


def shell_frame(model):
    adj, source = model["adj"], model["sources"]
    sites = model["shell"]["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    rotation = acb_mat([[int(i==j) for j in cycle] for i in range(len(sites))])
    h = acb(0,1)*adj(source)*(rotation-adj(rotation))*source
    pole = model["shell"]["soft"]["pole"]
    v = adj(acb_mat([[source[pole,j] for j in range(3)]]))
    basis = []
    # Krylov 内积是 H 的实矩，Gram–Schmidt 只需要实系数。
    for column in (v, h*v, h*h*v):
        for q in basis:
            column -= q*(adj(q)*column)[0,0].real
        length = sum(abs(z)**2 for z in column.entries()).sqrt()
        assert length > 0
        basis.append(column/length)
    frame = acb_mat([[q[i,0] for q in basis] for i in range(3)])
    matrix = adj(frame)*h*frame
    assert matrix[0,1].real > arb(".00058") and matrix[1,2].real > arb(".9")
    return frame, matrix


def eigenline(t, h, reciprocal):
    a, b = h[0,1].real, h[1,2].real
    d, e = h[1,1].real, h[2,2].real
    if reciprocal:
        return ((d*t-1)*(e*t-1)-(b*t)**2, -a*t*(e*t-1), a*b*t*t)
    return ((d-t)*(e-t)-b*b, -a*(e-t), a*b)


def form(matrix, vector):
    return sum(vector[i]*vector[j]*matrix[i,j].real for i in range(3) for j in range(3))


def verify(replay):
    orbit = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))
    model, up, down = orbit["calibration"]()
    adj = model["adj"]
    frame, h = shell_frame(model)
    factor = abs(down.det()/up.det())**(arb(2)/3)
    gu, gd = [adj(v*frame)*(v*frame) for v in (up, down)]
    bound = 34*gd-factor*gu
    # 纯 P 的零空间是二维复空间，直接用复厄米主子式。
    block = acb_mat([[bound[i,j] for j in (1,2)] for i in (1,2)])
    assert bound[0,0].real > 0 and block[0,0].real > 0 and block.det().real > 0

    half = arb(1)/2
    cells = [(arb(-2), arb(2), 0, False), (-half, half, 0, True)]
    leaves = deepest = 0
    while cells:
        left, right, depth, reciprocal = cells.pop()
        vector = eigenline(left.union(right), h, reciprocal)
        if form(bound, vector) > 0:
            leaves += 1
            deepest = max(deepest, depth)
        else:
            assert depth < 40
            mid = (left+right)/2
            cells.extend(((left,mid,depth+1,reciprocal), (mid,right,depth+1,reciprocal)))
    # 非零首坐标使这一向量成为某个实铅笔的谱线。
    vector = eigenline(arb(".9008063"), h, False)
    assert abs(vector[0]) > 0
    assert factor*form(gu,vector)/form(gd,vector) > arb("33.7")
    shape = run_path(str(replay/"ch39-two-zero-source-correlations.py"))["shape"]
    masses = [m+arb(0,t.upper()) for m,t in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    assert shape(masses[:3])/shape(masses[3:]) > arb("74.25")
    print("全部正铅笔谱源迹比的上确界在 (33.7,34)；六质量要求大于 74.25。",
          "终格与最大深度", leaves, deepest)


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)

"""穷尽固定十观测轨道上有向壳三循环读数的驻点，并检验实三次源法向。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, acb, acb_mat


def verify():
    root = Path(__file__).resolve().parent
    neighbours = run_path(str(root / "ch39-four-quadratic-source-neighbours.py"))
    load = neighbours["load"]
    orbit, words, gram = [load(k) for k in (
        "two-zero-orbit-condition", "two-zero-quadratic-source", "quadratic-gram-source")]
    model, _, _ = orbit["calibration"]()
    adj, source = model["adj"], model["sources"]
    sites = model["shell"]["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    rotation = acb_mat([[int(i==j) for j in cycle] for i in range(len(sites))])
    observable = acb(0,1)*adj(source)*(rotation-adj(rotation))*source
    up, down = [gram["at_zero"](model,q,n) for q,n in ((-4,2345666),(2,964890))]
    parts, _ = orbit["orbit"](model, up, down)
    u, d = [adj(m)*m for m in (up,down)]
    frame, _ = words["normal_frame"](u/u.trace().real, d/d.trace().real)
    coefficients = words["orbit_coefficients"](parts, adj, [observable, frame[0]])
    _, dcos, dsin, xcos, xsin, ycos, ysin = coefficients.tolist()[0]
    grad = arb_mat([[0,dsin,-dcos,xsin,-xcos,0,0],
                    [0,-dsin,dcos,0,0,ysin,-ycos]])
    columns = words["trigonometric_columns"]
    radius = arb(".0001")
    pad = arb(0, radius.upper())
    centers = [(arb(x), arb(y)) for x, y in (
        ("1.969471","3.630613"),("1.969470",".489020"),
        ("5.111063",".489020"),("5.111063","3.630613"))]
    intervals = (("3.0336","3.0338"),("3.3370","3.3373"),
                 ("3.2538","3.2542"),("3.1167","3.1171"))
    for k, ((x,y), (lo,hi)) in enumerate(zip(centers,intervals,strict=True)):
        values, deriv = columns(x,y)
        vbox, dbox = columns(x+pad,y+pad)
        h = grad*deriv
        inv = arb_mat([[h[i,j].mid() for j in range(2)] for i in range(2)]).inv()
        shift = inv*grad*values
        error = arb_mat([[1,0],[0,1]]) - inv*grad*dbox
        for i in range(2):
            assert abs(shift[i,0])/radius + sum(abs(error[i,j]) for j in range(2)) < arb(".1")
        h = grad*dbox
        if k < 2:
            assert h.det() > 0 and (1 if k==0 else -1)*h[0,0] > 0
        else:
            assert h.det() < 0
        value, normal = (coefficients*vbox).entries()
        assert arb(lo) < value < arb(hi) and abs(normal) > arb(".012")
        print("驻点读数、对实三次法向的配对",value.str(8),normal.str(8))
    cells = [(0,0,0)]
    leaves = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0, half.upper())
        x, y = (2*i+1)*half+pad, (2*j+1)*half+pad
        gradient = grad*columns(x,y)[0]
        if any(not z.contains(0) for z in gradient.entries()) or any(
                abs(x-a)<radius and abs(y-b)<radius for a,b in centers):
            leaves += 1
        else:
            assert depth < 25
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("两相位完整覆盖终格",leaves)

    # 谱函数 S=f(H) 必须与 H 对易。用厄米正交单位基逐项读 i[H,S]。
    axes = []
    for i in range(3):
        diagonal = acb_mat(3,3)
        diagonal[i,i] = 1
        axes.append(diagonal)
        for j in range(i):
            for phase in (acb(1), acb(0,1)):
                basis = acb_mat(3,3)
                basis[i,j] = phase/arb(2).sqrt()
                basis[j,i] = basis[i,j].conjugate()
                axes.append(basis)
    readings = [acb(0,1)*(axis*observable-observable*axis) for axis in axes]
    readings.append(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
    moments = words["orbit_coefficients"](parts, adj, readings)
    cells = [(0,0,0)]
    commutator_leaves = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0,half.upper())
        vals = moments*columns((2*i+1)*half+pad,(2*j+1)*half+pad)[0]
        lower = sum(vals[k,0].abs_lower()**2 for k in range(9))
        if lower > (arb(".64")*abs(vals[9,0]).upper())**2:
            commutator_leaves += 1
        else:
            assert depth < 17
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("交换子范数大于 0.64 倍协方差迹，终格数", commutator_leaves)



if __name__ == "__main__":
    verify()

"""穷尽固定十观测轨道上有向壳三循环读数的驻点，并检验实三次源法向。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, acb, acb_mat


def check_minimum_word_degree(h, pole, parts, adj, words, neighbours, boxes, equations):
    identity = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    symmetric = [identity,h,h*h,pole,h*pole+pole*h]
    cubic = h*pole*h
    orthogonal = []
    for word in symmetric+[cubic]:
        residual = word
        for axis in orthogonal:
            residual = residual - (axis*residual).trace().real*axis
        size = (residual*residual).trace().real.sqrt()
        assert size > arb("2.4e-7")
        orthogonal.append(residual/size)
    normal = orthogonal[-1]
    coeff = words["orbit_coefficients"](parts, adj, [normal])
    columns = words["trigonometric_columns"]
    for x,y in boxes:
        x,y = neighbours["refine"]((x,y,0),arb(".02"),equations,columns)
        value = (coeff*columns(x,y)[0])[0,0]
        assert value < arb("-.77")
        print("实三次基元法向在精确根处的配对",value.str(11))


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

    pole_row = acb_mat([[source[model["shell"]["soft"]["pole"],j] for j in range(3)]])
    pole = adj(pole_row)*pole_row
    pole_vector = adj(pole_row)
    columns_from_pole = (pole_vector, observable*pole_vector, observable**2*pole_vector)
    krylov = acb_mat([[col[i,0] for col in columns_from_pole] for i in range(3)])
    assert arb("3.08918e-7") < abs(krylov.det()) < arb("3.08919e-7")

    imag_axes = []
    for i,j in ((0,1),(0,2),(1,2)):
        outer = columns_from_pole[j]*adj(columns_from_pole[i])
        imag_axes.append((outer-adj(outer))/acb(0,2))
    imag_moments = words["orbit_coefficients"](parts, adj, imag_axes)
    cells = [(0,0,0)]
    real_leaves = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0,half.upper())
        vals = columns((2*i+1)*half+pad,(2*j+1)*half+pad)[0]
        residuals = imag_moments*vals
        if any(z.abs_lower()>arb("1e-14") for z in residuals.entries()):
            real_leaves += 1
        else:
            assert depth < 40
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("整轨道上极点三列的虚部至少一项超过 1e-14，终格数",real_leaves)

    commutator = acb(0,1)*(observable*pole-pole*observable)
    slopes = [(axis*commutator).trace().real for axis in imag_axes]
    assert abs(slopes[0]) > arb("3e-7")
    missing = [imag_axes[i] - (slopes[i]/slopes[0])*imag_axes[0] for i in (1,2)]
    equations = words["orbit_coefficients"](parts, adj, missing)
    phase_boxes = [(arb(x),arb(y)) for x,y in (
        ("1.657056","4.507081"),("5.450512","4.752972"),
        ("2.130443",".906317"),("5.058382",".640596"))]
    isolate = run_path(str(root/"ch39-real-cubic-udu-orbit.py"))["certify_four_roots"]
    leaves, deepest = isolate(equations, columns, phase_boxes)
    print("在两个壳算子的实代数上只加 i[H,P] 后恰有四个校准源，终格数与最高层",leaves,deepest)
    check_minimum_word_degree(observable, pole, parts, adj, words, neighbours, phase_boxes, equations)

    higher = acb(0,1)*(observable**2*pole-pole*observable**2)
    second_dir = [(axis*higher).trace().real for axis in imag_axes]
    direction = [slopes[(j+1)%3]*second_dir[(j+2)%3] - slopes[(j+2)%3]*second_dir[(j+1)%3]
                 for j in range(3)]
    divisor = max(z.abs_lower() for z in direction)
    assert divisor > arb("9.5e-14")
    cubic_harmonics = arb_mat([[z/divisor for z in direction]]) * imag_moments
    constant, dc, ds, xc, xs, yc, ys = cubic_harmonics.entries()
    y_size = (yc*yc+ys*ys).sqrt()
    beat = (dc*dc+ds*ds).sqrt()
    offset = abs(constant) + y_size
    amplitude = (xc*xc+xs*xs).sqrt() - beat
    assert xc - beat - offset > arb("7e-8")
    slope_floor = (amplitude**2-offset**2).sqrt()
    assert slope_floor > arb("1.2e-7")
    print("三次两条闭曲线的横向斜率下界", slope_floor.lower().str(10))




    pole_readings = [acb(0,1)*(axis*pole-pole*axis) for axis in axes]
    pole_moments = words["orbit_coefficients"](parts, adj, pole_readings)
    cells = [(0,0,0)]
    pencil_leaves = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0,half.upper())
        vals = columns((2*i+1)*half+pad,(2*j+1)*half+pad)[0]
        first, second = moments*vals, pole_moments*vals
        penalty = (arb(".16")*first[9,0].abs_upper())**2
        first_floor = sum(first[k,0].abs_lower()**2 for k in range(9))-penalty
        second_floor = sum(second[k,0].abs_lower()**2 for k in range(9))-penalty
        cross = sum(first[k,0]*second[k,0] for k in range(9))
        if first_floor>0 and second_floor>0 and first_floor*second_floor > cross.abs_upper()**2:
            pencil_leaves += 1
        else:
            assert depth < 18
            cells.extend((depth+1,2*i+di,2*j+dj) for di in (0,1) for dj in (0,1))
    print("任意归一壳循环与极点位置线性组合的交换子大于 0.16 倍源迹，终格数", pencil_leaves)




if __name__ == "__main__":
    verify()

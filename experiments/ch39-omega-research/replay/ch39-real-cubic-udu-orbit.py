"""固定十个夸克读数时，实二次词空间加入 UDU 后恰有四个正协方差。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def verify(replay):
    quadratic = run_path(str(replay / "ch39-two-zero-quadratic-source.py"))
    orbit = run_path(str(replay / "ch39-two-zero-orbit-condition.py"))
    model, up, down = orbit["calibration"]()
    parts, _ = orbit["orbit"](model, up, down)
    adjoint, trace = model["adj"], quadratic["trace"]
    up, down = [adjoint(v)*v for v in (up, down)]
    up, down = [v/trace(v).real for v in (up, down)]
    frame, _ = quadratic["normal_frame"](up, down)
    cubic = up*down*up
    # 与前两个法向的严格正交性来自迹循环。第三方向检验七维独立性。
    assert trace(frame[2]*cubic).real > arb(".00017")
    residuals = quadratic["orbit_coefficients"](parts, adjoint, frame[:2])
    columns = quadratic["trigonometric_columns"]
    centers = [(arb(x), arb(y)) for x, y in (
        ("1.93903", "3.46545"), ("1.97986", ".13985"),
        ("5.06480", "3.10882"), ("5.11728", ".22654"))]
    leaves, deepest = certify_four_roots(residuals, columns, centers)
    print("UDU 四解，原盒整条两相位轨道的终格数与最深层", leaves, deepest)


def certify_four_roots(residuals, columns, centers):
    assert len(centers) == 4
    radius = arb(1)/50
    pad = arb(0, radius.upper())
    assert all(radius < u < 2*arb.pi()-radius for p in centers for u in p)
    assert all(abs(x-u)>2*radius or abs(y-v)>2*radius
               for (x,y),(u,v) in combinations(centers, 2))
    identity = arb_mat([[1,0],[0,1]])
    for x, y in centers:
        values, derivative = columns(x, y)
        _, interval_derivative = columns(x+pad, y+pad)
        jacobian = residuals*derivative
        preconditioner = arb_mat([[jacobian[i,j].mid() for j in range(2)]
                                  for i in range(2)]).inv()
        shift = preconditioner*residuals*values
        error = identity-preconditioner*residuals*interval_derivative
        for i in range(2):
            assert (abs(shift[i,0])/radius + sum(abs(error[i,j]) for j in range(2))
                    < arb(".245"))

    cells = [(0, 0, 0)]
    leaves = deepest = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0, half.upper())
        x, y = (2*i+1)*half+pad, (2*j+1)*half+pad
        values = residuals*columns(x, y)[0]
        inside = any(abs(x-u)<radius and abs(y-v)<radius for u, v in centers)
        if inside or any(not z.contains(0) for z in values.entries()):
            leaves += 1
            deepest = max(deepest, depth)
        else:
            assert depth < 16
            cells.extend((depth+1, 2*i+di, 2*j+dj) for di in (0,1) for dj in (0,1))
    return leaves, deepest


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)

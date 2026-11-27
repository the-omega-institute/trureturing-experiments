"""在两个实际制备盒内检验固定十读数轨道与厄米二次源空间的交点。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def trace(matrix):
    return sum(matrix[i, i] for i in range(3))


def norm(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt()


def bracket(left, right):
    return acb(0, 1) * (left*right - right*left)


def normal_frame(up, down):
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    symmetric = [identity, up, down, up*up-up, down*down-down, (up-down)*(up-down)]

    def minor_row(matrix):
        return [matrix[0, 0].real, matrix[2, 2].real,
                matrix[0, 2].real, matrix[0, 2].imag,
                matrix[1, 2].real, matrix[1, 2].imag]

    minor = arb_mat([minor_row(matrix) for matrix in symmetric]).det()
    assert minor < arb("-8e-11")
    normals = [bracket(up, down), bracket(up*up, down), bracket(up, down*down)]
    frame = []
    for current in normals:
        for axis in frame:
            current -= trace(axis*current).real * axis
        length = norm(current)
        assert length > arb(".0001")
        frame.append(current / length)
    return frame, norm(normals[0])


def orbit_coefficients(parts, adjoint, frame):
    terms = [sum((f*adjoint(f) for f in parts), acb_mat(3, 3))]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        term = parts[i]*adjoint(parts[j])
        terms.extend((term+adjoint(term), acb(0, 1)*(term-adjoint(term))))
    return arb_mat([[trace(axis*term).real for term in terms] for axis in frame])


def trigonometric_columns(x, y):
    cx, sx = x.cos(), x.sin()
    cy, sy = y.cos(), y.sin()
    cd, sd = (x-y).cos(), (x-y).sin()
    value = arb_mat([[entry] for entry in (1, cd, sd, cx, sx, cy, sy)])
    derivative = arb_mat([[0, 0], [-sd, sd], [cd, -cd], [-sx, 0],
                          [cx, 0], [0, -sy], [0, cy]])
    return value, derivative


def root_boxes(coefficients, commutator_norm):
    centers = [(arb(x), arb(y), sign) for x, y, sign in (
        ("1.759345", "4.347800", -1), ("2.238043", "1.510924", 1),
        ("4.890319", "1.289542", 1), ("5.380130", "4.706125", -1))]
    radius = arb(1)/40
    pad = arb(0, radius.upper())
    assert all(radius < u < 2*arb.pi()-radius for x, y, _ in centers for u in (x, y))
    assert all(abs(x-u) > 2*radius or abs(y-v) > 2*radius
               for (x, y, _), (u, v, _) in combinations(centers, 2))
    identity = arb_mat([[1, 0], [0, 1]])
    residuals = arb_mat([[coefficients[i, j] for j in range(7)] for i in (1, 2)])
    for x, y, sign in centers:
        center_values, center_derivative = trigonometric_columns(x, y)
        box_values, box_derivative = trigonometric_columns(x+pad, y+pad)
        derivative = residuals*center_derivative
        preconditioner = arb_mat([[derivative[i, j].mid() for j in range(2)]
                                  for i in range(2)]).inv()
        shift = preconditioner*residuals*center_values
        error = identity-preconditioner*(residuals*box_derivative)
        for i in range(2):
            estimate = abs(shift[i, 0])/radius + sum(abs(error[i, j]) for j in range(2))
            assert estimate < arb(".46")
        gamma = (coefficients*box_values)[0, 0]/commutator_norm
        assert sign*gamma > 1
        print("盒中心与交换子系数", x.str(7), y.str(7), gamma.str(4))
    return centers, radius


def cover(coefficients, centers, radius):
    # 同一棵四分树同时认证到六维空间的距离和七维空间的交点数。
    cells = [(0, 0, 0)]
    leaves = deepest = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi()/2**depth
        pad = arb(0, half.upper())
        x = (2*i+1)*half + pad
        y = (2*j+1)*half + pad
        residual = coefficients*trigonometric_columns(x, y)[0]
        distance_ok = any(abs(v) > arb(1)/40 for v in residual.entries())
        root_free = any(not residual[k, 0].contains(0) for k in (1, 2))
        inside = any(abs(x-u) < radius and abs(y-v) < radius for u, v, _ in centers)
        if distance_ok and (root_free or inside):
            leaves += 1
            deepest = max(deepest, depth)
        else:
            assert depth < 15, (depth, i, j)
            cells.extend((depth+1, 2*i+di, 2*j+dj) for di in (0, 1) for dj in (0, 1))
    print("整环面到对称二次空间的距离大于 1/40，加入交换子后恰四解；终格数与最深层",
          leaves, deepest)


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    model, up, down = orbit["calibration"]()
    adjoint = model["adj"]
    parts, _ = orbit["orbit"](model, up, down)
    grams = [adjoint(matrix)*matrix for matrix in (up, down)]
    frame, commutator_norm = normal_frame(*(g/trace(g).real for g in grams))
    coefficients = orbit_coefficients(parts, adjoint, frame)
    centers, radius = root_boxes(coefficients, commutator_norm)
    cover(coefficients, centers, radius)


if __name__ == "__main__":
    verify()

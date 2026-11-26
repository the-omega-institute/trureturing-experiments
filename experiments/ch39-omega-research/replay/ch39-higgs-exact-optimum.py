"""Isolate the algebraic optimum of the maximum edge cost using six finite equations."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import acb, acb_mat, arb, arb_mat, ctx

previous = run_path(str(Path(__file__).with_name("ch39-higgs-maximum-edge.py")))
base = previous["base"]


def pencils():
    identity = base["identity"](5)
    c = (arb(5).sqrt()-1)/4
    root = acb(c, (1-c*c).sqrt())
    links = base["field"]()[:4]
    result = []
    for k in previous["momenta"]:
        parts = []
        for a, u in enumerate(links):
            f = u*root**k[a]
            parts.append(2*identity-f-base["adjoint"](f))
        result.append((parts[3], [part-parts[3] for part in parts[:3]] + [-identity]))
    return result


def partials(matrix, directions):
    # Determinants are multilinear in columns. Each replacement differentiates one column.
    columns = range(5)
    def replace(changes):
        value = acb_mat(matrix)
        for col, direction in changes:
            for row in columns:
                value[row, col] = direction[row, col]
        return value.det().real

    grad = [sum(replace([(j, op)]) for j in columns) for op in directions]
    hess = [[arb(0) for _ in directions] for _ in directions]
    for i, left in enumerate(directions):
        for j in range(i+1):
            right = directions[j]
            value = sum(replace([(a, left), (b, right)]) for a, b in permutations(columns, 2))
            hess[i][j] = hess[j][i] = value
    return matrix.det().real, grad, hess


def equations(point, data):
    values = []
    rows = []
    slopes = []
    derivatives = []
    for origin, dirs in data:
        matrix = sum((op*x for op, x in zip(dirs, point[:4])), origin)
        det, grad, hess = partials(matrix, dirs)
        values.append(det)
        rows.append(grad+[arb(0), arb(0)])
        slopes.append([-grad[a]/grad[3] for a in range(3)])
        derivatives.append([[(grad[a]*hess[3][j]/grad[3]-hess[a][j])/grad[3]
                             for j in range(4)] for a in range(3)])
    mixture = [*point[4:], 1-point[4]-point[5]]
    for a in range(3):
        values.append(sum(mixture[i]*slopes[i][a] for i in range(3)))
        row = [sum(mixture[i]*derivatives[i][a][j] for i in range(3)) for j in range(4)]
        rows.append(row+[slopes[i][a]-slopes[2][a] for i in range(2)])
    return arb_mat([[v] for v in values]), arb_mat(rows)


def stabilizer():
    links = []
    for columns, powers in zip(base["ROWS"], base["POWERS"]):
        link = sp.zeros(5)
        for row, (column, exponent) in enumerate(zip(columns, powers)):
            link[row, column] = sp.I**exponent
        links.append(link**5)
    eye = sp.eye(5)
    commutators = sp.Matrix.vstack(*[
        sp.kronecker_product(u, eye) - sp.kronecker_product(eye, u.T) for u in links
    ])
    assert commutators.rank() == 24


def verify():
    stabilizer()
    ctx.prec = 200
    data = pencils()
    center = [arb(v)/10**12 for v in (190228603135, 393680674736, 194536323233,
                                       508983786061, 701252351895, 57318735330)]
    # Work at exactly representable dyadic centers and with a dyadic box radius.
    center = [x.mid() for x in center]
    radius = arb(2)**-32
    box = [v+arb(0, radius) for v in center]
    f0, j0 = equations(center, data)
    _, jbox = equations(box, data)
    B = j0.inv()
    B = arb_mat([[B[i, j].mid() for j in range(6)] for i in range(6)])
    identity = arb_mat([[int(i==j) for j in range(6)] for i in range(6)])
    error = identity-B*jbox
    contraction = max(sum(abs(error[i, j]) for j in range(6)) for i in range(6))
    shift = max(abs(x) for x in (B*f0).entries())
    assert contraction < arb("0.00000037") and shift < arb("0.000000000000445")
    # Width of the isolated energy can be improved by applying the fixed-point map.
    enclosed = [abs(x) + contraction*radius for x in (B*f0).entries()]
    assert box[3] < arb("0.509") and 24*(radius+arb(2)**-190) < arb("0.000000006")
    assert min(box[0:3]+box[4:6]+[1-sum(box[:3]), 1-sum(box[4:])]) > 0
    energy = center[3]+arb(0, enclosed[3].upper())
    assert energy > arb("0.50898378606055") and energy < arb("0.50898378606145")
    print("Contraction, maximal first shift and optimal edge cost:", contraction.upper(),
          shift.upper(), energy)


if __name__ == "__main__":
    verify()

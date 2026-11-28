"""壳实协方差的两个带号 CP 镜像反例，在同一制备点保持九个精确读数。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat


pairs = tuple(combinations(range(3), 2))


def imaginary_equations(frame, moment, adj):
    diagonal = acb_mat([[moment[i,i].real if i == j else 0
                         for j in range(3)] for i in range(3)])
    constant = frame * diagonal * adj(frame)
    columns = [acb_mat([[frame[i,j]] for i in range(3)]) for j in range(3)]
    terms = []
    for i, j in pairs:
        term = (columns[i] * adj(columns[j])) * moment[i,j].conjugate()
        terms.extend((term + adj(term), acb(0,1) * (term - adj(term))))
    return [[term[i,j].imag for term in terms] for i, j in pairs], [
        [-constant[i,j].imag] for i, j in pairs
    ]

def exclusion(model, U, upper, lower, factor, pair_index):
    adj = model['adj']
    relative = upper * lower.inv()
    _, F = model["eig"](adj(relative) * relative)
    C = adj(U) * lower.inv() * F
    N = adj(U) * adj(lower) * F
    forward = adj(N) * factor
    backward = factor.inv() * C
    rows, constants = imaginary_equations(C, forward * adj(forward), adj)
    inv_rows, inv_constants = imaginary_equations(N, adj(backward) * backward, adj)
    Q, b = arb_mat(rows + inv_rows), arb_mat(constants + inv_constants)
    guess = Q.mid().solve(b.mid())
    radial = arb_mat([[guess[i,0].mid() if i in (pair_index, pair_index+1) else 0]
                      for i in range(6)])
    dual = Q.mid().transpose().solve(radial).transpose().mid()
    coefficients = (dual * Q).entries()
    cap = sum((abs(coefficients[j]).upper()**2 +
               abs(coefficients[j+1]).upper()**2).sqrt() for j in (0,2,4))
    margin = (dual*b)[0,0] - cap
    return Q, b, margin, C * forward.conjugate()

def readings(model, U, upper, lower, factor, scales):
    adj = model['adj']
    outputs = []
    for prepared, scale in zip((upper, lower), scales, strict=True):
        Y = arb(scale) * prepared * U * factor
        outputs.append(model["eig"](adj(Y) * Y))
    (up_levels, up_vectors), (down_levels, down_vectors) = outputs
    return [value.sqrt() for value in up_levels + down_levels] + model["mixing_angles"](up_vectors, down_vectors)



def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    model, up_box, down_box = load("two-zero-orbit-condition")["calibration"]()
    U, _ = load("shell-pencil-six-mass-window")["shell_frame"](model)
    at_zero = load("quadratic-gram-source")["at_zero"]
    adj = model["adj"]
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    eight_factor = acb_mat([[arb(s) for s in row.split()] for row in (
        "-1.0299955061154296 .0006739588194340176 -.00003734351146289303",
        "-.325480640747518 -.004401646063320289 .000005612213663910973",
        "-.7365532522102924 .0010026105016756317 .000049741117807932745",
    )])


    up_zero = at_zero(model, -4, 2345666)
    down_zero = at_zero(model, 2, 964890)
    for name, factor, scales, pair_index, cutoff, cp, square in (
        ("I", identity, ("1", "1"), 2, ".04", ("-.957063","-.957059"), ("1.0917804","1.0917805")),
        ("eight", eight_factor, ("19.604138072394083","2.774685099285122"),
         0, ".90", ("-.659669","-.659665"), ("20.91841","20.91842")),
    ):
        assert abs(factor.det()) > 0
        _, _, gap, _ = exclusion(model, U, up_box, down_box, factor, pair_index)
        assert gap > arb(cutoff)
        values = readings(model, U, up_box, down_box, factor, scales)
        assert all(0 < angle < arb.pi()/2 for angle in values[6:9])
        assert arb(cp[0]) < values[9] < arb(cp[1])
        if name == "eight":
            for i in (*range(6), 6, 8):
                assert abs(values[i] - model["ref"][i]) < arb(".096") * model["tol"][i]
        Q, b, _, mirror = exclusion(model, U, up_zero, down_zero, factor, pair_index)
        assert abs(Q.det()) > 0
        w = Q.solve(b)
        norm_sq = w[pair_index,0]**2 + w[pair_index+1,0]**2
        assert arb(square[0]) < norm_sq < arb(square[1])
        original = readings(model, U, up_zero, down_zero, factor, scales)
        conjugate = readings(model, U, up_zero, down_zero, mirror, scales)
        assert all(x.overlaps(y) for x,y in zip(original[:9],conjugate[:9],strict=True))
        assert original[9].overlaps(-conjugate[9]) and conjugate[9] > 0
        print(name, "实性分离余量下界", gap.lower().str(11), "原相位", values[9].str(12),
              "实际点的坐标平方和", norm_sq.str(15))
        print("实际点的十个读数", [value.str(13) for value in original])


if __name__ == '__main__':
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)

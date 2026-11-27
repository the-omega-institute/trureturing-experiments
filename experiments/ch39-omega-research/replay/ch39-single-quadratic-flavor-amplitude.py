"""单个二次振幅可精确恢复原源协方差，并在实际零点上给出有理复系数见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


def load(name):
    return run_path(str(Path(__file__).with_name(f"ch39-{name}.py")))


def cayley_numerator(t, identity):
    z, w = acb(t[0], t[1]), acb(t[2], t[3])
    h = acb_mat([[0, 0, -z.conjugate()], [0, 0, -w.conjugate()], [z, w, 0]])
    return (1 + sum(a*a for a in t)) * identity + 2 * (h + h*h)


def certify_exact(model, a, b, identity):
    quad = load("two-zero-quadratic-source")
    adj, trace, bracket = model["adj"], quad["trace"], quad["bracket"]
    up, down = [adj(v) * v for v in (a, b)]
    up, down = up / trace(up).real, down / trace(down).real
    words = (identity, up, down, up*up, down*down, up*down, down*up)
    # 若振幅的像在单一坐标轴上，落在另外两行的六个复元素必须全为零。
    for axis in range(3):
        rows = [(i, j) for i in range(3) if i != axis for j in range(3)]
        for omit in range(7):
            system = acb_mat([[w[i, j] for k, w in enumerate(words) if k != omit] for i, j in rows])
            assert abs(system.det()) > arb("7e-13")
    frame, _ = quad["normal_frame"](up, down)
    source = model["K"] * adj(model["K"])
    inverse_source = source.inv()
    spectral = [h / trace(h).real for h in (source, inverse_source)]
    minor = arb_mat([[trace(axis*h).real for h in spectral] for axis in frame[1:]]).det()
    assert minor < arb("-.00503")
    print("二次词空间与原协方差对易条件的二阶子式", minor.str(8))
    comm = bracket(up, down)
    normals = [trace(comm*comm).real * g - trace(comm*g).real * comm
               for g in (bracket(up*up, down), bracket(up, down*down))]
    test = [normal * model["K"] for normal in normals]

    def residual(t):
        value = cayley_numerator(t, identity)
        pair = [trace(s * value) for s in test]
        return arb_mat([[z.real] for z in pair] + [[z.imag] for z in pair])

    directions = [[int(i == j) for i in range(4)] for j in range(4)]

    def jacobian(t):
        # 每个剩余量为实二次多项式，单位步的中心差分恰等于导数。
        columns = []
        for e in directions:
            plus = residual([v + step for v, step in zip(t, e, strict=True)])
            minus = residual([v - step for v, step in zip(t, e, strict=True)])
            columns.append((plus - minus) / 2)
        return arb_mat([[column[i, 0] for column in columns] for i in range(4)])

    center = list(map(arb, ("1.0316517", "-.70582525", "-.4081726", ".7594895")))
    derivative = jacobian(center)
    midpoint = arb_mat([[derivative[i, j].mid() for j in range(4)] for i in range(4)])
    inverse = midpoint.inv()
    shift = max(abs(v).upper() for v in (inverse * residual(center)).entries())
    error = arb_mat(directions) - inverse * derivative
    slope = max(sum(abs(error[i, j]).upper() for j in range(4)) for i in range(4))
    hessian = [inverse * (jacobian([v + step for v, step in zip(center, e, strict=True)]) - derivative)
               for e in directions]
    curvature = max(sum(abs(h[i, j]).upper() for j in range(4) for h in hessian) for i in range(4))
    radius = arb(".001")
    lipschitz = slope + radius*curvature
    assert shift < arb(".000041") and lipschitz < arb(".00272")
    assert shift + radius*lipschitz < radius
    print("精确源协方差的四变量压缩，格心位移与导数上界", shift.str(8), lipschitz.str(8))


def certify_rational(model, identity):
    gram = load("quadratic-gram-source")
    adj = model["adj"]
    entries = (("-.20377732092973252", "3.0570597460114257"), ("-107.91386769430039", "67.170570375232302"),
               ("341.28168042780419", "-2541.1379299144978"), ("4689.5818512303358", "553.11029754044034"),
               ("1905.7257842147012", "16310.124089421444"), ("-3995.9544888158789", "-906.56188424282243"),
               ("-4011.5499464290533", "-984.75393463375042"))
    c = [acb(arb(x), arb(y)) for x, y in entries]
    pairs = ((-4, 2345666), (2, 964890))
    sites = [[gram["at_zero"](model, q, n+i) for i in range(-64, 65)] for q, n in pairs]
    grams = [[adj(a) * a for a in sector] for sector in sites]

    def amplitude(u, d):
        words = (identity, u, d, u*u, d*d, u*d, d*u)
        return sum((z*w for z, w in zip(c, words, strict=True)), acb_mat(3, 3))

    x = amplitude(grams[0][64], grams[1][64])
    yu, yd = sites[0][64]*x*arb("7.06539940386823251"), sites[1][64]*x
    (mu, lu), (md, ld) = [model["eig"](adj(y)*y) for y in (yu, yd)]
    readings = [v.sqrt() for v in mu+md] + model["mixing_angles"](lu, ld)
    errors = [abs(v-y)/w for v, y, w in zip(readings, model["ref"], model["tol"], strict=True)]
    maximum = max(e.upper() for e in errors)
    assert maximum < arb(".12352")
    floor = min(abs(amplitude(u, d).det()).lower() for u in grams[0] for d in grams[1])
    assert floor > arb("7.95e-6")
    print("十窗的最大半宽误差与 16641 处振幅行列式的模下界", maximum.str(8), floor.str(8))


if __name__ == "__main__":
    model, up, down = load("two-zero-orbit-condition")["calibration"]()
    identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    certify_exact(model, up, down, identity)
    certify_rational(model, identity)

"""固定两种荷的算术枢轴，按谱序是否相同分配两种插值，核验形式交汇的符号。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    probes = {}
    for q, levels, direct, _ in native["shell_weights"](root):
        if q in (-4, 2):
            pivot = 0 if q == -4 else 2
            nodes = [v/levels[-1] for v in levels[:-1]]+[arb(1)]
            probes[q] = list(zip(nodes, [v.real for v in direct[pivot]], strict=True))
    edges = []
    for (i, (s, u)), (j, (t, v)) in product(enumerate(probes[-4]), enumerate(probes[2])):
        edges.append((s*t, (s+t)/2-s*t, int(i == j), u*v))
    rebuild = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (central,) = rebuild(root, (2_345_666,), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    ctx.prec = 192
    pi = arb.pi()
    beta = [v/(2*pi) for v in model["beta"]]
    d = (-arb(115)/218, arb(333)/218, -1)
    dot = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))["dot"]
    gauge = run_path(str(root/"ch39-signed-gravity-unification.py"))
    coeff = gauge["gauge_weights"](model)
    projected = [dot(d, row) for row in coeff]
    heat = lib["expint"]

    def calculate(point, *, derivative=False):
        c, h, l = point
        jac = arb_mat(3, 3)
        k, w = arb(0), arb(0)
        intercept = list(initial)
        for base, difference, same, weight in edges:
            v = base+difference*(l if same else h)
            x = c*v
            if derivative:
                exp_values = [(-z/x).exp() for _, _, z in central]
                rates = [sum(a*y for (_, a, _), y in zip(central, exp_values, strict=True))/(2*pi),
                         -x*sum(g*y for (g, _, _), y in zip(central, exp_values, strict=True))/(4*pi),
                         -dot(projected, exp_values[1:])/(4*pi*x)]
                for i, rate in enumerate(rates):
                    jac[i, 0] += weight*rate*v
                    jac[i, 1+same] += weight*rate*c*difference
            else:
                kn, wn = single["moments"](central, heat, x)
                k += weight*kn
                w += weight*wn
                e1 = [heat(1, z/x) for _, _, z in central[1:]]
                for i in range(3):
                    intercept[i] -= weight*sum(b[i]*e for b, e in zip(coeff, e1, strict=True))/(4*pi)
        if derivative:
            return jac
        return arb_mat([[k-1], [w-arb("1e-122")], [dot(d, intercept)]]), intercept

    center = arb_mat([[arb(x)] for x in (".00560510475479", ".900144624274", ".15840391850")])
    for _ in range(2):
        y, _ = calculate(center.entries())
        J = calculate(center.entries(), derivative=True)
        center = arb_mat([[v.mid()] for v in (center-J.solve(y)).entries()])
    radius = arb("1e-13")
    box = [v+arb(0, radius.upper()) for v in center.entries()]
    J = calculate(center.entries(), derivative=True)
    assert not J.det().contains(0)
    inverse = J.inv()
    residual, _ = calculate(center.entries())
    shift = max(abs(v).upper() for v in (inverse*residual).entries())
    identity = arb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    defect = identity-inverse*calculate(box, derivative=True)
    norm = max(sum(abs(defect[i,j]) for j in range(3)).upper() for i in range(3))
    print("三参数根盒中心、压缩率与中心位移", *[v.str(18) for v in center.entries()],
          norm.str(10), shift.str(10),flush=True)
    assert norm*radius+shift < radius
    _, intercept = calculate(box)
    logscale = (intercept[0]-intercept[1])/(beta[0]-beta[1])
    common = intercept[0]-beta[0]*logscale
    print("形式交汇的 ln(M/M_Z)、统一逆耦合及三个截距", logscale.str(15), common.str(15),
          *[v.str(11) for v in intercept])
    assert 10000 < logscale < 10200 and -2500 < common < -2300


if __name__ == "__main__":
    verify()

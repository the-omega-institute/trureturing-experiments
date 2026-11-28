"""两个原生枢轴的已定实权，借共同谱标尺与和积插值匹配两个引力矩。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    _, levels, direct, _ = next(item for item in native["shell_weights"](root) if item[0] == -3)
    unit_levels = [v/levels[-1] for v in levels[:-1]]+[arb(1)]
    first = zip(unit_levels, direct[0], strict=True)
    second = zip(unit_levels, direct[1], strict=True)
    edges = [(s*t, (s+t)/2-s*t, u.real*v.real)
             for (s, u), (t, v) in product(first, second)]
    assert sum(w for _, _, w in edges).contains(1)

    rebuild = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    single, lib, (central, alternative) = rebuild(root, (2_345_666, 27_804_036), tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = run_path(str(root/"ch39-signed-gravity-unification.py"))
    ctx.prec = 192
    pi = arb.pi()
    d = (-arb(115)/218, arb(333)/218, -1)
    dot = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))["dot"]
    projected = [dot(d, g) for g in gauge["gauge_weights"](model)]
    heat = lib["expint"]
    def gravity(c, mix, rows):
        k, w = arb(0), arb(0)
        for base, difference, weight in edges:
            x = c*(base+mix*difference)
            kval, wval = single["moments"](rows, heat, x)
            k += weight*kval
            w += weight*wval
        return arb_mat([[k-1], [w-arb("1e-122")]])

    def jacobian(c, mix):
        dc_k, dt_k, dc_w, dt_w = (arb(0) for _ in range(4))
        for base, difference, weight in edges:
            node = base+mix*difference
            x = c*node
            ex = [(-z/x).exp() for _, _, z in central]
            dk = sum(a*y for (_, a, _), y in zip(central, ex, strict=True))/(2*pi)
            dw = -x*sum(g*y for (g, _, _), y in zip(central, ex, strict=True))/(4*pi)
            dc_k += weight*dk*node
            dt_k += weight*dk*c*difference
            dc_w += weight*dw*node
            dt_w += weight*dw*c*difference
        return arb_mat([[dc_k, dt_k], [dc_w, dt_w]])

    def gauge_residual(c, mix, rows):
        residual = dot(d, initial)
        for base, difference, weight in edges:
            x = c*(base+mix*difference)
            residual -= weight*sum(b*heat(1, z/x) for b, (_, _, z)
                                   in zip(projected, rows[1:], strict=True))/(4*pi)
        return residual

    c, t = arb(".0842439178631902714"), arb(".900565026742734474")
    radius = arb("1e-17")
    box = (c+arb(0, radius.upper()), t+arb(0, radius.upper()))
    center_matrix = jacobian(c, t)
    assert not center_matrix.det().contains(0)
    inverse = center_matrix.inv()
    defect = arb_mat([[1, 0], [0, 1]])-inverse*jacobian(*box)
    lipschitz = max(sum(abs(defect[i, j]) for j in range(2)).upper() for i in range(2))
    shift = max(abs(v).upper() for v in (inverse*gravity(c, t, central)).entries())
    assert lipschitz < arb("2e-14") and shift+lipschitz*radius < radius
    central_projection = gauge_residual(*box, central)
    dk, dw = gravity(*box, alternative).entries()
    alt_projection = gauge_residual(*box, alternative)
    assert arb("-18.95022") < central_projection < arb("-18.95020")
    assert arb("2.5069e-9") < dk < arb("2.5070e-9")
    assert arb("2.4305e-10") < dw < arb("2.4306e-10")
    assert arb("-18.4757") < alt_projection < arb("-18.4755")
    print("中央根盒中心与半宽、映射斜率上界", c.str(18), t.str(18), radius.str(6), lipschitz.str(10))
    print("中央规范投影；替代编号的牛顿、真空和规范残差",
          *[v.str(13) for v in (central_projection, dk, dw, alt_projection)])


if __name__ == "__main__":
    verify()

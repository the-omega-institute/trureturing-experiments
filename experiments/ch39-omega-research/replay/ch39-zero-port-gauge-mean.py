"""规范均值固定时，以一个实对角端口检验全谱非相干混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_poly

root = Path(__file__).parent
single = run_path(str(root/"ch39-zero-single-port.py"))
gauge = run_path(str(root/"ch39-high-zero-planck-window.py"))


def verify():
    ctx.prec = 260
    q, *nums = single["polynomials"]()
    offset = fmpq(gauge["OFFSET"].numerator, gauge["OFFSET"].denominator)
    k = offset-fmpq(1, 4)
    den = fmpq_poly([1, -2*k])
    domain = arb(-1).union(arb(0))/200
    assert acb.zeta_zero(1).imag > 10
    def evaluate(poly, x):
        out = arb(0)
        for coefficient in reversed(poly.coeffs()):
            out = out*x+arb(coefficient)
        return out
    assert evaluate(q, domain) > 0 and evaluate(den, domain) > 0
    heights = [acb.zeta_zero(n).imag for n in (11430, 11431, 11432)]
    xs = [-1/(2*(h*h+arb(1)/4)) for h in heights]
    gl, gy, gr = [1/(h*h+arb(offset)) for h in heights]
    weight = (gy-gr)/(gl-gr)
    assert 0 < weight < 1
    gaps = []
    for num in nums:
        slope = num.derivative()*q-num*q.derivative()
        curvature = den*(slope.derivative()*q-2*slope*q.derivative())-4*k*slope*q
        assert 0 not in evaluate(curvature, domain)
        l, y, r = [evaluate(num,x)/evaluate(q,x) for x in xs]
        gaps.append(abs(r-y+weight*(l-r)))
        if num == nums[6] and len(gaps) == 7:
            # 半个负二阶导数在规范变量 X 上正，最小值下界取自 x=0。
            product = -den**3*curvature
            assert evaluate(product.derivative()*q-3*product*q.derivative(),domain)<0
            scale = heights[1]**2+arb(offset)
            half_curvature = evaluate(product,arb(0))/(8*scale**2*evaluate(q,arb(0))**3)
            kappa = arb("8.59405615e-23")
            assert half_curvature > kappa
            tangent = -evaluate(den,xs[1])**2*evaluate(slope,xs[1])/(2*scale*evaluate(q,xs[1])**2)
            required = kappa*arb(".0068")
            gauge_coefficient = -tangent-kappa*arb(".01018")
            assert gauge_coefficient > 0 and required > arb("5.84395818e-25")
            assert required/gauge_coefficient > arb("1.66614681e-9")
            mass = (required-gauge_coefficient*arb("1e-9")-arb("1e-25"))/(kappa*arb(".013785"))
            assert mass>arb(".112")
            print("Flavor coefficient box: port deficit, relative gauge error, example share",
                  required.str(15), (required/gauge_coefficient).str(15), mass.str(15))
    best = max(range(19), key=lambda i: float(gaps[i]))
    tied = [j for j in range(19) if nums[j] == nums[best]]
    assert gaps[best] > 0 and all(gaps[i] < gaps[best] for i in range(19) if i not in tied)
    num = nums[6]
    left_val, target_val, right_val = [evaluate(num, x)/evaluate(q, x) for x in xs]
    x_left, x_right = gl/gy-1, gr/gy-1
    chord_slope = (left_val-right_val)/(x_left-x_right)
    y_at_zero = right_val-chord_slope*x_right
    signed_gap = target_val-y_at_zero
    gauge_cross = signed_gap/chord_slope
    assert x_right < gauge_cross < 0 < x_left and signed_gap > 0
    print("Relative gauge offset first masking a zero-free target",
          gauge_cross.str(23), "slope", chord_slope.str(23))
    print("Gauge-matched neighbor pair, left weight", weight.str(22),
          "best diagonal ports", tied, "sharp gap", gaps[best].str(22))


if __name__ == "__main__":
    verify()

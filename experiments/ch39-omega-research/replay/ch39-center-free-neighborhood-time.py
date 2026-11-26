"""由规范与端口读数求不使用 11431 时若干邻域占时的锐下界。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb

budget = run_path(str(Path(__file__).with_name("ch39-global-three-zero-clock.py")))
study, box, support, edges = (budget[key] for key in ("study", "box", "support", "edges"))
center = 11431


def model():
    cost, _, scale, *_ = support["prev"]["setup"]()

    @cache
    def point(n):
        height = acb.zeta_zero(n).imag
        x = arb(0) if n == center else study["center"] / (height**2 + study["offset"]) - 1
        return height, x, arb(0) if n == center else cost(x)[0]

    reference = [point(n) for n in (11431, 11432, 11415, 11430, 11372, 11304)]
    corner = [v + sign * arb(2)**-28 for v, sign in zip(box["centers"], (1, -1, 1, -1))]
    rows = [box["flavor_data"](study["native"]["plane_response"](h), corner)[0]
            for h, _, _ in reference]
    weights, _ = edges["candidate"]([x + scale*c/2 for _, x, c in reference], rows, arb(1)/3)
    eps = scale * support["dot"](weights, [c for _, _, c in reference]) / 2
    assert scale > 0 and eps > 0
    return point, scale, eps


def verify_radius(point, scale, eps, radius, expected, lower_digits):
    left, right = center - radius - 1, center + radius + 1
    _, a, ca = point(left)
    _, b, cb = point(right)
    slope = (ca - cb) / (a - b)
    assert a > b and 1 + scale*slope > 0

    def chord(x):
        return cb + slope*(x - b)

    def deficit(n):
        _, x, c = point(n)
        return chord(x) - c

    peak = deficit(expected)
    assert peak > 0
    assert all(peak > deficit(n) for n in range(left + 1, right) if n not in (center, expected))

    bound = (chord(-eps) - 2*eps/scale) / peak
    assert arb(lower_digits) / 10**12 < bound < arb(lower_digits + 1) / 10**12
    _, u, cu = point(expected)
    left_weight = (-eps - b - bound*(u - b)) / (a - b)
    right_weight = 1 - bound - left_weight
    assert left_weight > 0 and right_weight > 0
    gauge = left_weight*a + bound*u + right_weight*b
    port = gauge + scale*(left_weight*ca + bound*cu + right_weight*cb)
    assert 0 in gauge + eps and 0 in port - eps
    # 端点外由原端口代价的严格凸性统一覆盖。三点权重随两读数连续。
    print(radius, (left, expected, right), "mass", bound.str(20),
          "outer weights", left_weight.str(15), right_weight.str(15))


if __name__ == "__main__":
    inputs = model()
    for radius, label, digits in (
        (2, 11432, 221468071237),
        (3, 11432, 525512633705),
        (4, 11430, 658985022763),
        (5, 11430, 757645186302),
        (8, 11432, 917482341597),
    ):
        verify_radius(*inputs, radius, label, digits)

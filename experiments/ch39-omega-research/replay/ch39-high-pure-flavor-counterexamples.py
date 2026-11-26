"""直接按索引计算高零点，检验四个可独立扰动的味系数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def flavor_test():
    model = run_path(str(Path(__file__).with_name("ch39-fixed-local-flavor.py")))
    coefficients = [arb(n) / 2048 + arb(0, arb(2)**-28) for n in model["codes"]]
    radii = [model["shell"]["sine_radius"]] * 3 + [model["shell"]["phase_radius"]]

    def normalized_errors(zero):
        assert zero.real == arb(1) / 2
        response = model["native"]["plane_response"](zero.imag)
        errors = model["errors"](model["mixing"](response, coefficients))
        return [value / radius for value, radius in zip(errors, radii, strict=True)]

    return normalized_errors


def verify():
    ctx.prec = 160
    evaluate = flavor_test()
    root = Path(__file__).parent
    finite = run_path(str(root / "ch39-high-pure-rational-clocks.py"))
    phases = run_path(str(root / "ch39-arc-relative-determinant.py"))["phases"]
    for n, bound in ((11431, ".952614"), (79545, ".773225"), (97370, ".831528")):
        zero = acb.zeta_zero(n)
        worst = max(value.upper() for value in evaluate(zero))
        assert worst < arb(bound) < 1
        actual = phases(zero.imag)
        for distance, (sign, code) in zip((724, 804), finite["clock_codes"][n], strict=True):
            z = actual[distance]
            assert sign*z.real > 0
            slope = sign*z.imag / (1 + sign*z.real)
            assert abs(slope*(1 << 40) - code) < arb(1)/2
        print("Single zero", n, "height", zero.imag.str(15), "maximum error / tolerance", worst.str(12))


if __name__ == "__main__":
    verify()

"""用全高斯有理的固定门近似两组双壳解，并直接计算三个壳的引力读数。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb


DENOMINATOR = 2**20
SEEDS = (
    (-1255588, 572840, -3930214, -267745, -327121, -267745, 756531,
     -1633060, -1115725, 976850, 412747, -3961832),
    (-1500734, 684187, -6042665, -190470, -106864, -57401, 684187,
     -3193612, -1551825, 384618, 391683, -15792752),
)


def rational_circle(numerator):
    q = Fraction(numerator, DENOMINATOR)
    c, s = (1 - q*q)/(1 + q*q), 2*q/(1 + q*q)
    assert c*c + s*s == 1
    return acb(arb(c.numerator)/c.denominator, arb(s.numerator)/s.denominator)


def apply_gates(station, layers, seed):
    gates = []
    for k in range(3):
        mix, alpha, beta, scalar = map(rational_circle, seed[4*k:4*k+4])
        gates.append((scalar*alpha*mix.real, -scalar*beta.conjugate()*mix.imag,
                      scalar*beta*mix.imag, scalar*alpha.conjugate()*mix.real))
    result = station * 1
    for k, pairs in layers:
        a, b, c, d = gates[k]
        for near, far in pairs:
            for j in range(2):
                u, v = result[near, j], result[far, j]
                result[near, j], result[far, j] = a*u+b*v, c*u+d*v
    return result


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    layers = lib("two-shell-ordered-bond-gates")["LAYERS"]
    geometry = lib("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    for number in (2_306_888_100_092, 2_306_888_100_093, 2_306_888_100_094):
        _, _, _, basis, station, weights, heat, field = geometry(
            root, fifth, columns, levels, number=number
        )
        for index, seed in enumerate(SEEDS):
            projection = basis.conjugate().transpose() * apply_gates(station, layers, seed)
            residual, _ = gravity(projection, [], levels, weights, heat, field, True)
            k, w = residual[0, 0] + 1, residual[1, 0] + arb("1e-122")
            print(number, index, "K、W、最大绝对误差上界",
                  k.str(13), w.str(13), max(v.abs_upper() for v in residual.entries()).str(14),
                  flush=True)
            if number == 2_306_888_100_094:
                assert (k > 0 and w < 0) if index == 0 else (k < 0 and w > 0)
            else:
                assert max(v.abs_upper() for v in residual.entries()) < arb(".000032")


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))

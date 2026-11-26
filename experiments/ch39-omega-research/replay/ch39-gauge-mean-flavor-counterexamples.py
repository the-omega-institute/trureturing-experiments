"""按谱权重平均规范靶值时，两个零点可同时重现规范均值并通过味容差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 192
    root = Path(__file__).parent
    flavor = run_path(str(root / "ch39-fixed-local-flavor.py"))
    gauge = run_path(str(root / "ch39-high-zero-planck-window.py"))
    coefficients = [arb(k)/2048 + arb(0, arb(2)**-28) for k in flavor["codes"]]
    offset = gauge["q"](gauge["OFFSET"])
    zeros = {n: acb.zeta_zero(n) for n in (1, 11431, 79545, 97370)}
    reduced = {n: 1/(zero.imag**2 + offset) for n, zero in zeros.items()}
    shell = flavor["shell"]
    tolerances = [shell["sine_radius"]]*3 + [shell["phase_radius"]]

    def observables(n):
        response = flavor["native"]["plane_response"](zeros[n].imag)
        matrix = flavor["mixing"](response, coefficients)
        j = (matrix[0,0]*matrix[1,0].conjugate()*matrix[1,1]*matrix[0,1].conjugate()).imag
        return [abs(matrix[i,k])**2 for i,k in ((0,1), (1,2), (0,2), (1,0))] + [j]

    low = observables(1)
    for n, error_ceiling in ((79545, ".773255"), (97370, ".831559")):
        weight = (reduced[11431] - reduced[n]) / (reduced[1] - reduced[n])
        assert 0 < weight < 1
        a, b, u, d, j = [weight*x + (1-weight)*y for x, y in zip(low, observables(n), strict=True)]
        ss, tt = a/(1-u), b/(1-u)
        real = (1-u)*(d-ss*(1-tt)-(1-ss)*tt*u)/2
        values = [ss.sqrt(), tt.sqrt(), u.sqrt(), acb(real, j).arg()]
        worst = max((abs(value-center)/radius).upper() for value, center, radius in
                    zip(values, shell["centers"], tolerances, strict=True))
        variance = (reduced[1]-reduced[11431])*(reduced[11431]-reduced[n])/reduced[11431]**2
        print("Two-point mean", 1, n, "low weight", weight.str(12), "max CKM fraction", worst.str(12),
              "relative variance", variance.str(12), flush=True)
        assert worst < arb(error_ceiling) < 1
        assert variance > 516689


if __name__ == "__main__":
    verify()

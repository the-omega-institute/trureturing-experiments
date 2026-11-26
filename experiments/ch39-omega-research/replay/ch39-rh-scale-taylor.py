"""Bound the change in a finite five-tick calibration after a common pole rescaling."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def row_norm(matrix):
    return max(sum(abs(matrix[i, j]).upper() for j in range(matrix.ncols())).upper()
               for i in range(matrix.nrows()))


def main():
    setup = run_path(str(Path(__file__).with_name("ch39-rh-shell-smallest-tick.py")))
    ctx.prec, ctx.cap = 1536, 3
    heights = [acb.zeta_zero(i).imag for i in setup["support"]]
    values = [[1/(1+j*j/(g*g)) for g in heights] for j in range(42)]
    logs = [arb(x).log() for x in (2, 3, arb(3)/2, 6)]
    clock = [[(5*g*log).cos() for g in heights] for log in logs]
    inverse = arb_mat(values + clock).inv()
    target = setup["prior"]["read"](42) + [arb(0)]*4
    weights = (inverse*arb_mat([[v] for v in target])).entries()
    radius = arb(2)**-166
    # For v=1/(1+j^2/g^2), the scaled reading is v/[1+(1-v)*delta].
    # Finite geometric division gives the remainder below. Clock rows stay fixed.
    rows = values
    bound = arb(0)
    for power in range(1, 17):
        rows = [[a*(v-1) for a, v in zip(row, original)]
                for row, original in zip(rows, values)]
        coeff = inverse*arb_mat(rows + [[0]*46 for _ in range(4)])
        bound += row_norm(coeff)*radius**power
    bound += 46*row_norm(inverse)*radius**17/(1-radius)
    assert bound < arb(".002492")
    assert bound*max(weights) < min(weights)*(1-bound)
    print("Exact 42 readings and four phases have positive weights throughout |scale-1| <= 2^-166.")


if __name__ == "__main__":
    main()


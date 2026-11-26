"""仅用指定交叉块排除五维与六维 Gram 补全。"""
from pathlib import Path
from runpy import run_path

from flint import fmpq

previous = run_path(str(Path(__file__).with_name("ch39-local-memory-rank.py")))
submatrix = previous["submatrix"]


def verify():
    _, _, _, g0, g1, g2 = previous["MEM"]["coefficients"]()
    delta = g2 - g1 * g0.solve(g1)
    local = previous["LOCAL"]
    points = sorted(local["shell_grid"](1) | {(17, 17, 17)})

    def forced(cutoff_exclusive, i, j, rows, cols):
        blocks = [(rows, cols), ([i], cols), (rows, [j])]
        assert all(local["squared"](points[u], points[v]) >= cutoff_exclusive
                   for left, right in blocks for u in left for v in right)
        return (submatrix(delta, [i], cols)
                * submatrix(delta, rows, cols).solve(submatrix(delta, rows, [j])))[0, 0]

    diagonal = forced(292, 0, 0, (1, 3, 5, 7, 9), (10, 13, 16, 17, 18))
    assert diagonal < 0
    a = forced(216, 6, 10, (0, 1, 2, 3, 5, 7), (13, 14, 15, 16, 17, 18))
    b = forced(216, 6, 10, (0, 1, 2, 3, 17, 18), (5, 9, 12, 13, 14, 15))
    assert fmpq(336, 10**11) < a < fmpq(337, 10**11)
    assert fmpq(182, 10**11) < b < fmpq(183, 10**11)
    print("Below 292, rank five forces a negative diagonal", float(diagonal))
    print("Below 216, rank six forces incompatible entries", float(a), float(b))


if __name__ == "__main__":
    verify()

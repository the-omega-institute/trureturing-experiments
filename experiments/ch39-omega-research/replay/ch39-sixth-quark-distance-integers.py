"""仅用有理数与整数惯性核验第六壳四个轻夸克窗的射影距离下界。"""
from fractions import Fraction
from math import floor
from pathlib import Path
from runpy import run_path


def verify(root):
    rows = [line.split() for line in (root/"ch39-sixth-quark-distance-cells.tsv").read_text().splitlines()]
    scale = 10**14
    assert rows[:4] == [["scale", str(scale)], ["shell", "2306888100091"], ["sites", "14", "8", "7", "6"],
                        ["dw_upper", "281000", str(scale)]] and len(rows) == 82
    dw = Fraction(int(rows[3][1]), scale)
    inertia = run_path(str(root/"ch39-fourth-shell-strange-integers.py"))
    shifts = [Fraction(2)**k*m/4 for k in range(-32, 10) for m in (4, 5, 6, 7)]
    indices = [(i, j) for i in range(4) for j in range(i, 4)]
    tests = ((-4, 2, (".0669", ".00000629"), (".0059", ".00316")),
             (2, 1, (".0118", ".0000137"), (".0632", ".000271")))
    for block, (charge, power, *windows) in enumerate(tests):
        data = rows[4+39*block:4+39*(block+1)]
        assert data[0][:2] == ["train", str(charge)] and len(data[0]) == 21
        source = [Fraction(int(v), scale) for v in data[0][2:]]
        assert all(0 <= v < 1 for v in source)
        resolvents = []
        projectors = []
        for i, (prob, proj) in enumerate(zip(data[1:20], data[20:], strict=True)):
            assert prob[:3] == ["prob", str(charge), str(i)] and len(prob) == 22
            upper = [Fraction(int(n), scale) for n in prob[3:]]
            assert min(upper) >= 0 and sum(upper) >= 1
            resolvents.append([sum(a/(t+s) for a, t in zip(upper, source, strict=True)) for s in shifts])
            assert proj[:3] == ["projector", str(charge), str(i)] and len(proj) == 23
            numbers = list(map(int, proj[3:]))
            pairs = list(zip(numbers[::2], numbers[1::2], strict=True))
            assert all(v[1] == 0 for (j, k), v in zip(indices, pairs, strict=True) if j == k)
            projectors.append(pairs)

        for j, (radius, ceiling) in enumerate(windows):
            distance = power*(Fraction(radius)+2*dw)
            overlap = 1-distance**2/2
            assert overlap > 0
            lower = [max(Fraction(0), max(overlap**2/r-s for r, s in zip(row, shifts, strict=True)))
                     for row in resolvents]
            integers = [floor(a*scale) for a in lower]
            assert all(0 <= n <= scale for n in integers)
            cells = {key: tuple(sum(w*r[k][part] for w, r in zip(integers, projectors, strict=True))
                                for part in (0, 1)) for k, key in enumerate(indices)}
            # 每个投影的实部、虚部逐项误差 <1/S，十九项的算子误差 <114/S。
            threshold = (Fraction(ceiling)+Fraction(114, scale))*scale**2
            assert threshold.denominator == 1
            count = inertia["negative_count"](inertia["realify"](4, cells), int(threshold))
            assert count <= 2*j
            print("荷、三代能序、射影排除半径、窗顶、实化负惯性", charge, j+1, radius, ceiling, count)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)

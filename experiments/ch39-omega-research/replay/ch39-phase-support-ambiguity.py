"""Positive off-roster alternative with all four shell means and 67 exact time rows."""
from fractions import Fraction
from itertools import permutations
from math import comb
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
upper = run_path(str(root / "ch39-phase-support-identification.py"))
base, B, multiply, rounded = [upper[k] for k in ("base", "B", "multiply", "rounded")]
outside = (30, 33, 42, 46, 57, 59, 77, 82, 85, 100, 117, 124, 171, 196, 235, 447,
           453, 469, 483, 486, 2318, 2327)


def determinant(matrix):
    result = B(0)
    for columns in permutations(range(4)):
        inversions = sum(j > k for i, j in enumerate(columns) for k in columns[i+1:])
        term = multiply(*(row[k] for row, k in zip(matrix, columns)))
        result += term if inversions % 2 == 0 else B(-term.hi, -term.lo)
    return rounded(result)


def direction():
    roster = set(upper["time"]["roster"]()[0])
    ids = sorted(roster | set(outside))
    offsets = [*base["OFFSETS"], *range(43, 66)]
    assert len(ids) == 72 and len(roster) == 50 and len(set(offsets)) == 67
    prefix = base["data"]("ch39-tick-separators.json")["prefix"]
    x = [upper["height_cell"](prefix[i-1][0]) for i in ids]
    phases = base["observations"](ids)[0][-4:]
    h = [rounded((v-x[-1])/x[-1]) for v in x]
    beta = [rounded(multiply(*(B(1)+B(j*j)*v for j in offsets if j)) /
                    multiply(*(v-y for k, y in enumerate(x) if k != i))) for i, v in enumerate(x)]
    matrix = []
    for row in phases:
        f = [B(Fraction(lo, base["SCALE"]), Fraction(hi, base["SCALE"])) for lo, hi in row]
        powers = list(beta)
        moments = []
        for degree in range(5):
            moments.append(rounded(sum((multiply(p, y) for p, y in zip(powers, f)), B(0))))
            powers = [multiply(p, v) for p, v in zip(powers, h)]
        matrix.append(moments)
    # Exact signed maximal minors annihilate all four phase rows.
    coefs = [B((-1)**j)*determinant([r[:j]+r[j+1:] for r in matrix]) for j in range(5)]
    # At x=-1/a, the local coordinate is -1-t for t>0. Positivity of the
    # translated coefficients tests every new positive spacing at once.
    for k in range(5):
        coefficient = sum((B((-1)**j * comb(j,k))*coefs[j] for j in range(k,5)), B(0))
        assert coefficient.lo > 0
    signs = {}
    weights = {}
    for i, v, b in zip(ids, h, beta):
        poly = B(0)
        for c in reversed(coefs):
            poly = rounded(poly*v + c)
        w = b*poly
        weights[i]=w
        assert w.lo*w.hi > 0
        signs[i] = 1 if w.lo > 0 else -1
    assert len({signs[i] for i in outside}) == 1
    sign = signs[outside[0]]
    assert any(sign * signs[i] < 0 for i in roster)
    total = sum((weights[i] for i in outside), B(0))
    return {i: rounded(w/total) for i,w in weights.items()}, x, coefs, total, offsets


def main():
    weights, _, _, _, _ = direction()
    print("67 time readings and all four shell means agree with a positive perturbation.")
    print("New sites:", len(outside), "highest index:", max(outside))
    print("Every positive time offset square outside the original set distinguishes this pair.")
    print("Decreased old weights:", sum(w.hi < 0 for w in weights.values()))


if __name__ == "__main__":
    main()

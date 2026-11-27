"""在严格三循环对称下，以有理线性不等式排除最后一种六壳选线。"""

import gzip
import json
from fractions import Fraction
from math import ceil, floor
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb

BOX_SCALE = 2**24
ROW_SCALE = 10**6
DUAL_SCALE = 10**11
DIM = 36


def spectral_rows(model, cuts):
    vectors = [
        [shell["main"][i, k] for k in range(6)] for shell in model for i in shell["big"]
    ]
    projectors = [
        acb_mat([[a.conjugate() * b for b in vector] for a in vector])
        for vector in vectors
    ]
    rows = []
    for coefficients, rhs in cuts:
        assert len(coefficients) == DIM and all(type(c) is int for c in coefficients)
        assert type(rhs) is int
        matrix = sum(
            (c * p for c, p in zip(coefficients, projectors, strict=True) if c),
            acb_mat(6, 6),
        )
        for n in range(1, 7):
            minor = acb_mat(
                [[int(i == j) * rhs - matrix[i, j] for j in range(n)] for i in range(n)]
            )
            assert minor.det().real > 0
        rows.append((coefficients, rhs))
    for j in range(6):
        for sign in (1, -1):
            rows.append(
                (
                    [sign * ROW_SCALE if i // 6 == j else 0 for i in range(DIM)],
                    sign * ROW_SCALE,
                )
            )
    return rows


def moment_planes():
    for shell in range(6):
        for power, index in ((2, 0), (4, 2)):
            for tick in range(5):
                yield shell, "tangent", power, index, tick
            yield shell, "secant", power, index, 0
        for numerator in (9, 16, 25, 36, 50, 65):
            for sign in (1, -1):
                yield shell, "joint", numerator, sign, 0


def moment_row(plane, weights, targets, lower, upper):
    shell, kind, power, index, tick = plane
    start = 6 * shell
    bounds = targets[shell]
    rhs = (
        (bounds[index + 1] if kind == "tangent" else -bounds[index])
        if kind != "joint"
        else 0
    )
    coefficients = [0] * DIM
    for i, weight in enumerate(weights[shell], start):
        a, b = Fraction(lower[i], BOX_SCALE), Fraction(upper[i], BOX_SCALE)
        if kind == "tangent":
            x = a + Fraction(tick, 4) * (b - a)
            slope = power * x ** (power - 1)
            rhs += weight * (power - 1) * x**power
        elif kind == "secant":
            secant = sum(a**h * b ** (power - 1 - h) for h in range(power))
            slope = -secant
            rhs += weight * (a**power - a * secant)
        else:
            g = Fraction(power, 100)
            sign = index
            midpoint, radius = (a + b) / 2, (b - a) / 2
            derivative = 4 * midpoint**3 - 2 * g * midpoint
            intercept = midpoint**4 - g * midpoint**2 - derivative * midpoint
            curvature = min(0, 6 * a**2 - g) if sign == 1 else max(0, 6 * b**2 - g)
            slope = sign * derivative
            rhs -= weight * sign * (intercept + curvature * radius**2)
        coefficients[i] = floor(weight * slope)
    if kind == "joint":
        rhs += bounds[3] - g * bounds[0] if sign == 1 else g * bounds[1] - bounds[2]
    return coefficients, ceil(rhs)


def verify_tree(tree, fixed, weights, targets):
    planes = list(moment_planes())
    splits = leaves = contractions = 0

    def prove_lower(lower, upper, dual, objective=None):
        # 若 Ax≤b，则任意非负乘子 y 给出 c·x≥(c+yA)·x−yb。
        coefficients = [0] * DIM
        if objective is not None:
            i, sign = objective
            coefficients[i] = sign * DUAL_SCALE * ROW_SCALE
        numerator = 0
        for k, multiplier in dual:
            assert type(k) is int and 0 <= k < len(fixed) + len(planes)
            assert type(multiplier) is int and multiplier > 0
            row, rhs = (
                fixed[k]
                if k < len(fixed)
                else moment_row(planes[k - len(fixed)], weights, targets, lower, upper)
            )
            numerator -= multiplier * rhs * BOX_SCALE
            for i, coefficient in enumerate(row):
                coefficients[i] += multiplier * coefficient
        numerator += sum(
            c * (l if c >= 0 else h)
            for c, l, h in zip(coefficients, lower, upper, strict=True)
        )
        return Fraction(numerator, DUAL_SCALE * ROW_SCALE * BOX_SCALE)

    def visit(node, lower, upper):
        nonlocal splits, leaves, contractions
        if node[0] == "empty":
            assert len(node) == 2 and prove_lower(lower, upper, node[1]) > 0
            leaves += 1
        elif node[0] == "bound":
            _, changes, child = node
            for i, sign, dual in changes:
                assert type(i) is int and 0 <= i < DIM and sign in (-1, 1)
                bound = BOX_SCALE * prove_lower(lower, upper, dual, (i, sign))
                if sign == 1:
                    lower[i] = max(lower[i], floor(bound))
                else:
                    upper[i] = min(upper[i], ceil(-bound))
                assert lower[i] <= upper[i]
                contractions += 1
            visit(child, lower, upper)
        else:
            kind, i, mid, left, right = node
            assert kind == "split" and type(i) is int and 0 <= i < DIM
            assert type(mid) is int and lower[i] < mid < upper[i]
            left_upper, right_lower = upper.copy(), lower.copy()
            left_upper[i] = right_lower[i] = mid
            splits += 1
            visit(left, lower.copy(), left_upper)
            visit(right, right_lower, upper.copy())

    visit(tree, [0] * DIM, [BOX_SCALE] * DIM)
    assert leaves == splits + 1
    return splits, leaves, contractions


def exclude(root, model, moment_bounds):
    certificate = root / "ch39-six-cyclic-exclusion.json.gz"
    with gzip.open(certificate, "rt", encoding="utf-8") as data:
        cuts, tree = json.load(data)
    rows = spectral_rows(model, cuts)

    def below(x):
        return int((ROW_SCALE * x).lower().floor().fmpz())

    def above(x):
        return int((ROW_SCALE * x).upper().ceil().fmpz())

    weights = []
    for shell in model:
        row = [-shell["w"][i] for i in shell["big"]]
        values = [below(x) for x in row]
        assert len(values) == 6 and all(
            0 < w and above(x) <= w + 1 for w, x in zip(values, row, strict=True)
        )
        weights.append(values)
    # 向下取整每个权重误差至多一单位，每个读数在 [0,1]。六行合计至多六单位。
    targets = [
        [n for x in pair for n in (below(x) - 6, above(x))] for pair in moment_bounds
    ]
    counts = verify_tree(tree, rows, weights, targets)
    assert counts == (3600, 3601, 5223)
    print(
        "最后选线的联合谱界、二分次数、严格空叶及坐标收缩",
        len(cuts),
        *counts,
        flush=True,
    )


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    sieve = run_path(str(root / "ch39-six-cyclic-line-sieve.py"))
    exclude(
        root,
        *sieve["experiment"](
            root,
            *source["verify"](root, return_operator="spectral"),
            tolerance=arb("1e-6"),
        ),
    )

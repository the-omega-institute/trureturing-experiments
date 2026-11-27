"""六质量与最小混合角在完整参考窗内时，排除全部实系数二次源协方差。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def largest_component(moment, eigenvalues):
    """单位向量在最大本征线上的平方分量的两侧界。"""
    low, middle, high = eigenvalues
    lower = (moment - middle) / (high - middle)
    upper = (moment - low) / (high - low)
    return lower.lower().union(upper.upper())


def component_bounds(model, relative, eigenvalues):
    readings = [value + arb(0, width.upper()) for value, width in zip(model["ref"], model["tol"], strict=True)]
    up, down = readings[:3], readings[3:6]
    scale_squared = (prod(up) / (prod(down) * abs(relative.det())))**(arb(2) / 3)
    target_eigenvalues = [scale_squared * lam for lam in eigenvalues]
    # 用 T 的最大奇异值控制最重上型一行的前两项。
    deficit = target_eigenvalues[-1] * (down[0]**2 + down[1]**2) / up[2]**2
    assert 0 < deficit < arb(".0356")
    column = (readings[8].sin(), arb(0).union(deficit.sqrt().upper()),
              (1 - deficit).sqrt().lower().union(arb(1)))
    moment = sum((mass * modulus / down[2])**2 for mass, modulus in zip(up, column, strict=True))
    inverse_moment = sum((down[2] * modulus / mass)**2 for mass, modulus in zip(up, column, strict=True))
    high_weight = largest_component(moment, target_eigenvalues)
    low_weight = largest_component(inverse_moment, [1 / lam for lam in reversed(target_eigenvalues)])
    weights = (low_weight, 1 - low_weight - high_weight, high_weight)
    assert all(weight > 0 for weight in weights)
    low, middle, high = (v.sqrt() for v in weights)
    return (arb(1), high / middle, high / low), (down[1] / down[2])**2 / (low * middle)


def cover(rows, ratios, noise):
    errors = [noise * sum(abs(z) for z in row) for row in rows]
    cells = [(0, 0, 0)]
    leaves = deepest = 0
    while cells:
        depth, i, j = cells.pop()
        half = arb.pi() / 2**depth
        pad = arb(0, half.upper())
        x, y = (2 * i + 1) * half + pad, (2 * j + 1) * half + pad
        phases = [acb(0, angle).exp() for angle in (x-y, x, y)]
        def away_from_zero(row, error):
            value = sum(r * (z * phase).real for r, z, phase in zip(ratios, row, phases, strict=True))
            return abs(value) > error
        if any(away_from_zero(row, err) for row, err in zip(rows, errors, strict=True)):
            leaves += 1
            deepest = max(deepest, depth)
        else:
            assert depth < 12, (depth, i, j)
            cells.extend((depth+1, 2*i+di, 2*j+dj) for di in (0, 1) for dj in (0, 1))
    print("六质量和一个混合角的整窗排除，终格和最大深度", leaves, deepest)


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    quad = run_path(str(Path(__file__).with_name("ch39-two-zero-quadratic-source.py")))
    model, up, down = orbit["calibration"]()
    adj, trace = model["adj"], quad["trace"]
    grams = (adj(up) * up, adj(down) * down)
    normals, _ = quad["normal_frame"](*(g / trace(g).real for g in grams))
    relative = up * down.inv()
    eigenvalues, frame = model["eig"](adj(relative) * relative)
    generalized = down.inv() * frame
    axes = [adj(generalized) * normal * generalized for normal in normals]
    assert all(axis[i, i].contains(0) for axis in axes for i in range(3))
    ratios, noise = component_bounds(model, relative, eigenvalues)
    pairs = ((0, 1), (0, 2), (1, 2))
    rows = [[axis[i, j] for i, j in pairs] for axis in axes]
    cover(rows, ratios, noise)


if __name__ == "__main__":
    verify()

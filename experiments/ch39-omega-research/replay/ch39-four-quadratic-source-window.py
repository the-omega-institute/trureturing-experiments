"""在两个实际零点各前后六十四号内，检验四个冻结二次源的正定性。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def positive_covariance(words, coefficients):
    source = sum((word * c for word, c in zip(words, coefficients, strict=True)), acb_mat(3, 3))
    # 精确矩阵厄米。正的顺序主子式认证正定，任一负主子式排除半正定。
    first = source[0, 0].real
    if first < 0:
        return False
    second = first * source[1, 1].real - abs(source[0, 1])**2
    if second < 0:
        return False
    third = source.det().real
    if third < 0:
        return False
    assert first > 0 and second > 0 and third > 0, (first, second, third)
    return True


def verify():
    neighbours = run_path(str(Path(__file__).with_name("ch39-four-quadratic-source-neighbours.py")))
    load = neighbours["load"]
    quad = load("quadratic-gram-source")
    model, _, _ = load("two-zero-orbit-condition")["calibration"]()
    adj = model["adj"]
    trace = load("two-zero-quadratic-source")["trace"]
    centers = [quad["at_zero"](model, q, n) for q, n in ((-4, 2345666), (2, 964890))]
    branches = list(neighbours["coefficients_at_roots"](model, *centers))
    radius = 64
    offsets = range(-radius, radius + 1)
    sectors = []
    for center, (charge, index) in zip(centers, ((-4, 2345666), (2, 964890)), strict=True):
        grams = []
        for offset in offsets:
            site = quad["at_zero"](model, charge, index + offset) if offset else center
            grams.append(adj(site) * site)
        central_trace = trace(grams[radius]).real
        sectors.append((grams, [g * (central_trace / trace(g).real) for g in grams]))

    for convention in range(2):
        passes = [set() for _ in branches]
        for i, up in zip(offsets, sectors[0][convention], strict=True):
            for j, down in zip(offsets, sectors[1][convention], strict=True):
                if (i, j) == (0, 0):
                    # 每个中央协方差为可逆轨道矩阵与其伴随的乘积。
                    for passing in passes:
                        passing.add((0, 0))
                    continue
                words = neighbours["words"](up, down, quad)
                for passing, (_, coeff) in zip(passes, branches, strict=True):
                    if positive_covariance(words, coeff):
                        passing.add((i, j))
        upper_axes = [{i for i, j in passing if j == 0 and i} for passing in passes]
        lower_axes = [{j for i, j in passing if i == 0 and j} for passing in passes]
        assert all((0, 0) in passing for passing in passes)
        assert not any(lower_axes)
        expected_upper = ((3, 81, 120, 11), (3, 82, 124, 10))[convention]
        assert tuple(map(len, upper_axes)) == expected_upper
        first_failure = (7, 14)[convention]
        assert set(range(1, first_failure)) | set(range(-first_failure + 1, 0)) <= upper_axes[2]
        assert first_failure not in upper_axes[2]
        assert all(passing.isdisjoint({-1, 1}) for k, passing in enumerate(upper_axes) if k != 2)
        print("按原始或单位迹 Gram 的第", convention + 1, "种延拓")
        mixed_counts = ((0, 184, 687, 0), (0, 176, 840, 0))[convention]
        for k, passing in enumerate(passes):
            other = sorted((i, j) for i, j in passing if i and j)
            assert len(other) == mixed_counts[k]
            example = min(other, key=lambda p: (max(map(abs, p)), sum(map(abs, p)), p), default=None)
            print("枝、上型轴、下型轴、两部门同变的正定数、最近同变例",
                  k + 1, len(upper_axes[k]), len(lower_axes[k]), len(other), example)


if __name__ == "__main__":
    verify()

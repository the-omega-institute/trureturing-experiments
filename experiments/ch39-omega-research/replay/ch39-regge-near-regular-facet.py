"""Upper bound on the local Freudenthal Regge score near a maximum-size regular facet."""

from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


ctx.prec = 80
FACES = run_path(str(Path(__file__).with_name("ch39-regge-shape-support.py")))["FACES"]
HEIGHT_FLOOR = arb(113) / 1000
SCORE_LIMIT = arb(9) / 64
QUARTER = arb(1) / 4
BASE_AREA = arb(3).sqrt() / 4


def square(value):
    lower, upper = value.abs_lower(), value.abs_upper()
    return (lower * lower).union(upper * upper)


def score_and_gradient(edges, height2, barycentric):
    score = arb(0)
    gradient = [arb(0) for _ in edges]
    height = height2.sqrt()
    for i, j, a, b, c, valence in FACES:
        shift = arb.pi() * (arb(2) / valence - arb(1) / 2)
        if j == 4:
            cotangent = barycentric[i] / (3 * height2 / 2).sqrt()
            score += BASE_AREA * (shift + cotangent.atan())
            continue

        assert c == 4
        difference = edges[a] - edges[b]
        # 16 A^2 is at least 4 height^2 by orthogonal projection.
        area = (2 * (edges[a] + edges[b]) - 1 - square(difference)).max(4 * height2).sqrt() / 4
        li, lj = barycentric[i], barycentric[j]
        denominator = height * (2 * height2 + square(li) + square(lj) + square(li + lj) / 2).sqrt()
        deficit = shift + ((height2 / 2 - li * lj) / denominator).atan()
        score += area * deficit
        gradient[a] += deficit * (1 - difference) / (16 * area)
        gradient[b] += deficit * (1 + difference) / (16 * area)
    return score, gradient


def classify(box):
    centers = [arb((lo + hi) / 2) for lo, hi in box]
    radii = [arb((hi - lo) / 2) for lo, hi in box]
    intervals = [arb(c, r) for c, r in zip(centers, radii)]
    deviations = [c - sum(centers) / 4 for c in centers]
    lambdas = [QUARTER - d for d in deviations]
    hc = sum(centers) / 4 - arb(3) / 8 - sum(square(d) for d in deviations) / 2
    drift = sum(li.abs_upper() * r for li, r in zip(lambdas, radii))
    upper = hc + drift
    if upper < HEIGHT_FLOOR:
        return "empty"

    # The quadratic remainder of h^2 is -||delta - mean(delta)||^2/2.
    lower = hc - drift - sum(r * r for r in radii) / 2
    height2 = lower.union(upper).max(HEIGHT_FLOOR)
    ls = [arb(li, (sum(radii) + 2 * r) / 4) for li, r in zip(lambdas, radii)]
    score, gradient = score_and_gradient(intervals, height2, ls)
    if score < SCORE_LIMIT:
        return "direct"
    if hc >= HEIGHT_FLOOR:
        center_score, _ = score_and_gradient(centers, hc, lambdas)
        remainder = sum(g.abs_upper() * r for g, r in zip(gradient, radii))
        if center_score + remainder < SCORE_LIMIT:
            return "differential"
    return "split"


def main():
    queue = [(((0., 1.),) * 4, 0)]
    counts = Counter()
    while queue:
        box, depth = queue.pop()
        result = classify(box)
        counts[result] += 1
        if result == "split":
            # All endpoints and midpoints remain exact binary64 dyadics.
            assert depth < 52
            i = depth % 4
            lo, hi = box[i]
            mid = (lo + hi) / 2
            for interval in ((lo, mid), (mid, hi)):
                queue.append((box[:i] + (interval,) + box[i+1:], depth + 1))
    delta = arb(1) / 1000
    determinant_floor = arb(1) / 16 - 6 * delta
    upper = arb(96).sqrt() * (SCORE_LIMIT + 3 * arb.pi() * delta / determinant_floor.sqrt())
    assert upper < arb(178) / 100
    print("Proof leaves:", dict(counts))
    print("Near-facet upper bound:", upper)


if __name__ == "__main__":
    main()

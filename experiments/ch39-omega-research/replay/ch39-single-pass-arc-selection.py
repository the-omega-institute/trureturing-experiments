"""在每个高度的通过率不超过一时，联合读数容许的单次选择概率。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
entropy = run_path(str(root / "ch39-neutral-entropy-weights.py"))
ranges = entropy["ranges"]
scale = 2_903_600_000


def distinguish_heights(xs, first_acceptance):
    geometry = run_path(str(root / "ch39-arc-relative-determinant.py"))
    ctx.prec = 432
    clocks = [geometry["phases"](scale/x.sqrt()) for x in xs]
    flights = []
    for left, right in combinations(clocks, 2):
        angles = [abs((left[d]/right[d]).arg()) for d in (724, 804)]
        delta = max(angles, key=float)
        winding = int((arb.pi()/(2*delta)).ceil().unique_fmpz())
        cosine = (winding*delta).cos()
        assert cosine < 0
        one_side = 1/(2*(1-cosine))
        middle = -cosine/(1-cosine)
        assert 0 < one_side < 1 and 0 < middle < 1
        flights.append(winding)
    assert len(flights) == 55 and sum(flights) == 5454 and max(flights) == 1642
    pair_angle = max(abs((clocks[0][d]/clocks[1][d]).arg()) for d in (724, 804))
    required = first_acceptance.sqrt().acos()/pair_angle
    assert 242 < required < 243
    print("辨认十一枚的最少查询下界、构造查询数和原边数", max(flights), sum(flights), 3*sum(flights))
    print("仅实现最优放行概率所需的查询数下界", required.str(12))


def verify():
    rows, times, _, _, _, reference = entropy["joint"]["columns_at"](scale, ranges["indices"])
    ctx.prec = 640
    xs = [row[1] for row in rows]
    particle = ranges["indices"].index(143921293989490)
    candidates = []
    for support, subset, weights in ranges["feasible_vertices"](rows, times):
        ratio = [p*sum(xs)/xs[i] for p, i in zip(weights, support, strict=True)]
        stage = support.index(particle)
        assert all(ratio[stage] > v for j, v in enumerate(ratio) if j != stage)
        candidates.append((weights[stage], support, ranges["cp_loss"](subset, weights, reference)))
    assert len(candidates) == 10
    minimum, support, loss = min(candidates, key=lambda row: float(row[0]))
    assert support == tuple(range(8))
    assert all(minimum < p for p, other, _ in candidates if other != support)
    optimum = xs[particle]/(sum(xs)*minimum)
    assert arb("1.79981117188e-19") < optimum < arb("1.79981117189e-19")
    print("十一枚的最优单次通过率及该解的 CP 亏额", optimum.str(15), loss.str(14))

    original_weights = entropy["joint"]["solve"](rows[:8], times)
    first_acceptance = original_weights[0]*xs[particle]/(minimum*xs[0])
    distinguish_heights(xs, first_acceptance)
    ctx.prec = 176
    short = run_path(str(root / "ch39-short-arc-born-weights.py"))
    first_four = run_path(str(root / "ch39-first-window-required-weight.py"))["hits"]
    lower_z = sum(1/acb.zeta_zero(n).imag**2 for n in short["indices"])
    heights = [acb.zeta_zero(n).imag for n in first_four]
    upper_z = sum(1/h**2 for h in heights)+entropy["inverse_square_tail"](heights[-1].lower())
    lower = xs[particle]/(scale**2*upper_z*minimum)
    upper = entropy["inverse_square_tail"](15*10**12)/(lower_z*arb(".23893"))
    assert lower > arb("3.42589e-23") and upper < arb("4.42870e-5")
    print("全部三窗候选的最优通过率界", lower.str(15), upper.str(15))


if __name__ == "__main__":
    verify()

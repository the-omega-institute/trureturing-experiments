"""共用时间矩与味驻点时，实际零点能否同时保持两个下型质量矩。"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
scale_model = run_path(str(root / "ch39-physical-shared-scale.py"))
ensemble = run_path(str(root / "ch39-time-flavor-ensemble.py"))


def inspect_block(block):
    ctx.prec = 128
    i, first, end = block
    _, period, origins = ensemble["frequencies"]()
    width = arb("1.4e-7")
    reference_height = acb.zeta_zero(11661).imag
    exceptional = []
    for j in range(first, end):
        center = origins[i] + j*period
        if center < 13 or center > 1000001:
            continue
        band = center + arb(0, width)
        if abs(acb(arb(1)/2, band).zeta()) > 0:
            continue
        left = (center-2*width).zeta_nzeros().unique_fmpz()
        right = (center+2*width).zeta_nzeros().unique_fmpz()
        assert left is not None and right is not None
        if left != right:
            assert (left, right) == (11660, 11661) and band.contains(reference_height)
            exceptional.append((i, j))
    return exceptional


def verify_low_heights():
    ctx.prec = 192
    speed, period, origins = ensemble["frequencies"]()
    flavor = ensemble["flavor"]["flavor"]
    c0 = flavor["singlet"](11661)[1]
    axes = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    tolerance = arb("1e-9")
    slope = abs(speed)*(sum(v*v for v in axes)-(abs(c0)+tolerance)**2).sqrt()
    assert slope*arb("1.4e-7") > tolerance and acb.zeta_zero(1).imag > 14
    blocks = []
    for i, origin in enumerate(origins):
        first = int(((14-origin)/period).floor().unique_fmpz())
        last = int(((1000000-origin)/period).ceil().unique_fmpz())
        blocks.extend((i, n, min(last+1, n+512)) for n in range(first, last+1, 512))
    with ProcessPoolExecutor(max_workers=4) as pool:
        exceptional = [hit for group in pool.map(inspect_block, blocks) for hit in group]
    assert len(exceptional) == 1
    _, _, curvature = scale_model["physical"]["observables"](11661)
    assert curvature < arb(".98")
    print("百万高度内唯一命中编号与四层曲率比", 11661, curvature.str(12))


def verify():
    verify_low_heights()
    # 百万至四百万的排除由紧邻的已归档相位窗证书给出。
    ctx.prec, ctx.cap = 768, 3
    shared = scale_model["shared"]
    (e1, e2, _), _, lambdas, slopes = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    assert e1.degree() == 1 and e1[1] < 0 and e2.degree() == 2 and e2[2] < 0
    for limit, cutoff, probability_lower in (
            ("20000000", 17888543819, "8.36e-11"),
            ("23834887.62", 19528394760, "6.49e-20")):
        x = arb(limit)
        assert cutoff < 4000000*x.sqrt()
        score, (left, middle), terms = scale_model["upper_polynomial"](x, times, lambdas, slopes)
        # 在端点 X 的右侧有 (x-r)^2(x-s)^2/Q(x) <= 1/576.
        assert 0 < left < middle < x and score < 0
        excess = sum(d/z for d, z in zip(terms, lambdas, strict=True) if d > 0)/576
        forced = -score/excess
        assert forced > arb(probability_lower)
        drift = -e2[2]*forced*arb("1e-18")
        print("尺度上限、偏离参考相位的最少权重、质量矩组合亏额", cutoff, forced.str(12), drift.str(12))


if __name__ == "__main__":
    verify()

"""荷五传播子加 pI 后，八站相消分支在零动量处的边界导数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    endpoint = run_path(str(root / "ch39-proton-endpoint-speed.py"))
    ctx.prec = 192
    zero = acb.zeta_zero(177).imag
    area = endpoint["area"]
    cases = (
        (0, 5, "-.0459972646597", "-.0459972646594", "-.014004600481", "-.014004600477"),
        (2, 3, ".1064730230790", ".1064730230793", ".012489891297", ".012489891301"),
    )
    for r, k, lower, upper, limit_lo, limit_hi in cases:
        lo, hi = arb(lower), arb(upper)
        for _ in range(24):
            mid = (lo+hi)/2
            green = shell["laplacian"](5, zero+mid).inv()
            signed_area = area(green[6,r], green[6,k])
            if signed_area > 0:
                hi = mid
            else:
                assert signed_area < 0
                lo = mid
        green, dh = endpoint["green_and_slope"](shell, zero+lo.union(hi))
        dp = -green*green
        x, y = green[6,r], green[6,k]
        assert (x.conjugate()*y).real < 0
        h_slope = area(dh[6,r], y) + area(x, dh[6,k])
        p_slope = area(dp[6,r], y) + area(x, dp[6,k])
        assert h_slope > 0
        boundary_slope = -p_slope/h_slope
        assert arb(limit_lo) < boundary_slope < arb(limit_hi)
        print("边界站对、位移对 p 的局部导数", (r,k), boundary_slope.str(14))


if __name__ == "__main__":
    verify()

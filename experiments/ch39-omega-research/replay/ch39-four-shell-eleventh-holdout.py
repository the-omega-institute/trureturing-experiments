"""第七至十壳校准盒在第十一壳的留出读数和邻域排除。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def verify(root):
    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    here = Path(__file__).resolve().parent
    calibration = run_path(str(here / "ch39-four-shell-final-sign.py"))
    centered = run_path(str(here / "ch39-three-shell-final-sign.py"))["centered_residual"]
    fifth, columns, levels = load("fifth-full-positive-joint")["verify"](root, return_operator="spectral")
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    _, _, _, basis, _, weights, heat, field = geometry(
        root, fifth, columns, levels, number=2_306_888_100_096
    )
    ordered = load("two-shell-ordered-bond-gates")
    left = basis.conjugate().transpose()
    for edge in ordered["LAYERS"][-1][1]:
        for site in edge:
            for i in range(left.nrows()):
                left[i,site] = -left[i,site]
    station = fifth * columns
    angles, active = calibration["ANGLES"], calibration["ACTIVE"]
    center = [arb(angles[j]) for j in active]
    radii = list(map(arb, calibration["RADII"]))
    box = [c+arb(0,r.upper()) for c,r in zip(center,radii,strict=True)]
    moved, *_ = ordered["transport"](station, center, parameters=angles, active=active)
    residual = centered(left*moved, levels, weights, heat, field)
    projected, *derivatives = [left*m for m in ordered["transport"](station, box, parameters=angles, active=active)]
    _, jacobian = load("shortest-common-complex-gate")["gravity_system"](
        projected, derivatives, levels, weights, heat, field, False
    )
    error = [sum(z.abs_upper()*r for z,r in zip(row,radii,strict=True)) for row in jacobian.tolist()]
    newton = 1 + residual[0,0] + arb(0,error[0].upper())
    vacuum = arb("1e-122") + residual[1,0] + arb(0,error[1].upper())
    assert arb(".37799") < newton < arb(".37885")
    assert arb(".06957") < vacuum < arb(".06970")

    # A Hermitian perturbation of norm rho changes each unit-vector Rayleigh quotient by at most rho.
    rho = arb(".009")
    low, heavy, top = levels
    variation = arb(0)
    for i, weight in enumerate(weights):
        amplitude = low + (heavy-low)*abs(projected[i,0])**2 + (top-low)*abs(projected[i,1])**2
        lower = max(low.lower(), (amplitude.lower()-rho).lower())
        upper = min(top.upper(), (amplitude.upper()+rho).upper())
        assert 0 < lower <= upper
        exposure_sum = arb(0)
        for multiplicity, _, mass in field:
            assert mass >= 0
            exposure = (-mass/lower**2).exp().union((-mass/upper**2).exp())
            exposure_sum += multiplicity * exposure
        slope = weight * lower.union(upper)**3 * exposure_sum / (2*arb.pi())
        variation += rho * slope.abs_upper()
    assert variation < arb(".0631")
    assert vacuum.lower() - variation.upper() > arb("1e-122")
    print("留出牛顿、真空读数及 0.009 邻域中的真空变化上限",
          newton.str(12), vacuum.str(12), variation.str(12))


if __name__ == "__main__":
    ctx.prec = 192
    verify(Path(sys.argv[1]).resolve())

"""共享三循环等变正势时，第 11431 个零点所需的迹预算。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

winding = run_path(str(Path(__file__).with_name("ch39-common-loop-winding.py")))
point = winding["spectrum"]
study, native = point["study"], point["native"]


def verify():
    ctx.prec = 192
    height = acb.zeta_zero(11431).imag
    gain = arb(2)**-10
    floor = arb("1.3")-4*arb.pi()/32
    a, b = winding["counts"]
    ray = acb(3, 13*arb(3).sqrt())
    phase = (ray/abs(ray))**-493
    vector = acb_mat([[(arb(a)/(a+b)).sqrt()], [phase*(arb(b)/(a+b)).sqrt()]])
    central_scale = gain**2/study["center"]

    # A positive source of bounded trace can only lower a competitor's energy.
    bare = point["negative_energy"](height, vector, gain, gauge_weight=0)
    competitor = acb.zeta_zero(11139064).imag
    other = point["negative_energy"](competitor, vector, gain, gauge_weight=0)
    advantage = (other-bare)/central_scale
    assert advantage > arb(".87218")
    budget = arb("2.60")
    mixing = (gain*arb(2).sqrt()+central_scale*budget)/(floor-central_scale*budget)
    center_drop = budget*(1/arb(3).sqrt()+mixing)**2
    assert advantage > center_drop
    print("Bare advantage, maximum allowed central decrease at trace", budget,
          advantage.str(13), center_drop.str(13))

    prove_orbit_gap(height, vector, gain, arb(23)/25, arb(1)/50000, arb("7.44e-9"),
                    arb("3.02e-17"))
    prove_orbit_gap(height, vector, gain/arb(3).sqrt(), arb(1), arb(1)/10000,
                    arb("7.8e-9"), arb("1e-16"))


def prove_orbit_gap(height, vector, gain, weight, curvature, allowance, required):
    floor = arb("1.3")-4*arb.pi()/32
    # Start with the marked point only. Add the other two orbit sites below.
    gap = point["full_spectrum"](height, vector, gain, curvature=curvature,
                                 allowance=allowance, gauge_weight=weight)
    shift = point["negative_energy"](height, vector, gain, gauge_weight=weight)
    cutoff = 512
    high_scale = weight*gain**2/(cutoff**2+study["offset"])
    correction = high_scale*2*gain**2/(floor-high_scale)**2
    assert gap-correction > required

    # For smaller h compare each new Schur complement to the old center root.
    count = int(arb(cutoff).zeta_nzeros().unique_fmpz())
    remaining = shift
    indices = [native["inside"].index(i) for i in (3, 4)]
    for zero in acb.zeta_zeros(1, count):
        h = zero.imag
        D, J = point["internal"](h)
        source = J*vector
        scale = weight*gain**2/(h*h+study["offset"])
        for i in indices:
            D[i, i] -= scale
        schur = shift-scale-gain**2*(point["adjoint"](source)*point["solve"](D,source,shift))[0, 0].real
        slope = 1+2*gain**2/(floor-scale)**2
        remaining = min(remaining.lower(), (schur/slope).lower())
    assert remaining > required
    print("Per-site weight, flavor coupling, low count and gap, high-zero corrected gap",
          weight.str(7), gain.str(9),
          count, remaining.str(12), (gap-correction).str(12))


if __name__ == "__main__":
    verify()

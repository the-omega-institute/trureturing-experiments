"""三条原生辐边的共同权重、常量零模和破缺矢量阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
charged = run_path(str(root / "ch39-clocked-pole-current-defect.py"))
winding = run_path(str(root / "ch39-common-loop-winding.py"))
point = winding["spectrum"]
geo, study = point["geo"], point["study"]
orbit, pole = charged["orbits"][0], charged["pole"]
gain = arb(2)**-10


def flavor_port():
    a, b = winding["counts"]
    ray = acb(3, 13*arb(3).sqrt())
    phase = (ray/abs(ray))**-493
    return acb_mat([[(arb(a)/(a+b)).sqrt()], [phase*(arb(b)/(a+b)).sqrt()]])


def flavor_checks(h, weight, max_scale, center_scale):
    native_edges = {(i, j) for i, j, _, _ in geo["edges"]}
    assert all((i, j) not in native_edges for i in orbit for j in orbit if i < j)
    floor = arb("1.3")-4*arb.pi()/32
    port = flavor_port()
    bare = point["negative_energy"](h, port, gain, gauge_weight=0)
    other = point["negative_energy"](acb.zeta_zero(11139064).imag, port, gain, gauge_weight=0)
    advantage = (other-bare)/center_scale
    budget = arb("5.19")
    mixing = (gain*arb(2).sqrt()+center_scale*budget)/(floor-center_scale*budget)
    assert advantage-budget*(1/arb(6).sqrt()+mixing)**2 > arb(".00072")

    # Variational comparison of the star against a single marked site,
    # valid on the actual ground vectors of either matrix.
    mixing = (gain*arb(2).sqrt()+weight*max_scale)/(floor-4*weight*max_scale)
    error = 2*mixing+4*mixing**2
    upper = arb(".924")
    assert weight*(1+error) < upper
    gap = point["full_spectrum"](h, port, gain, curvature=arb(1)/50000,
                                 allowance=arb("7.44e-9"), gauge_weight=upper)
    loss = (upper-weight*(1-error))*center_scale
    assert gap-loss > arb("1.5e-17")
    print("Native edge trace, comparison error, full flavor gap lower",
          (6*weight).str(8), error.str(9), (gap-loss).str(13))


def vector_threshold(h, weight, max_scale, center_scale):
    L = charged["charge_blocks"](h)
    neutral = acb_mat(19, 19)
    for i, j, _, _ in geo["edges"]:
        neutral[i, i] += 1
        neutral[j, j] += 1
        neutral[i, j] = neutral[j, i] = -1
    neutral_lift = neutral+acb_mat([[arb(1)/19]*19 for _ in range(19)])
    spokes = acb_mat([[int(i == v)-int(i == pole) for v in orbit] for i in range(19)])
    star = spokes*point["adjoint"](spokes)
    assert star.trace().real == 6
    eye = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
    changes = []
    for block in (L, neutral_lift):
        gram = point["adjoint"](spokes)*block.solve(spokes)
        ratio = (eye-weight*center_scale*gram).det().real
        direct = (block-weight*center_scale*star).det().real/block.det().real
        assert abs(ratio-direct) < arb("1e-43")
        changes.append(ratio.log())
    shift = -7*(changes[0]-changes[1])/(8*arb.pi())
    bound = -21*(1-weight*max_scale).log()/(8*arb.pi())
    print("Static alpha2^-1-alpha3^-1 shift, every positive zero/cutoff upper",
          shift.str(16), bound.str(13))
    assert arb("5.92617673619e-18") < shift < arb("5.92617673620e-18")
    assert bound < arb("3.186e-9")


def verify():
    ctx.prec = 224
    charged["geometry"]()
    h = acb.zeta_zero(11431).imag
    scale = lambda height: gain**2/(height**2+study["offset"])
    weight = arb(23)/25
    max_scale = scale(acb.zeta_zero(1).imag)
    center_scale = scale(h)
    flavor_checks(h, weight, max_scale, center_scale)
    vector_threshold(h, weight, max_scale, center_scale)


if __name__ == "__main__":
    verify()

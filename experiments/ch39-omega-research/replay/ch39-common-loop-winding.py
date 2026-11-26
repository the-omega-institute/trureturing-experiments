"""固定重子闭路，同一闭路在味端选中 11431 所需的最小绝对绕数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

root = Path(__file__).parent
joint = run_path(str(root / "ch39-joint-path-minimum.py"))
direction = run_path(str(root / "ch39-pure-port-direction-window.py"))
study, spectrum = direction["study"], joint["spectrum"]
counts = joint["a"], joint["b"]
half_loops = (
    ((14, 0, 10), 14), ((14, 13, 0), -6), ((14, 7, 6), -21),
    ((14, 6, 9), 6), ((14, 6, 15), 1), ((14, 9, 13), -1),
    ((14, 18, 9), 8),
)


def response(height):
    angles = study["geo"]["phases"](height)
    return direction["shifted_response"](angles[724], angles[804])


def exclude_shorter(units, limit, reference):
    a, b = counts
    imbalance = arb(a-b)/(a+b)
    transverse = 2*arb(a*b).sqrt()/(a+b)
    phases = {(i, k): unit**k for i, unit in enumerate(units)
              for k in range(1-limit, limit)}
    center, derivative = response(reference)
    zmax = arb(2)**-20
    assert 1/study["center"]-2*center[0] < 1
    # The center's negative root has magnitude at most g^2*(inv_* - tr(T_*)).
    # The resolvent series through first order has per-block error <= 2*z^2/d^3.
    shift = arb(0).union(zmax)
    tail = 4*zmax**2/(arb("1.3")-4*arb.pi()/32)**3
    for n in (4955890, 7345613):
        height = acb.zeta_zero(n).imag
        values, slopes = response(height)
        diff = [v-w for v, w in zip(values, center)]
        ds = [v-w for v, w in zip(slopes, derivative)]
        gauge_difference = 1/study["center"]-1/(height**2+study["offset"])

        def evaluate(row, phase):
            return (row[0]+imbalance*row[1]
                    +transverse*(row[2]*phase.real-row[3]*phase.imag))
        def beaten(phase):
            upper = evaluate(diff, phase)+shift*evaluate(ds, phase)+tail+gauge_difference
            return upper < 0

        phases = {k: ph for k, ph in phases.items() if not beaten(ph)}

    assert not phases, list(phases)
    print("Two real zero indices exclude all signed windings of magnitude below", limit)


def certify_port(phase, height):
    a, b = counts
    port = acb_mat([[(arb(a)/(a+b)).sqrt()], [phase*(arb(b)/(a+b)).sqrt()]])
    gap = spectrum["full_spectrum"](height, port, arb(2)**-10,
                                    curvature=arb(1)/10000, allowance=arb("7.8e-9"))
    assert gap > arb("3.14e-16")


def verify():
    ctx.prec = 192
    unit = acb(-83, 13*arb(3).sqrt())/86
    height = acb.zeta_zero(11431).imag
    assert joint["word"](joint["baryon_loops"]) == (-35, 13)
    exclude_shorter([unit], 3163, height)
    phase = unit**3163
    assert arb(".89744863760") < phase.arg() < arb(".89744863761")
    certify_port(phase, height)

    ray = joint["word"](half_loops)
    assert ray == (8, 13) and joint["multiply"](ray, ray) == (-35, 13)
    half = acb(3, 13*arb(3).sqrt())
    half /= abs(half)
    exclude_shorter([half, -half], 493, height)
    phase = half**-493
    assert arb(".89746764319") < phase.arg() < arb(".89746764320")
    certify_port(phase, height)


if __name__ == "__main__":
    verify()

"""Read a floating-window weight from a calibrated three-hop return and bound its precision."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
gauge = run_path(str(root / "ch39-floating-scale-winner.py"))
local = run_path(str(root / "ch39-three-step-zero-readout.py"))
q = gauge["q"]


def competitor_witness():
    # Identity quark mixing is an invariant submodel of the Yukawa flow.
    time = q("26.859794675")
    common = gauge["target"](10871, time)
    logs = list(map(q, ("-25.4446915838 -12.9928374813 -1.5528331663 -23.8189018241 "
                        "-17.8423094899 -10.1208976879 -25.6505068564 -14.9499484326 "
                        "-9.2850641599").split()))
    start, y, _ = gauge["evolve"]([common*r for r in gauge["ratios"]], logs, -time, 131072,
                                  diagonal_quarks=True)
    factors = [1/(x.sqrt()*c) for x, c in zip(start, gauge["centers"])]
    assert all(q(".999") < x < q("1.001") for x in factors)
    assert all(q(".99") < (v/2).exp()/c < q("1.01") for v, c in zip(y, gauge["yukawa_centers"]))
    print("10871 feasible, initial weak and strong ratios", [x.str(9) for x in factors[1:]])


def calibrated_return():
    nodes, edges = local["shell"]["geometry"]()
    graph = [{} for _ in nodes]
    for i, j, amplitude, charge in edges:
        graph[i][j] = amplitude, charge
        graph[j][i] = amplitude, -charge
    index = nodes.index((33, 33, 75))
    constant, slope, counts = local["triangle_terms"](graph, index)
    assert counts == {0: 12, 1: 4} and slope == Fraction(3517425, 42342346661888)
    reference_plus, reference_minus = q(constant + slope), q(constant - slope)
    for n in (10871, 10875):
        rho = acb.zeta_zero(n)
        matrix = local["shell"]["matrices"](nodes, edges, (rho-1)/rho)
        observed = (matrix**3)[index, index]
        assert observed.imag.contains(0)
        p = (reference_plus - observed.real)/(reference_plus - reference_minus)
        assert p.overlaps(1/(4*rho.imag**2+1)) and p.rad() < arb("1e-30")
        time = q("26.859482438")
        fixed = gauge["fixed"]
        prefactor = (fixed["geo"]["BOOK_PLANCK"]/1000)**2 * 3 * arb(5).sqrt() / (4*arb.pi())
        predicted = prefactor*(-2*time).exp() * 4*p / (1+4*p*(q(fixed["OFFSET"])-q("1/4")))
        assert predicted.overlaps(gauge["target"](n, time)) and predicted.rad() < arb("1e-19")
    radius = Fraction(1, 2**48)
    # Three raw readings with independent absolute error e imply |delta p| <= e/(slope-e).
    assert Fraction(29, 10**23) < slope*radius/(1+radius)
    print("Sufficient raw return error per setting: 2.9e-22")


def rank():
    root_two, log_two = arb(2).sqrt(), arb(2).log()
    first, last, winner = 10855, 10964, 10875
    p = [1/(4*z.imag**2+1) for z in acb.zeta_zeros(first, last-first+2)]
    radius = arb(2)**-48
    # The absolute four-tap modulus error follows from the 64-degree coefficient certificate.
    error = q(".000092")
    assert q(".00001264")*(3+3*root_two) < error

    def modulus(u):
        return (3-2*u)*(2*(113+72*u)).sqrt()/13

    def weight(x, uncertainty):
        u = root_two * (((1/x-1).sqrt()/2)*log_two).cos()
        return 4*x*(modulus(u)+uncertainty)

    jitter = arb(0, error)
    partition = [arb(2*j+1-64, 1)*radius/64 for j in range(64)]
    weights = []
    for value in p[:-1]:
        parts = [weight(value+delta, jitter) for delta in partition]
        weights.append((min(x.lower() for x in parts), max(x.upper() for x in parts)))
    floor = weights[winner-first][0]
    rivals = [hi for n, (_, hi) in enumerate(weights, first) if n != winner]
    # The global modulus maximum has u=-59/108, independent of the instrument error.
    tail = 4*(p[-1]+radius)*(modulus(-arb(59)/108)+error)
    assert floor-max(tail, max(rivals)) > q("2.5e-12")
    # A strictly positive input is guaranteed throughout the necessary feasible band.
    assert 1/(4*acb.zeta_zero(12092).imag**2+1) > radius

    double_radius = 2*radius
    reversed_gap = weight(p[10871-first]+double_radius, 0) - weight(p[winner-first]-double_radius, 0)
    assert reversed_gap > q("2.1e-12")
    print("2^-48 worst-case winning gap", (floor-max(tail, max(rivals))).str(10),
          "2^-47 reversal gap", reversed_gap.str(10))


def controlled_triangle():
    # Two coherent paths traverse the same charged triangle in opposite senses.
    theta_min = 2*(1/(2*arb(11644))).atan()
    theta_max = 2*(1/(2*arb(10570))).atan()
    levels = range(17, 32)
    trials = 2048
    angle_error = (arb(2).sqrt()/8).asin()
    assert angle_error < q(".18")
    assert 2**17*(theta_max-theta_min) + angle_error < arb.pi()
    assert 3*angle_error < arb.pi()
    theta_error = angle_error / 2**32
    # Project the final angle onto the prior interval before converting to p.
    assert theta_max.sin()*theta_error/2 < arb(2)**-48
    candidates = 12092-10855+1
    failure = 4*len(levels)*candidates*(-arb(trials)/128).exp()
    assert failure < q(".0084")

    # Include both endpoints, both sides of the ranking counterexample, and all
    # four extremal combinations of quadrature errors in deterministic replays.
    angles = [theta_min, theta_max]
    for n in (10871, 10875):
        angles.append(2*(1/(2*acb.zeta_zero(n).imag)).atan())
    for theta in angles:
        for re_sign, im_sign in ((1,1), (1,-1), (-1,1), (-1,-1)):
            estimate = (theta_min+theta_max)/2
            for j in levels:
                factor = 2**(j+1)
                phase = factor*theta
                measured = acb(phase.cos()+q(re_sign)/8, phase.sin()+q(im_sign)/8)
                principal = measured.arg()
                turns = ((factor*estimate-principal)/(2*arb.pi())+q("1/2")).floor().unique_fmpz()
                assert turns is not None
                estimate = (principal + 2*arb.pi()*int(turns))/factor
                assert abs(estimate-theta) < q(".18")/factor
    hops = 3 * 2 * trials * sum(2**j for j in levels)
    print("Controlled triangle per mode:", 2*len(levels)*trials, "trials,", hops,
          "edge traversals, simultaneous failure <", failure.upper().str(6))


if __name__ == "__main__":
    ctx.prec = 160
    calibrated_return()
    rank()
    controlled_triangle()
    competitor_witness()

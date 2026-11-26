"""Recover a zeta phase from a projected shell response and separate six candidates."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx, fmpq, fmpq_mat

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = shell["ball"]
cut_denominators = (22271, 22273, 22274, 22275, 22277)
edge_weight = F(12, 56)**3 / 19


def inverse_port_entries(matrix, clock, source, corner):
    pencil = -matrix
    for i in range(matrix.nrows()):
        pencil[i, i] += clock
    resolvent = pencil.inv()
    observed = acb_mat([[resolvent[i, j] for j in range(19)] for i in range(19)])
    inverse = observed.inv()
    return inverse[source, corner], inverse[corner, source]


def relay_reference(nodes, edges, source, corner):
    count = len(nodes)-19
    pencil = fmpq_mat([[int(i == j) for j in range(count)] for i in range(count)])
    incident = fmpq_mat(2, count)
    for i, j, weight, exponent in edges:
        if j < 19:
            continue
        assert exponent == 0
        a = fmpq(weight.numerator, weight.denominator)
        if i >= 19:
            pencil[i-19, j-19] = pencil[j-19, i-19] = -a
        elif i in (source, corner):
            incident[int(i == corner), j-19] = a
    t = (incident*pencil.solve(incident.transpose()))[0, 1]
    return F(int(t.p), int(t.q))


def verify():
    ctx.prec = 192
    nodes, edges = shell["geometry"]()
    source = nodes.index((21, 51, 69))
    corner = nodes.index((33, 33, 75))
    assert (source, corner, edge_weight, 1) in edges
    assert all(exponent == 0 for _, j, _, exponent in edges if j >= 19)
    clocks = (acb(1), acb(-1, arb(3).sqrt())/2)
    cuts = [ball(F(2*q, 1+q*q)) for q in cut_denominators]
    relay = relay_reference(nodes, edges, source, corner)
    assert F("6.66418e-6") < relay < F("6.66419e-6")
    intensity_cuts = [
        ball(edge_weight**2 + relay**2 + 2*edge_weight*relay*F(q*q-1, q*q+1))
        for q in cut_denominators
    ]
    margins, readings, heights, intensity_margins = [], [], [], []
    for n in range(11485, 11491):
        rho = acb.zeta_zero(n)
        phase = (rho-1)/rho
        matrix = shell["matrices"](nodes, edges, phase)
        for clock in clocks:
            forward, backward = inverse_port_entries(matrix, clock, source, corner)
            observed = (backward-forward)/acb(0, 2*ball(edge_weight))
            assert observed.imag.contains(0) and observed.imag.rad() < arb("1e-40")
            s = observed.real
            assert s.overlaps(phase.imag) and s.rad() < arb("1e-40")
            height = (1+(1-s*s).sqrt())/(2*s)
            assert height.overlaps(rho.imag) and height.rad() < arb("1e-34")
            assert all(abs(s-cut) > arb("1.6057e-9") for cut in cuts)
            assert 11485 + sum(s < cut for cut in cuts) == n
            if clock == 1:
                intensity = forward.real**2 + forward.imag**2
                x = (intensity-ball(edge_weight**2+relay**2))/ball(2*edge_weight*relay)
                height_squared = (1+x)/(4*(1-x))
                assert height_squared.overlaps(rho.imag**2) and height_squared.rad() < arb("1e-28")
                assert all(abs(intensity-cut) > arb("9.9e-22") for cut in intensity_cuts)
                assert 11485 + sum(intensity > cut for cut in intensity_cuts) == n
                intensity_margins.extend(abs(intensity-cut) for cut in intensity_cuts)
        heights.append(2*height)
        readings.append(s)
        margins.extend(abs(s-cut) for cut in cuts)
        print(n, "twice recovered height", (2*height).str(20), "phase ordinate", s.str(17))

    for left, q, right in zip(heights, cut_denominators, heights[1:]):
        assert left < q < right
    assert 22276 < heights[-2] < 22277 < heights[-1] < 22278
    gaps = [u-v for u, v in zip(readings, readings[1:])]
    assert all(g > gaps[-1] for g in gaps[:-1])
    assert arb("3.72975e-9") < gaps[-1] < arb("3.72976e-9")
    assert abs(readings[-1]-cuts[-1]) < arb("1.6058e-9")
    print("minimum rational-cut margin", min(margins).str(14), "closest phase gap", gaps[-1].str(14))
    assert len(intensity_margins) == 30
    print("minimum squared-modulus margin", min(intensity_margins).str(14))

    # Transfer an operator-norm error on the entire projected matrix to the scalar.
    norm_bound = F(98, 1000)
    inverse_bound = 1 + norm_bound + norm_bound**2/(1-norm_bound)
    error = F(6, 10**13)
    propagated = inverse_bound**2 * error / ((1-inverse_bound*error)*edge_weight)
    assert inverse_bound*error < 1 and propagated < F("1.6057e-9")


if __name__ == "__main__":
    verify()

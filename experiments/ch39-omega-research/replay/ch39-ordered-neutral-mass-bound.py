"""用全中性谱的第二个正质量与两处族流排除动能正定的局部矢量。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb_poly

budget = run_path(str(Path(__file__).with_name("ch39-positive-kinetic-spectral-budget.py")))
poles, mat, q = (budget[key] for key in ("poles", "mat", "q"))


def variations_in(polynomial, left, right):
    shifted = polynomial(arb_poly([left, 1]))
    degree = shifted.degree()
    scaled = [shifted[i]*(right-left)**i for i in range(degree+1)]
    sequence = [sum(scaled[i]*comb(degree-i, j-i) for i in range(j+1))
                for j in range(degree+1)]
    assert all(not v.contains(0) for v in sequence)
    signs = [v > 0 for v in sequence]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def count_old_poles(ports, z, b, d, cut):
    # The cleared determinant acts on the current-coupled tree realization.
    unit = mat["unit"]
    (pc, vc), (pr, vr) = ([arb_poly(a) for a in pair] for pair in ports)
    ac = vc*z-pc*(9/(24*unit))
    ar = -pr*(1/(4*unit))
    polynomial = (ac*vr+ar*vc)*arb_poly([b, -d*unit])+4*vc*vr
    assert polynomial[0] == 0 and not polynomial[1].contains(0)
    positive = arb_poly(polynomial.coeffs()[1:])
    assert variations_in(positive, q(0), cut/unit) == 1


def check_dark_masses(ports, hyper_ratio, triplet_squared, cut):
    (pc, vc), (pr, vr) = ([arb_poly(a) for a in pair] for pair in ports)
    det = arb_poly([0, 1])*pc*pr+hyper_ratio*(vc*pr*q("2/3")+pc*vr)
    assert det[0] == 0
    restricted = arb_poly(det.coeffs()[1:])
    endpoints = [q(v)**2*triplet_squared/mat["unit"] for v in
                 ("2", "2.087", "2.091", "2.432", "2.437")]
    endpoints.append(cut/mat["unit"])
    for left, right, expected in zip(endpoints[:-1], endpoints[1:], (0, 1, 0, 1, 0), strict=True):
        assert right > left and variations_in(restricted, left, right) == expected


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"]
    targets = [target(n, time) for n in range(11423, 11431)]
    *_, colour, right = poles["neutral"]["graph_constants"](port_builder=poles["neutral"]["pole_polynomials"])
    for case_index, logs in enumerate(poles["retune"]["log_sets"]):
        case = budget["reference"](logs, time)
        (scale, z, rho, eta, c, _), _, _ = poles["retune"]["predictions"](logs, time)
        b, d = c/rho, eta/rho
        e = poles["reference_readings"](logs)[3]**2/scale
        e2 = mat["geometry"]["second_neutral_mass"](scale, z, b, d)**2/scale
        relative_loss = q("0.00001")
        cut = e2*(1-relative_loss)**2
        assert 0 < e < cut < e2
        count_old_poles((colour, right), z, b, d, cut)
        fraction = 1/(4*case.mass_ratio_squared+q(1).exp())
        pole_tree = mat["massive_exchange"](-cut, z, b, d)
        for label, goal in enumerate(targets, 11423):
            metric, matrices = budget["spectral_forms"](case, goal, budget["tolerance"])
            (first, second), _ = budget["at_mass"](case, matrices, q(2))
            k = 1/(poles["offset"]+goal/(15*scale))
            check_dark_masses((colour, right), case.hypercharge/k-q("1/24"),
                              case.mass_ratio_squared*e, cut)
            pole_form = -(pole_tree*k+mat["hypercharge"]*case.hypercharge)*fraction
            # Necessary upper bounds on both the second exchange and the
            # number of positive masses below the chosen cut.
            good = any(mat["positive"](metric-first+(second-first)*z2+(pole_form-first)*zm)
                       for z2 in (6, 8, 10, 12, 16, 1, 3, 20, 30)
                       for zm in (q(".0001"),q(".001"),q(".01"),q(".1"),q(1),q(3),q(10)))
            if not good:
                raise ArithmeticError(f"no dual witness for {case_index=} {label=}")
        print("case",case_index,"all eight labels force a second mass below",
              (1-relative_loss).str(10),"times reference; 2M3/M2",
              (4*case.mass_ratio_squared*e/e2).sqrt().str(8),flush=True)


if __name__ == "__main__":
    verify()

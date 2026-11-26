"""Classify the absolute normalization left free by tree-level resonance masses."""
from pathlib import Path
from runpy import run_path

from flint import arb

poles = run_path(str(Path(__file__).with_name("ch39-charged-pole-zero-index.py")))
q, unit, roots = (poles[key] for key in ("q", "unit", "roots"))
green = poles["cross"]["fit"]["green"]
lepton_width = poles["cross"]["fit"]["widths"]["normalized_lepton_width"]
physical_masses = poles["four"]["physical_masses"]


def second_neutral_mass(scale_squared, z, b, d):
    # Use the reduced integer polynomial independently of the inverse's graph evaluation.
    _, _, (pc, vc), (pr, vr) = poles["neutral"]["graph_constants"](
        port_builder=poles["neutral"]["pole_polynomials"])
    num, den = 3*pc*vr + 2*pr*vc, 8*vc*vr
    common = num.gcd(den)
    num, den = num // common, den // common
    bare = []
    for value, multiplicity in den.complex_roots():
        assert multiplicity == 1 and value.imag == 0
        bare.append(value.real * unit)
    bare.sort(key=float)
    lo, hi = q(".10885"), q(".11")
    assert bare[0] == 0 and q(".003") < bare[1] < lo < hi < min(bare[2], b/d)
    assert all(a < b_next for a, b_next in zip(bare, bare[1:]))

    def secular(energy):
        x = energy / unit
        return z - num(x)/(unit*den(x)) + 4/(b-d*energy)

    return (scale_squared * roots(secular, lo, hi)).sqrt()


def neutral_shape(scale_squared, z, neutral_masses):
    first, second = (m**2/scale_squared for m in neutral_masses)
    delta1, delta2 = (-4/(z+green(e)) for e in (first, second))
    d = (delta1-delta2)/(second-first)
    b = delta1+d*first
    assert min(d, b) > 0
    return b, d


def check_neighbours(targets, model, neutral_mass, width):
    scale_squared, z, rho, eta, c, masses = model
    b, d = c/rho, eta/rho
    for target in targets:
        alternative = poles["offset"] + target/(15*scale_squared)
        assert alternative > 0 and not alternative.overlaps(rho)
        scaled_masses = physical_masses(b*alternative, alternative, z)
        assert all(x.overlaps(y) for x, y in zip(masses, scaled_masses, strict=True))
        scaled_width = lepton_width(alternative, d*alternative, b*alternative, z)
        assert (scaled_width/width).overlaps(rho/alternative)
        # Retain the observed lepton width but drop the second neutral mass.
        eta_new, c_new, _ = poles["recover_neutral"](neutral_mass, width, scale_squared, z, alternative)
        assert lepton_width(alternative, eta_new, c_new, z).overlaps(width)
        ratio = (alternative/eta_new).sqrt()
        assert abs(ratio-1) > q(".005")
        print("neighbour width shift (ppm)", (10**6*(rho/alternative-1)).str(10),
              "coupling ratio after dropping the second neutral mass", ratio.str(12))


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"]
    neighbours = [target(n, time) for n in (11430, 11431, 11432)]
    offset = poles["offset"]
    assert offset > 0
    poles["verify_ratio_monotonicity"]()
    graph_gap = poles["window"]["graph_gap"]()
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        model, _, width = poles["retune"]["predictions"](logs, time)
        scale_squared, z, rho, eta, c, _ = model
        charged1, charged2, _, neutral1, _ = poles["reference_readings"](logs)
        neutral2 = second_neutral_mass(scale_squared, z, c/rho, eta/rho)
        print("case", case, "second neutral M/U", neutral2.str(14))
        check_neighbours((neighbours[0], neighbours[2]), model, neutral1, width)

        mass_error = arb(0, q("1e-8").upper())
        box = [m*(1+mass_error) for m in (charged1, charged2, neutral1, neutral2)]
        fitted_scale, fitted_z, _ = poles["charged_geometry"](*box[:2])
        b, d = neutral_shape(fitted_scale, fitted_z, box[2:])
        assert abs(1/d.sqrt()-1) < q(".00001824")
        observed_width = width*(1+arb(0, q("1e-5").upper()))
        fitted_rho = lepton_width(arb(1), d, b, fitted_z)/observed_width
        answer = 15*fitted_scale*(fitted_rho-offset)
        print("arithmetic interval", answer.lower().str(15), answer.upper().str(15), flush=True)
        assert neighbours[0] > answer > neighbours[2] and answer.contains(neighbours[1])
        assert q("3.29402227") < answer < q("3.29413040")
        fitted_rho = fitted_rho.intersection(neighbours[1]/(15*fitted_scale)+offset)
        fitted_masses = physical_masses(b, arb(1), fitted_z)
        fitted = (fitted_scale, fitted_z, fitted_rho, d*fitted_rho, b*fitted_rho, fitted_masses)
        windows = poles["window"]["window_data"](fitted, graph_gap)
        poles["cross"]["exclude_other_orbits"](case, windows, q(logs[2]))


if __name__ == "__main__":
    verify()

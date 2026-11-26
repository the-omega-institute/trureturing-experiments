"""Test leading final-state QCD radiation in the three retuned charged-current witnesses."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
poles = run_path(str(root / "ch39-charged-pole-zero-index.py"))
retune, window, study, q = (poles[k] for k in ("retune", "window", "study", "q"))
comp, adj = study["comp"], study["adjoint"]


def gluon_radiator(index, windows):
    onset, h3, end, _, first, second = windows
    endpoints = []
    assert adj["ratios"][2] == 1
    for upper in (False, True):
        state = window["flavour"]["before"](index, upper)
        side = 2 if upper else 3
        sm = adj["flow"](False, False)
        segments = (
            (q(23), onset, sm),
            (onset, h3, {**sm, "NEUTRAL_YUKAWA": first[side]}),
            (h3, end, {**adj["flow"](True, False), "NEUTRAL_YUKAWA": second[side]}),
        )
        for a, b, flow in segments:
            assert b > a
            for _ in range(16):
                state = comp["advance"](state, (b - a) / 16, upper, flow=flow)
        inverse_colour = 1 / state[2] - first[0][2] - second[0][2]
        assert inverse_colour > 0
        endpoints.append(1 / (4 * arb.pi()**2 * inverse_colour))
    return endpoints[0].union(endpoints[1])


def verify():
    reference_time = study["nominal"][2]
    target = comp["model"]["target"]
    left, centre, right = (target(n, reference_time) for n in (11430, 11431, 11432))
    gap = window["graph_gap"]()
    poles["verify_ratio_monotonicity"]()
    common_correction = q(".00806458").union(q(".00806466"))
    mass_error = arb(0, q("1e-8").upper())
    width_error = arb(0, q("1e-6").upper())
    for index, logs in enumerate(retune["log_sets"]):
        model, _, _ = retune["predictions"](logs, reference_time)
        correction = gluon_radiator(index, window["window_data"](model, gap))
        assert common_correction.contains(correction)
        m1, m2, tree_width, _, _ = poles["reference_readings"](logs)
        # Massless quarks and one colour-singlet real bidoublet at the quark port.
        # Widths are partial widths of a single charged pole divided by its mass.
        L, z, rho = model[:3]
        e1 = m1**2 / L
        component = 16 * z**2 / poles["right"](e1, poles["unit"], True)
        scalar_width = component / (96 * arb.pi() * rho)
        assert (18 * scalar_width).overlaps(tree_width)
        observed_width = tree_width * (1 + correction) * (1 + width_error)
        masses = (m1 * (1 + mass_error), m2 * (1 + mass_error))
        naive = poles["recover_charged"](*masses, observed_width)[0]
        corrected = poles["recover_charged"](*masses, observed_width / (1 + common_correction))[0]
        assert q("3.26600") < naive < q("3.26608") < centre
        assert left > corrected > right and corrected.contains(centre)
        print("case", index, "alpha_s/pi", correction.str(11),
              "quark/scalar width ratio", (tree_width * (1 + correction) / scalar_width).str(12))
        for label, value in (("naive target", naive), ("corrected target", corrected)):
            print(label, value.lower().str(14), value.upper().str(14))


if __name__ == "__main__":
    verify()

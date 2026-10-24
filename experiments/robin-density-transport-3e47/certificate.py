#!/usr/bin/env python3
"""Rational certificates for the density-weighted Robin transport application.

The owning Nicolas comparison consumes Polak W5/X1/X2 and the existing GAK
whole-law transport and reserve. This checks only the changed numerical
enclosures, not those analytic proofs, verified zeros, or the full RH claim.
"""

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path


Q = Fraction


def log_enclosure(value, terms=64):
    """Use the existing GAK25 positive-series enclosure for rational value >= 1."""
    value = Q(value)
    if value < 1:
        raise ValueError("positive-series input must be at least one")

    def reduced_enclosure(reduced):
        z = (reduced - 1) / (reduced + 1)
        lower = 2 * sum((z ** (2*j+1) / (2*j+1) for j in range(terms)), Q(0))
        upper = lower + 2*z**(2*terms+1) / ((2*terms+1)*(1-z*z))
        return lower, upper

    exponent = 0
    while value >= 2:
        value /= 2
        exponent += 1
    low, high = reduced_enclosure(value)
    low2, high2 = reduced_enclosure(Q(2))
    return low + exponent*low2, high + exponent*high2


def encode(value):
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value


def certificate(scenario="baseline"):
    checks = []
    capped = scenario == "capped-6e47"
    if scenario not in ("baseline", "capped-6e47"):
        raise ValueError("unknown scenario")

    def positive(name, margin):
        if margin <= 0:
            raise ArithmeticError((name, margin))
        checks.append({"name": name, "positive_margin": margin})

    height = 3 * 10**12
    clock_lower = 3 * 10**47 if capped else 8 * 10**46
    clock_upper = 6 * 10**47 if capped else 3 * 10**47
    root_upper = 8 * 10**23 if capped else 5478 * 10**20
    log3_low, log3_high = log_enclosure(Q(3))
    log10_low, log10_high = log_enclosure(Q(10))
    positive("log3_upper", Q(1099, 1000) - log3_high)
    positive("log10_upper", Q(2303, 1000) - log10_high)
    h_low = log3_low + 12*log10_low
    h_high = log3_high + 12*log10_high
    positive("height_log_lower", h_low - 28)
    positive("height_log_upper", Q(1437, 50) - h_high)
    exponential_prefix = sum((Q(84, 25)**j/factorial(j) for j in range(11)), Q(0))
    positive("log2874_upper_from_positive_exp_prefix", exponential_prefix - Q(1437, 50))
    positive("clock_log_lower", log_enclosure(Q(clock_lower))[0] - 107)
    derivative_factor = Q(64, 27) if capped else Q(3)
    positive("exponent_derivative_lower", 107 - derivative_factor*Q(321, 10))
    # Baseline L/(L-3s) decreases with L and increases with s. At the enclosing
    # endpoints it is exactly ten; the actual log inequalities are strict.
    if not capped and Q(107)/(107-3*Q(321, 10)) != 10:
        raise ArithmeticError("log-ratio endpoint is not ten")
    positive("upper_square_root", root_upper**2 - clock_upper)
    positive("upper_fourth_root", 10**48 - clock_upper)
    positive("lower_square_root", clock_lower - 10**46)
    s2 = Q(37, 25*10**12)
    if capped:
        ratio = Q(107)/(107-derivative_factor*Q(321, 10))
        if ratio != Q(45, 13):
            raise ArithmeticError("capped log-ratio transcription mismatch")
        positive("capped_log_ratio", Q(7, 2)-ratio)
        a_three_eighths_upper = Q(85, 10)*10**17
        height_log_five_thirds_lower = Q(16, 10)*10**23
        positive("upper_three_eighths_power", a_three_eighths_upper**4-root_upper**3)
        positive("lower_height_log_five_thirds_power", (28*height)**5-height_log_five_thirds_lower**3)
        density_bound = (
            2*10**12*s2
            + 144*Q(7, 2)*29**4*a_three_eighths_upper/height_log_five_thirds_lower
            + 30000*a_three_eighths_upper/height**2
        )
        density_upper = Q(1900)
        expected_density = Q(2276049253, 1200000)
    else:
        density_bound = (
            2*10**12*s2
            + Q(144*29**2*root_upper*10, height**2)
            + Q(30000*root_upper, height**2)
        )
        density_upper = Q(76000)
        expected_density = Q(9442616, 125)
    if density_bound != expected_density:
        raise ArithmeticError("density bound transcription mismatch")
    positive("complete_density_mass", density_upper-density_bound)

    count = 7 if capped else 9
    cutoff = Q(51, 20)
    a = cutoff/height
    x = a*a
    q0 = Q(7, 32)  # Existing all-frequency GAK41 bound, not re-proved here.
    q = Q(219, 1000)
    positive("inverse_moment_denominator", 1-2*count*x)
    positive("variance_denominator", 3-(count+3)*x)
    inverse_bound = 1/(1-2*count*x)
    variance_bound = count*x/(3-(count+3)*x)
    q_star = (q0 + a/(cutoff*(1-x)))/(1-x)
    positive("actual_real_part_complex_normalization", q-q_star)
    epsilon = count*a + count*x/(1-x)
    positive("whole_support_log", Q(1, 1000)-epsilon)
    positive("whole_support_lower", 1-epsilon-Q(99, 100))
    positive("whole_support_upper", Q(101, 100)-1/(1-epsilon))
    eta_bound = Q(1, 100)/107**2 + Q(100007, 100000)/10**23 + Q(178, 100)/10**30
    positive("existing_whole_future_prime_error", Q(1, 10**6)-eta_bound)
    head = Q(108, 107)*Q(527, 500)*Q(1, 20)*inverse_bound
    tail = Q(108, 107)*Q(10001, 10000)*q**count*density_upper
    transport = Q(260001, 500000)*root_upper*variance_bound
    elementary = 5*inverse_bound/10**23
    allowance = head+tail+transport+elementary
    for name, bound, upper in [
        ("verified_head", head, Q(532, 10000)),
        ("complete_weighted_tail", tail, Q(464, 10000) if capped else Q(889, 10000)),
        ("directed_transport", transport, Q(7014, 10000) if capped else Q(6175, 10000)),
        ("elementary_correction", elementary, Q(6, 10**23)),
    ]:
        positive(name, upper-bound)
    signed_upper = Q(161, 200) if capped else Q(19, 25)
    positive("signed_allowance", signed_upper-allowance)
    reserve = Q(67838284064811, 83593750000000)  # Existing GAK37.
    final_margin = Q(1, 200) if capped else Q(1, 20)
    positive("all_integer_margin", reserve-signed_upper-final_margin)
    # The selected-core difference is exactly 9/500. Its source reserve is
    # strictly larger than 0.778, so the analytic margin is strict.
    if not capped and Q(778, 1000)-Q(19, 25) != Q(9, 500):
        raise ArithmeticError("selected-source difference mismatch")
    return {
        "scope": "Changed rational enclosures only; external analytic premises retained.",
        "scenario": scenario,
        "clock": {"lower_excluded": clock_lower, "upper_included": clock_upper},
        "parameters": {"height": height, "count": count, "a": a, "q": q},
        "bounds": {
            "density_mass_upper": density_bound, "head": head, "tail": tail,
            "transport": transport, "elementary": elementary,
            "signed_allowance": allowance, "reserve": reserve,
            "reserve_minus_signed_allowance": reserve-allowance,
        },
        "checks": checks,
        "check_count": len(checks),
        "boundaries": [
            "Positive-series log/exp, published zero density, finite zero verification, GAK coefficient/transport and reserve are mathematical inputs.",
            "This program neither evaluates I_psi nor certifies the full analytic inference in Lean.",
            "The application controls a bounded real-clock interval, not every smaller integer, the least counterexample, or RH.",
            *(["The capped scenario additionally assumes ASC1 for every actual nontrivial zero; its upstream proof closure is not certified here."] if capped else []),
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scenario", choices=("baseline", "capped-6e47"), default="baseline")
    args = parser.parse_args()
    result = certificate(args.scenario)
    args.output.write_text(json.dumps(encode(result), indent=2)+"\n")
    print(json.dumps({
        "check_count": result["check_count"],
        "signed_allowance_display": float(result["bounds"]["signed_allowance"]),
        "reserve_minus_allowance_display": float(result["bounds"]["reserve_minus_signed_allowance"]),
        "output": str(args.output),
    }))


if __name__ == "__main__":
    main()

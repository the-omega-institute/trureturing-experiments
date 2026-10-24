#!/usr/bin/env python3
"""Exact scalar certificate for a truncated Fourier-law Robin application.

Fourier identities and the GAK analytic inputs are mathematical premises;
this program certifies their changed rational allowances, not RH or Lean.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from runpy import run_path

Q = Fraction
source = Path(__file__).resolve().parents[1] / "robin-density-transport-3e47" / "certificate.py"
log_enclosure = run_path(str(source))["log_enclosure"]


def encode(value):
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, list):
        return [encode(item) for item in value]
    return value


def certificate():
    checks = []
    def positive(name, value):
        if value <= 0:
            raise ArithmeticError((name, value))
        checks.append({"name": name, "positive_margin": value})

    height = 3 * 10**12
    epsilon = Q(1, 2000)
    radius = height * epsilon
    pi_upper = Q(22, 7)  # Classical pi < 22/7 is an analytic premise.
    positive("pi_squared_below_ten", 10-pi_upper**2)
    positive("density_tail_coefficient", 23-Q(64, 9)*pi_upper)
    positive("radius_above_two_pi", radius-2*pi_upper)
    positive("radius_below_two_to_31", 2**31-radius)
    log2_upper = log_enclosure(Q(2))[1]
    positive("log2_below_seven_tenths", Q(7, 10)-log2_upper)
    positive("log_radius_below_22", 22-31*log2_upper)
    positive("truncated_third_moment", 1100-(2*pi_upper**3+46*22))

    tau0 = Q(46, 3)/radius**3
    tau1 = Q(23)/radius**2
    tau2 = Q(46)/radius
    positive("retained_probability_mass", 1-tau0)
    v0 = Q(10, height**2)/(1-tau0)
    positive("inverse_moment_denominator", 1-2*v0)
    inverse = 1/(1-2*v0)
    variance = v0/(1-2*epsilon)
    # m <= 1+v0, so m^2 <= inverse and log(m) <= v0.
    positive("inverse_moment_bound", inverse-(1+v0)**2)
    positive("whole_log_support", Q(1, 1000)-epsilon-v0)
    positive("whole_support_lower", 1-epsilon-v0-Q(99, 100))
    positive("whole_support_upper", Q(101, 100)-1/(1-epsilon-v0))

    transform_numerator = (tau0+tau1/height+tau2/(2*height**2)
        +Q(1100)/(6*height**3*(1-epsilon)))
    transform_upper = transform_numerator/((1-tau0)*(1-v0))
    q = Q(1, 10**26)
    positive("complete_complex_tail", q-transform_upper)

    clock_lower = 8*10**46
    clock_upper = 17*10**47
    root_upper = Q(1304, 1000)*10**24
    positive("upper_square_root", root_upper**2-clock_upper)
    positive("lower_square_root", clock_lower-10**46)
    positive("clock_log_lower", log_enclosure(Q(clock_lower))[0]-107)
    eta_upper = Q(1,100)/107**2+Q(100007,100000)/10**23+Q(178,100)/10**30
    positive("existing_transport_eta", Q(1,10**6)-eta_upper)
    head = Q(108,107)*Q(527,500)*Q(1,20)*inverse
    tail = Q(108,107)*Q(10001,10000)*q*root_upper*Q(74,25*10**12)
    transport = Q(260001,500000)*root_upper*variance
    elementary = 5*inverse/10**23
    allowance = head+tail+transport+elementary
    reserve = Q(67838284064811,83593750000000)
    positive("verified_head", Q(532,10000)-head)
    positive("full_infinite_tail", Q(4,10**14)-tail)
    positive("directed_transport", Q(755,1000)-transport)
    positive("elementary_response", Q(6,10**23)-elementary)
    positive("complete_signed_allowance", Q(81,100)-allowance)
    positive("reserve_minus_rounded_allowance", reserve-Q(81,100)-Q(1,1000))
    positive("reserve_minus_exact_allowance", reserve-allowance-Q(3,1000))
    return {
        "scope": "Rational enclosures for BFR1-BFR10; analytic inputs are retained premises.",
        "clock": {"lower_excluded": clock_lower, "upper_included": clock_upper},
        "parameters": {"height": height, "epsilon": epsilon, "radius": radius,
            "root_upper": root_upper, "complex_tail_budget": q},
        "bounds": {"removed_mass_upper": tau0, "second_moment_upper": v0,
            "inverse_moment_upper": inverse, "variance_upper": variance,
            "complex_tail_upper": transform_upper, "verified_head": head,
            "full_tail": tail, "directed_transport": transport,
            "elementary": elementary, "signed_allowance": allowance,
            "reserve": reserve, "reserve_minus_allowance": reserve-allowance},
        "checks": checks, "check_count": len(checks),
        "boundaries": [
            "The cosine Fourier density, its first three Fourier cancellations, and Plancherel are mathematical inputs, not numerical integrations here.",
            "All five GAK analytic premises, finite zero verification, full response transport, and the unrestricted reserve are reused without rerunning their proofs.",
            "The law is truncated before exponential moments are taken. No exponential moment of the untruncated polynomial-tail law is asserted.",
            "This proves no unbounded signed estimate, least-counterexample bound, new zero-height verification, or Lean-certified Robin theorem.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = certificate()
    args.output.write_text(json.dumps(encode(result), indent=2)+"\n")
    print(json.dumps({"check_count": result["check_count"],
        "signed_allowance_display": float(result["bounds"]["signed_allowance"]),
        "reserve_minus_allowance_display": float(result["bounds"]["reserve_minus_allowance"]),
        "output": str(args.output)}))


if __name__ == "__main__":
    main()

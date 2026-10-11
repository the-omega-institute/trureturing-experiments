#!/usr/bin/env python3
"""Exact rational certificates for the bounded Robin smoothing application.

The analytic inputs and full-response transport are those of the owning
Library/ArithSums/nicolas2025comparison.md GAK application. This program
checks the changed sine-polynomial, support, normalization and allowance
inequalities; it does not certify those analytic inputs or prove RH.
"""

import argparse
import json
from fractions import Fraction
from math import factorial
from pathlib import Path


Q = Fraction


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

    def positive(name, margin):
        assert margin > 0, (name, margin)
        checks.append({"name": name, "positive_margin": margin})

    height = 3 * 10**12
    cutoff = Q(51, 20)
    count = 21
    a = cutoff / height
    x = a * a
    q0 = Q(7, 32)
    q = Q(219, 1000)
    p13 = sum(
        ((-1)**j * cutoff**(2*j + 1) / factorial(2*j + 1) for j in range(7)),
        Q(0),
    )
    positive("initial_sine_P13_upper", q0 * cutoff - p13 - Q(1, 10000))
    positive("alternating_tail_first_ratio", 1 - cutoff**2 / (16 * 17))

    pi_lower = Q(157, 50)
    t0 = Q(34, 25)
    t_upper = 1 / q0 - pi_lower
    p5 = t0 - t0**3 / 6 + t0**5 / 120
    p5_prime = 1 - t0**2 / 2 + t0**4 / 24
    positive("middle_segment_P5_concavity", 6 - t_upper**2)
    positive("tangent_point_in_middle_segment", t_upper - t0)
    positive("tangent_slope", q0 - p5_prime - Q(1, 1000))
    positive("tangent_intercept", q0 * pi_lower - (p5 - t0*p5_prime) - Q(3, 1000))
    positive("sine_sign_before_second_zero", 2*pi_lower - 1/q0)
    positive("cutoff_below_pi", pi_lower - cutoff)

    positive("inverse_moment_denominator", 1 - 2*count*x)
    positive("variance_denominator", 3 - (count + 3)*x)
    inverse_bound = 1 / (1 - 2*count*x)
    variance_bound = count*x / (3 - (count + 3)*x)
    q_star = (q0 + a/(cutoff*(1-x))) / (1-x)
    positive("all_real_parts_complex_tail", q - q_star)
    epsilon = count*a + count*x/(1-x)
    positive("whole_support_log_bound", Q(1, 1000) - epsilon)
    positive("whole_support_lower", 1 - epsilon - Q(99, 100))
    positive("whole_support_upper", Q(101, 100) - 1/(1-epsilon))

    clock_upper = 8 * 10**46
    root_upper = Q(2829, 1000) * 10**23
    positive("upper_clock_square_root", root_upper**2 - clock_upper)
    reserve = Q(67838284064811, 83593750000000)
    head = Q(108, 107) * Q(527, 500) * Q(1, 20) * inverse_bound
    tail = (
        2 * Q(108, 107) * Q(10001, 10000) * root_upper
        * Q(37, 25*10**12) * q**count
    )
    transport = Q(260001, 500000) * root_upper * variance_bound
    elementary = 5 * inverse_bound / 10**23
    allowance = head + tail + transport + elementary
    positive("complete_allowance_vs_reserve", reserve - allowance - Q(2, 1000))
    positive("displayed_margin_lower", reserve - allowance - Q(2407, 10**6))
    positive("displayed_margin_upper", Q(2408, 10**6) - (reserve - allowance))

    old_count = 23
    old_x = (Q(5, 2)/height)**2
    old_inverse = 1/(1-2*old_count*old_x)
    old_variance = old_count*old_x/(3-(old_count+3)*old_x)
    old_allowance = (
        Q(108, 107)*Q(527, 500)*Q(1, 20)*old_inverse
        + 2*Q(108, 107)*Q(10001, 10000)*root_upper
        * Q(37, 25*10**12)*Q(2501, 10000)**old_count
        + Q(260001, 500000)*root_upper*old_variance
        + 5*old_inverse/10**23
    )
    positive("old_allowance_does_not_pay_new_endpoint", old_allowance - reserve)
    positive("new_law_strictly_reduces_variance", old_variance - variance_bound)
    return {
        "scope": "Exact changed-parameter certificates only; external analytic premises retained.",
        "clock": {"lower_excluded": Q(7, 2)*10**46, "upper_included": clock_upper},
        "parameters": {"height": height, "count": count, "a": a, "q0": q0, "q": q},
        "bounds": {
            "inverse_moment": inverse_bound, "variance": variance_bound,
            "complex_tail_base": q_star, "whole_support_log": epsilon,
            "head": head, "tail": tail, "transport": transport,
            "elementary": elementary, "allowance": allowance,
            "reserve": reserve, "reserve_minus_allowance": reserve - allowance,
        },
        "checks": checks,
        "check_count": len(checks),
        "boundaries": [
            "Sine Taylor and concavity, pi lower bound, GAK analytic transport and reserve are mathematical inputs to the application.",
            "All actual zero real parts and multiplicities, both signs of ordinate, and the complete improper response are retained by GAK.",
            "This is not a kernel proof, an independent primary-proof audit, a least-counterexample bound or an unbounded RH estimate.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = certificate()
    args.output.write_text(json.dumps(encode(result), indent=2) + "\n")
    print(json.dumps({
        "check_count": result["check_count"],
        "reserve_minus_allowance_display": float(result["bounds"]["reserve_minus_allowance"]),
        "strict_margin": "2/1000",
        "output": str(args.output),
    }))


if __name__ == "__main__":
    main()

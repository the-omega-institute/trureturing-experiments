"""One fixed rational application of accepted theta rows; no source acquisition."""
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess

INPUT_COMMIT = "e5c5d8bf2b157e5c94af05ba724e18f4131b4c88"
PREFIX = "docs/reports/theta-mixed-matrix/"
N, H, EPS = 2**40, 3_000_000_000_000, F(1, 32)


def dyad(value):
    sig, exp = value["dyadic"]
    return F(int(sig)) * F(2) ** int(exp)


def encoded(value):
    value = F(value)
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def load(name):
    raw = subprocess.run(
        ["git", "show", INPUT_COMMIT + ":" + PREFIX + name],
        check=True, capture_output=True, text=True, close_fds=False,
    ).stdout
    return json.loads(raw)


def main():
    output = Path(__file__).with_name("verified-height-unwindowed-gap-result.json")
    if output.exists():
        raise RuntimeError("REFUSE_REPLAY: the fixed application already has a saved result")
    deriv = load("derivative-bandwidth-result.json")
    rows = load("joint-high-floor-result.json")
    exterior = load("sharper-exterior-result.json")
    if not (deriv["epsilon"] == "1/4" and deriv["bandwidth_N"] == 64
            and rows["bandwidth_N"] == 64 and rows["matrix_c"] == "3/8"
            and rows["grid_boxes"] == 128 and rows["radius"] == "3/2"
            and exterior["bandwidth_N"] == 64 and exterior["exterior_radius"] == "3/2"):
        raise ValueError("incompatible inherited fixed supplier metadata")
    checks = {}
    norms = [dyad(v) for v in deriv["derivative_norm_squared_upper"]]
    checks["saved_derivative_caps"] = len(norms) == 4 and all(
        0 < v < cap for v, cap in zip(norms, [F(1, 2), F(2), F(42), F(28000)]))
    sob = norms[0] + 3 * norms[1] + 3 * norms[2] + norms[3]
    kp = 64 * F(16, 15) * sob + 10 * (8 * 2432 + 2 * 3648)
    eta = 6 * kp / EPS * 29 / N
    checks["sobolev_cap"] = sob < 28133
    checks["commutator_cap"] = kp < 9 * 2**18
    checks["whole_projector_error_cap"] = eta < F(1, 80)
    checks["symbol_increment"] = dyad(deriv["m_N_over_2_upper"]) < F(22, 3)
    margins = []
    expected_left = F(0)
    for cell in rows["cells"]:
        left, right = F(cell["left"]), F(cell["right"])
        if left != expected_left or right - left != F(3, 256):
            raise ValueError("incomplete or reordered accepted row partition")
        s2 = dyad(cell["s2_lower"])
        if s2 <= 0:
            raise ValueError("nonpositive accepted theta squared lower endpoint")
        margins.append(dyad(cell["joint_floor_lower"]) + 6 * s2 - F(31, 64))
        expected_left = right
    checks["complete_partition"] = len(margins) == 128 and expected_left == F(3, 2)
    checks["new_fixed_row_floor"] = min(margins) > 0
    ext = dyad(exterior["exterior_W_floor_lower"])
    checks["same_exterior_floor"] = ext > F(31, 64)
    leakage = F(128, 9) * norms[2] / N**2
    checks["high_leakage_cap"] = leakage < F(1, 2**70)
    checks["high_form_gap"] = F(31, 64) - leakage > F(61, 128)
    log_ceiling = (H + 4).bit_length()
    r4 = F(256, H) * ((F(2, H**3) + F(1, 2 * H**2)) * log_ceiling
                      + (F(2, H**2) + F(1, H)) / (H + 4))
    checks["ordinate_tail_cap"] = log_ceiling == 42 and r4 < F(1, 2**109)
    checks["fixed_separation"] = 3 * H > 8 * N
    tail = N * F(2**42) * (1 - F(2 * N, H))**(-4) * r4
    checks["shifted_zero_tail_cap"] = tail < F(1, 2**19)
    contradiction_margin = EPS - F(197, 6400) - tail
    checks["strict_contradiction"] = contradiction_margin > 0
    checks["improves_owned_gap"] = F(15, 32) > F(921, 2000)
    result = {
        "scope": "Conditional paper scalar supplier for original even minimal form; no Lean, RH or cofinal certificate",
        "input_commit": INPUT_COMMIT,
        "inputs": ["derivative-bandwidth-result.json", "joint-high-floor-result.json", "sharper-exterior-result.json"],
        "epsilon": "1/32", "bandwidth_N": N, "verified_height": H,
        "source_acquisitions": 0, "old_producer_executions": 0, "parameter_scan": False,
        "row_comparisons": len(margins), "checks": checks,
        "minimum_row_margin": encoded(min(margins)),
        "sobolev_squared_cap": encoded(sob), "commutator_cap": encoded(kp),
        "whole_projector_error_cap": encoded(eta), "high_leakage_cap": encoded(leakage),
        "fourth_tail_cap": encoded(r4), "shifted_zero_tail_cap": encoded(tail),
        "strict_contradiction_margin": encoded(contradiction_margin),
        "required_paper_bridges": ["original same-form and accepted scalar supplier premises",
            "unwindowed actual minimal-form and paired zero-sum cutoff passage",
            "owned principal coefficient64 and fixed smoothed-trapezoid bounds",
            "whole-projector commutator application to the new high restriction",
            "published Platt-Trudgian inclusive verified head"],
        "passed_fixed_supplier": all(checks.values()),
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed_fixed_supplier": result["passed_fixed_supplier"],
                      "row_comparisons": len(margins), "failed": [k for k,v in checks.items() if not v]}))
    return 0 if result["passed_fixed_supplier"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

"""Certify the original-theta E=4 RV numerator obstruction with Arb."""
import argparse
import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import flint
from flint import acb, acb_poly, arb, ctx, fmpq


ROOT = Path(__file__).resolve().parents[2]
SUPPLIER = ROOT / "docs/reports/theta-mixed-matrix"


def definitions(path, names, scope):
    data = path.read_bytes()
    nodes = [node for node in ast.parse(data).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    if len(nodes) != len(names) or {node.name for node in nodes} != names:
        raise ValueError("Existing supplier definitions must be unique")
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), scope)
    return hashlib.sha256(data).hexdigest()


def endpoint(value):
    if not value.is_finite() or not value.is_exact():
        raise ValueError("Finite exact directed endpoint required")
    return [str(part) for part in value.man_exp()]


def enclosure(value):
    return {"display": str(value), "lower_dyadic": endpoint(value.lower()),
            "upper_dyadic": endpoint(value.upper())}


def coarse_bounds(value):
    denominator = 10**12
    lower = fmpq(math.floor(float(value.lower()) * denominator) - 1, denominator)
    upper = fmpq(math.ceil(float(value.upper()) * denominator) + 1, denominator)
    if not arb(lower) < value or not value < arb(upper):
        raise ValueError("Coarse rational enclosure failed")
    return lower, upper


def produce(precision):
    if precision < 192:
        raise ValueError("At least 192-bit precision required")
    ctx.prec = precision
    scope = {"acb": acb, "acb_poly": acb_poly, "arb": arb, "ctx": ctx,
             "fmpq": fmpq}
    hashes = {}
    hashes["derivative_bandwidth.py"] = definitions(
        SUPPLIER / "derivative_bandwidth.py", {"polys"}, scope)
    scope["POLYS"] = scope["polys"]
    hashes["theta_direct_action.py"] = definitions(
        SUPPLIER / "theta_direct_action.py", {"horner", "theta_jets"}, scope)
    path = SUPPLIER / "theta_translation_bounds.py"
    spec = importlib.util.spec_from_file_location("existing_theta_bounds", path)
    bounds = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bounds)
    hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    _, constants, _ = bounds.derivative_supplier(SUPPLIER / "derivative_bandwidth.py")

    cap = arb(3)
    start = (2 * cap).exp()
    # Phi(v) <= C0 exp(9v/2 - pi exp(2v)); the profile multiplier
    # 2 cosh(v/2)/(1+h sinh(v/2)^2) <= 2 exp(v/2) on v>=0.
    # The substituted tail is <= C0 int_start^inf y^2 exp(-pi*y)dy.
    tail = (2 * constants[0] * start.sqrt()
            * bounds.exponential_integral_upper(2, arb.pi(), start)).upper()
    if not tail > 0 or not tail < arb(2)**-1000:
        raise ValueError("Complete original-theta tail not certified")

    def sample(x):
        h = arb(4 * x * (1-x))

        def integrand(v, analytic):
            try:
                jets, factor = scope["theta_jets"](v, order=0, terms=12)
                return (2 * factor * jets[0] * (v/2).cosh()
                        / (1 + h * (v/2).sinh()**2))
            except ValueError:
                return acb("nan")

        total = acb(0)
        for j in range(12):
            total += acb.integral(integrand, acb(fmpq(j, 4)),
                                  acb(fmpq(j+1, 4)),
                                  abs_tol=arb(2)**-160,
                                  rel_tol=arb(2)**-160)
        if not total.is_finite() or not total.imag.contains(0):
            raise ValueError("Original real profile integral not enclosed")
        value = total.real + arb(0, tail)
        if not value > 0:
            raise ValueError("Positive profile sample not enclosed")
        return value

    # An independent normalization check includes the entire positive tail.
    p0_integral = sample(fmpq(0))
    if not p0_integral.contains(fmpq(1, 2)):
        raise ValueError("Original xi normalization P(0)=xi(1)=1/2 failed")
    p0 = arb(fmpq(1, 2))
    a, b = sample(fmpq(1, 4)), sample(fmpq(1, 2))
    alo, ahi = coarse_bounds(a)
    blo, bhi = coarse_bounds(b)
    # U4(q)/q^2 = P0(t^2-2)-4a*t+6b, t=q+1/q.
    discriminant = 16*a*a - 24*p0*b + 8*p0*p0
    lower = 16*alo*alo - 12*bhi + 2
    upper = 16*ahi*ahi - 12*blo + 2
    if not lower <= upper < 0 or not discriminant < 0:
        raise ValueError("Strict negative quartic discriminant not certified")
    if not 2*p0 - 8*a + 6*b > 0:
        raise ValueError("RV hypothesis U4(1) != 0 not certified")
    return {
        "source": "arXiv:2609.25564v1 equations (9), (10), (16)",
        "original_kernel": "Phi(v)=sum(4*pi^2*n^4*exp(9*v/2)-6*pi*n^2*exp(5*v/2))*exp(-pi*n^2*exp(2*v))",
        "runtime": {"python_flint": flint.__version__, "precision_bits": precision},
        "supplier_sha256": hashes,
        "finite_theta_terms": 12,
        "infinite_theta_series": "included by existing theta_jets tail enclosure",
        "integral_cap": 3,
        "integral_tail": enclosure(tail),
        "profile_endpoint_identity": "P(0)=xi(1)=1/2",
        "profile_endpoint_integral": enclosure(p0_integral),
        "profile_quarter": enclosure(a),
        "profile_center": enclosure(b),
        "coarse_rational_bounds": {
            "profile_quarter": [str(alo), str(ahi)],
            "profile_center": [str(blo), str(bhi)]},
        "discriminant": enclosure(discriminant),
        "discriminant_rational_bounds": [str(lower), str(upper)],
        "U4_at_one": enclosure(2*p0 - 8*a + 6*b),
        "consequence": "U4 has no unit-circle zero; every-even-degree unit-circle premise fails at E=4",
        "not_established": ["large-degree or subsequence obstruction", "RH or its negation", "Robin", "Lean proof"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=256)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = produce(args.precision)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()

"""Exact exceptional-layer certificates for the fixed FIB parent.

SPDX-License-Identifier: Apache-2.0
The original source and response identities are reused from bundle 177.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import flint
from flint import fmpq, fmpq_poly, nmod_mat, nmod_poly


ROOT = Path(__file__).resolve().parent
SUPPLIER = ROOT.parent / "fib_shared_response_177" / "verify.py"
SPEC = importlib.util.spec_from_file_location("fib177", SUPPLIER)
prior = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prior)
CHECKS = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)


def coefficient_field_reduction(poly, factor):
    terms = poly.to_dict()
    check("no e in source polynomial", all(k[2] == 0 for k in terms))
    coefficients = []
    for i in range(int(poly.degrees()[0]) + 1):
        row = {k[1]: v for k, v in terms.items() if k[0] == i}
        coefficients.append(fmpq_poly([row.get(j, 0)
                                      for j in range(max(row, default=0) + 1)]) % factor)
    while len(coefficients) > 1 and not coefficients[-1]:
        coefficients.pop()
    return coefficients


def remove_minus_one(coefficients, factor):
    quotient = [fmpq_poly([])] * (len(coefficients) - 1)
    quotient[-1] = coefficients[-1]
    for j in range(len(coefficients) - 2, 0, -1):
        quotient[j - 1] = (coefficients[j] - quotient[j]) % factor
    check("exact u+1 division", not (coefficients[0] - quotient[0]) % factor)
    return quotient


def mod_rational(value, prime):
    return int(value.numerator) * pow(int(value.denominator), -1, prime) % prime


def mod_poly(poly, beta, prime):
    out = [0] * (int(poly.degrees()[0]) + 1)
    for (i, j, k), value in poly.to_dict().items():
        check("no e in specialized polynomial", k == 0)
        out[i] = (out[i] + mod_rational(value, prime) * pow(beta, j, prime)) % prime
    return nmod_poly(out, prime)


def mod_coefficients(coefficients, beta, prime):
    return nmod_poly([sum(mod_rational(value, prime) * pow(beta, j, prime)
                         for j, value in enumerate(row)) % prime
                      for row in coefficients], prime)


def integers(poly):
    return [int(x) for x in poly]


def gaussian_determinant(rows, prime):
    matrix = [row[:] for row in rows]
    determinant = 1
    for j in range(len(matrix)):
        pivot = next((i for i in range(j, len(matrix)) if matrix[i][j] % prime), None)
        if pivot is None:
            return 0
        if pivot != j:
            matrix[j], matrix[pivot] = matrix[pivot], matrix[j]
            determinant = -determinant
        value = matrix[j][j] % prime
        determinant = determinant * value % prime
        inverse = pow(value, -1, prime)
        for i in range(j + 1, len(matrix)):
            scale = matrix[i][j] * inverse % prime
            matrix[i] = [(x - scale * y) % prime
                         for x, y in zip(matrix[i], matrix[j])]
    return determinant % prime


def bezout_unit(label, polynomial, source):
    gcd, left, right = polynomial.xgcd(source)
    check(label + " unit", gcd == 1 and left * polynomial + right * source == 1)
    return {"polynomial": integers(polynomial), "left": integers(left), "right": integers(right)}


def separation(index, factor, source, A, D, phi, eta, beta, prime):
    check("prime", prime in (101, 103))
    check("field specialization", sum(mod_rational(x, prime) * pow(beta, j, prime)
                                      for j, x in enumerate(factor)) % prime == 0)
    r = mod_coefficients(source, beta, prime)
    check("source degree preserved", r.degree() == 14)
    leading = int(r.leading_coefficient())
    r = r / r.leading_coefficient()
    units = {label: bezout_unit(label, mod_poly(poly, beta, prime), r)
             for label, poly in (("A", A), ("D", D), ("phi", phi))}
    squarefree = bezout_unit("source derivative", r.derivative(), r)
    denominator = mod_poly(phi, beta, prime)
    numerator = mod_poly(eta, beta, prime)
    inverse = nmod_poly(units["phi"]["left"], prime)
    response = numerator * inverse % r
    power, columns = nmod_poly([1], prime), []
    for _ in range(14):
        columns.append([int(power[i]) for i in range(14)])
        power = power * response % r
    rows = [[columns[j][i] for j in range(14)] for i in range(14)]
    determinant = gaussian_determinant(rows, prime)
    check("cyclic determinant", determinant != 0 and determinant == int(nmod_mat(rows, prime).det()))
    check("advertised cyclic determinant", determinant == {1: 98, 2: 54, 3: 58}[index])
    return {"prime": prime, "beta": beta, "source_leading": leading,
            "source_monic": integers(r), "units": units, "squarefree": squarefree,
            "eta_numerator": integers(numerator), "eta": integers(response),
            "cyclic_columns": columns, "cyclic_determinant": determinant}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    check("Python >=3.12", sys.version_info >= (3, 12))
    check("python-flint 0.9.0", flint.__version__ == "0.9.0")
    R, H, image_discriminant, phi, _, _, _, exceptional = prior.algebra()
    _, h_factors = H.factor()
    check("irreducible response H", len(h_factors) == 1 and h_factors[0][1] == 1
          and h_factors[0][0].degrees() == H.degrees())
    A, D, E = prior.read("A"), prior.read("D"), prior.read("E")
    eta = E.derivative(2) * prior.FE - E
    factors = {i: prior.univ(prior.read("support_factor_" + str(i))) for i in range(1, 7)}
    r_at_zero = prior.CTX.from_dict({k: value for k, value in R.to_dict().items() if k[0] == 0})
    check("no complex specialization collapse", exceptional["lc_R"].gcd(prior.univ(r_at_zero)).degree() == 0)
    s, u, b = prior.s, prior.u, prior.b
    bb = s.Rational(29, 23)
    ap = 23 * ((10 - 3 * bb) * u*u + (8 - 3 * bb) * u + 10 * bb - 32)
    bp = 23 * ((10 - 8 * bb) * u*u + (19 - 18 * bb) * u + 30 * bb - 36)
    check("zero child line excluded", s.resultant(ap, bp, u) == -378559400)
    check("parent not identically zero on source", s.Rational(-4, 7) * (-40)
          + s.Rational(175, 29) * (-4) == s.Rational(-260, 203))
    source_discriminant = prior.univ(R.resultant(R.derivative(0), 0))
    image_discriminant = prior.univ(image_discriminant)
    profile_polynomials = [exceptional["lc_R"], exceptional["content"], exceptional["A"],
                           exceptional["D"], exceptional["phi"], image_discriminant, source_discriminant]
    expected_profile = [[2, 2, 0, 2, 2, 2, 4], [0, 2, 0, 2, 2, 2, 2],
                        [0, 2, 0, 2, 2, 2, 2], [0, 1, 0, 2, 1, 0, 0],
                        [0, 0, 1, 0, 0, 0, 0], [0, 2, 0, 2, 2, 0, 2]]
    profiles = {}
    for i, factor in factors.items():
        _, decomposition = factor.factor()
        check("irreducible F" + str(i), len(decomposition) == 1 and decomposition[0][1] == 1
              and decomposition[0][0].degree() == factor.degree())
        row = [prior.valuation(p, factor) for p in profile_polynomials]
        check("exact profile F" + str(i), row == expected_profile[i - 1])
        profiles[str(i)] = {"degree": factor.degree(), "real_roots": prior.sturm(factor),
                            "valuations": row}
    check("real root counts", [profiles[str(i)]["real_roots"] for i in range(1, 7)]
          == [2, 2, 2, 4, 20, 24])
    check("48 smooth-image candidates", sum(profiles[str(i)]["real_roots"] for i in (4, 5, 6)) == 48
          and all(image_discriminant.gcd(factors[i]).degree() == 0 for i in (4, 5, 6)))
    source_certificates = {}
    for i, beta, prime in ((1, 59, 101), (2, 10, 101), (3, 44, 103)):
        factor = factors[i]
        coefficients = coefficient_field_reduction(R, factor)
        degree = len(coefficients) - 1
        check("actual source degree F" + str(i), degree == (14 if i == 1 else 16))
        if i == 2:
            coefficients = remove_minus_one(remove_minus_one(coefficients, factor), factor)
        elif i == 3:
            check("exact u squared divisor", not coefficients[0] and not coefficients[1])
            coefficients = coefficients[2:]
        check("reduced degree F" + str(i), len(coefficients) == 15 and bool(coefficients[-1]))
        check("quadratic coefficient field", all(row.degree() <= 1 for row in coefficients))
        certificate = separation(i, factor, coefficients, A, D, phi, eta, beta, prime)
        certificate["actual_source_degree"] = degree
        certificate["removed_factor"] = {1: "1", 2: "(u+1)^2", 3: "u^2"}[i]
        certificate["coefficient_field_source"] = [[str(x) for x in row] for row in coefficients]
        source_certificates[str(i)] = certificate
    result = {"input_commit": "7c45a21dbea5cdc0974f878a96dc64281cc5db07",
              "supplier_commit": "74d7549ac300dc45620ca7ae0a699e4d6b95845e",
              "supplier_sha256": hashlib.sha256(SUPPLIER.read_bytes()).hexdigest(),
              "input_sha256": prior.INPUT_HASHES,
              "profile_columns": ["lc_R", "content", "Res_R_A", "Res_R_D", "Res_R_phi",
                                  "image_discriminant", "source_discriminant"],
              "factor_profiles": profiles, "quadratic_separation": source_certificates,
              "arithmetic": "exact rational coefficient reduction, Sturm counts and finite-field nonzero certificates",
              "not_established_by_computation": ["coverage at all original source chart exceptions",
                  "original-root eight/nine-read decision", "Lean verification", "independent review", "remote CI"]}
    path = ROOT / "certificate.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        check("certificate matches", json.loads(path.read_text()) == result)
    print(json.dumps({"checks": len(CHECKS), "supplier_checks": len(prior.CHECKS),
                      "smooth_image_candidates": 48, "quadratic_candidates": 6,
                      "MM_cyclic_determinants": [98, 54, 58]}))


if __name__ == "__main__":
    main()

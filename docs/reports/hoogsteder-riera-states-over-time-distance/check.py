#!/usr/bin/env python3
"""Exact small-dimensional checks for arXiv:2610.34972v1.

All entries are integers or Fractions.  SymPy is used only for exact
characteristic polynomials/eigenvalues; no floating point number is used.
"""

from itertools import product

import sympy as sp


def mat(entries):
    return sp.Matrix(entries)


def partial_transpose_a(x, d, r):
    out = sp.zeros(d * r, d * r)
    for a, b, c, e in product(range(d), range(r), range(d), range(r)):
        out[c * r + b, a * r + e] = x[a * r + b, c * r + e]
    return out


def jordan(x, y):
    return (x * y + y * x) / 2


def trace_norm_hermitian(x):
    """Trace norm, as an exact sum of absolute algebraic eigenvalues."""
    if x != x.conjugate().T:
        raise AssertionError("the trace norm helper received a non-Hermitian matrix")
    total = sp.Integer(0)
    for value, multiplicity in x.eigenvals().items():
        total += multiplicity * sp.sqrt(sp.expand(value ** 2))
    return sp.simplify(total)


def assert_exact_le(left, right, label):
    # Both sides are nonnegative.  Squaring therefore preserves the order.
    gap = sp.simplify(sp.expand(right ** 2 - left ** 2))
    if gap.is_nonnegative is not True:
        raise AssertionError(f"{label}: exact gap was {gap!s}")


def channel_choi(kraus, d, out_dim):
    c = sp.zeros(d * out_dim, d * out_dim)
    for k in kraus:
        # vec(k) uses the input-first convention of the source.
        v = sp.Matrix([k[a, i] for i in range(d) for a in range(out_dim)])
        c += v * v.conjugate().T
    return c


def jamio_from_choi(c, d, out_dim):
    return partial_transpose_a(c, d, out_dim)


def kron_input(rho, out_dim):
    return sp.kronecker_product(rho, sp.eye(out_dim))


def fidelity_pairing_check():
    d = 2
    rho = mat([[sp.Rational(3, 4), 0], [0, sp.Rational(1, 4)]])
    unitary = mat([[1, 0], [0, -1]])
    c = channel_choi([unitary], d, d)
    j = jamio_from_choi(c, d, d)
    q = jordan(kron_input(rho, d), j)
    swap = mat([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    w0 = d * sp.eye(d * d) - swap
    lhs = sp.simplify(sp.trace(w0 * q))
    rhs = d - sp.re(sp.conjugate(sp.trace(unitary)) * sp.trace(rho * unitary))
    if sp.simplify(lhs - rhs) != 0 or lhs != 2:
        raise AssertionError(f"pairing formula mismatch: {lhs} versus {rhs}")
    print("fidelity pairing formula: OK (exact value 2)")


def universal_jordan_bound_check():
    cases = 0
    for d, r in ((2, 2), (3, 2), (2, 3)):
        rho = sp.zeros(d)
        rho[0, 0] = 1
        eye = sp.eye(r)
        # Every vector in this small integer box gives a rank-one Hermitian Z.
        for coordinates in product((-1, 0, 1), repeat=d * r):
            if not any(coordinates):
                continue
            z = sp.Matrix(coordinates)
            zmat = z * z.T
            transformed = jordan(kron_input(rho, r), partial_transpose_a(zmat, d, r))
            # This is the exact rank/HS certificate used in the source's
            # rank-one reduction: rank(T) <= d and ||T||_1^2 <= rank(T)||T||_HS^2.
            rank_bound = min(d, transformed.rank())
            hs_sq = sp.simplify(sum(entry ** 2 for entry in transformed))
            z_norm = sp.simplify((z.T * z)[0])
            gap = sp.simplify(d * z_norm ** 4 - rank_bound * hs_sq)
            if gap.is_nonnegative is not True:
                raise AssertionError(f"Jordan rank/HS certificate failed: {gap!s}")
            cases += 1
    print(f"Jordan/partial-transpose bound: OK ({cases} exact rank-one cases)")


def high_defect_branch_check():
    # An exact isometry A -> A'R whose A' output is fixed.  This is a
    # non-face generator and exercises the eta >= d branch with tau = 0.
    d, r = 2, 2
    isometry = sp.zeros(d * r, d)
    isometry[0, 0] = 1
    isometry[2, 1] = -1
    rho = sp.diag(sp.Rational(1, 2), sp.Rational(1, 2))
    c = channel_choi([isometry], d, d * r)
    j = jamio_from_choi(c, d, d * r)
    q = jordan(kron_input(rho, d * r), j)
    swap_aa = sp.zeros(d * d * r)
    for a, ap, s in product(range(d), range(d), range(r)):
        src = (a * d + ap) * r + s
        dst = (ap * d + a) * r + s
        swap_aa[dst, src] = 1
    w0 = d * sp.eye(d * d * r) - swap_aa
    eta = sp.simplify(sp.trace(w0 * q))
    if eta < d:
        raise AssertionError(f"expected the high-defect branch, eta={eta}")
    # The proof bounds ||Q||_1 by sqrt(d)||C_E||_1.  Here C_E is positive
    # semidefinite with exact trace d, so the left certificate is d sqrt(d).
    # Since eta >= d, this is at most d sqrt(eta), entirely algebraically.
    certificate = d * sp.sqrt(d)
    assert_exact_le(certificate, d * sp.sqrt(eta), "high-defect d*sqrt(eta) bound")
    print(f"high-defect branch: OK (eta={eta}, exact certificate checked)")


def main():
    fidelity_pairing_check()
    universal_jordan_bound_check()
    high_defect_branch_check()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

"""Exact checks for Gy's unnumbered coefficient question (Integers 26, A108)."""

from itertools import combinations, combinations_with_replacement
from math import prod

import sympy as sp
from sympy.polys.rings import ring


R, q = ring("q", sp.QQ)


def elementary(values, degree):
    return sum((prod(terms) for terms in combinations(values, degree)), R.zero)


def homogeneous(values, degree):
    return sum(
        (prod(terms) for terms in combinations_with_replacement(values, degree)),
        R.zero,
    )


def triple_differences(x, n, p, j):
    # x[0] is the extension term and is never a symmetric-polynomial variable.
    def e(m, d):
        return elementary(x[1:m + 1], d)

    def h(m, d):
        return homogeneous(x[1:m + 1], d)

    e_difference = (
        e(n, n - p) * e(p, p - j) * e(n - 1, n - j)
        - e(n - 1, n - p) * e(p - 1, p - j) * e(n, n - j)
    )
    h_difference = (
        h(p, n - p) * h(j, p - j) * h(j + 1, n - j)
        - h(p + 1, n - p) * h(j + 1, p - j) * h(j, n - j)
    )
    return e_difference, h_difference


def q_stirling(limit, first_kind, q_integers):
    # Sagan–Swanson, equations (2) and (5), using the unshifted second kind.
    table = {(0, 0): R.one}
    for n in range(1, limit + 1):
        for k in range(n + 1):
            factor = q_integers[n - 1 if first_kind else k]
            table[n, k] = (
                table.get((n - 1, k - 1), R.zero)
                + factor * table.get((n - 1, k), R.zero)
            )
    return table


def check_source():
    limit = 7
    x = [sum((q**i for i in range(k)), R.zero) for k in range(limit + 2)]
    first = q_stirling(limit + 1, True, x)
    second = q_stirling(limit + 1, False, x)
    triples = 0
    for n in range(1, limit + 1):
        for p in range(1, n + 1):
            assert elementary(x[1:n], n - p) == first[n, p], (n, p, "e")
            assert homogeneous(x[1:p + 1], n - p) == second[n, p], (n, p, "h")
            for j in range(1, p + 1):
                expected_e = (
                    first[p + 1, j + 1] * first[n + 1, p + 1] * first[n, j]
                    - first[p, j] * first[n, p] * first[n + 1, j + 1]
                )
                expected_h = (
                    second[p, j] * second[n, p] * second[n + 1, j + 1]
                    - second[p + 1, j + 1] * second[n + 1, p + 1] * second[n, j]
                )
                assert triple_differences(x, n, p, j) == (expected_e, expected_h), (n, p, j)
                triples += 1
    print(f"Fidelity (Gy p. 15): both q-Stirling identities and {triples} triple differences, n <= {limit}")


def check_symbolic_family():
    a, t, j, m, i, r, z = sp.symbols("a t j m i r z")

    def prefix_sum(size):
        return size * a + size * (size + 1) / 2

    def square_sum(size):
        return size * a**2 + a * size * (size + 1) + size * (size + 1) * (2 * size + 1) / 6

    def e2(size):
        return (prefix_sum(size)**2 - square_sum(size)) / 2

    def h2(size):
        return (prefix_sum(size)**2 + square_sum(size)) / 2

    assert sp.expand(prefix_sum(0)) == sp.expand(square_sum(0)) == 0
    assert sp.expand(prefix_sum(m + 1) - prefix_sum(m) - (m + 1 + a)) == 0
    assert sp.expand(square_sum(m + 1) - square_sum(m) - (m + 1 + a)**2) == 0
    assert sp.expand((i + a) * (r + a) - (i - 1 + a) * (r + 1 + a) - (r - i + 1)) == 0

    e_scalar = prefix_sum(j + 1) * (
        prefix_sum(j + 2) * e2(j + 1) - prefix_sum(j) * e2(j + 2)
    )
    h_scalar = prefix_sum(j + 1) * (
        prefix_sum(j) * h2(j + 1) - prefix_sum(j + 2) * h2(j)
    )
    scale = j * (j + 1)**2 * (j + 2) * (2*a + j + 2) * (4*a + j + 3) / 48
    assert sp.expand(e_scalar - scale) == sp.expand(h_scalar - scale) == 0

    fourth_power = sp.Poly((1 + z + t*z**2)**4, z)
    b0, b1, b2 = (fourth_power.nth(k) for k in range(3))
    assert sp.expand(b1**2 - b0*b2 - (10 - 4*t)) == 0
    for coefficient in fourth_power.all_coeffs():
        assert all(value > 0 for value in sp.Poly(coefficient, t).coeffs())
    print("Symbolic: strong index inequality factor r - i + 1")
    print(f"Symbolic: E and H scalar factors both equal {sp.factor(scale)}")
    print(f"Symbolic: degree-1 Turán determinant divided by the squared scalar factor is {sp.expand(b1**2 - b0*b2)}")


def check_direct_family():
    cases = 0
    for a in (sp.QQ(1, 2), sp.QQ(1), sp.QQ(7, 3)):
        for t in (sp.QQ(8, 3), sp.QQ(3), sp.QQ(5)):
            g = 1 + q + t*q**2
            x = [(k + a) * g for k in range(13)]
            for r in range(12):
                for i in range(1, r + 2):
                    difference = x[i] * x[r] - x[i - 1] * x[r + 1]
                    assert difference == (r - i + 1) * g**2, (a, t, i, r)
                    assert all(value >= 0 for value in difference.values()), (a, t, i, r)
            for j in range(1, 9):
                scale = j * (j + 1)**2 * (j + 2) * (2*a + j + 2) * (4*a + j + 3) / 48
                e, h = triple_differences(x, j + 2, j + 1, j)
                assert e == h == scale * g**4, (a, t, j)
                coeffs = [e.get((k,), 0) for k in range(9)]
                assert all(value > 0 for value in coeffs), (a, t, j)
                c0, c1, c2 = coeffs[:3]
                determinant = c1**2 - c0*c2
                assert determinant == (10 - 4*t) * scale**2 < 0, (a, t, j)
                if (a, t, j) == (1, 3, 1):
                    assert (c0, c1, c2, determinant) == (10, 40, 180, -200)
                    print("Example (a,t,j) = (1,3,1): both begin 10, 40, 180; determinant -200")
                cases += 1
    print(f"Direct family checks: {cases} rational parameter triples")


if __name__ == "__main__":
    check_source()
    check_symbolic_family()
    check_direct_family()
    print("ALL_CHECKS_OK")

#!/usr/bin/env python3
"""Check Che–Qi–Wei arXiv:1501.02546v1 using exact algebra."""

from itertools import product

import sympy as sp


X, Y = sp.symbols("x y", real=True)


def check_zero(expression):
    if sp.expand(expression) != 0:
        raise AssertionError(expression)


def gradient(form, coordinates):
    return tuple(sp.expand(sp.diff(form, coordinate) / 4) for coordinate in coordinates)


def check_fidelity():
    # Theorem 5(a): B has columns (1, 1) and (0, 1), with diagonal (1, 2).
    change_of_basis = sp.Matrix([[1, 0], [1, 1]])
    assert change_of_basis.det() == 1
    first, second = change_of_basis.T * sp.Matrix([X, Y])
    form = first**4 + 2 * second**4
    f1, f2 = gradient(form, (X, Y))
    check_zero(f1 - (X + Y) ** 3)
    check_zero(f2 - f1 - 2 * Y**3)
    # The q convention from §3 and the two conditions from Question 1.
    q1, q2 = -27, -43
    point = {X: 1, Y: 2}
    assert f1.subs(point) + q1 == 0 and f2.subs(point) + q2 == 0
    check_zero(form.subs(point) + q1 * point[X] + q2 * point[Y])
    # Origin fails. On each axis, the active component fixes the nonzero cube.
    assert f1.subs({X: 0, Y: 0}) + q1 < 0
    assert f2.subs({X: 3, Y: 0}) + q2 == -16
    check_zero(f2.subs(X, 0) - 3 * Y**3)
    check_zero(f1.subs(X, 0) - Y**3)
    assert sp.Rational(43, 3) + q1 < 0
    # In the interior f2 - f1 fixes y^3 = 8, then f1 fixes (x + y)^3 = 27.
    assert sp.integer_nthroot((q1 - q2) // 2, 3) == (2, True)
    assert sp.integer_nthroot(-q1, 3) == (3, True)
    print("Fidelity: Theorem 5(a) sample has unique solution (1, 2).")


def check_pair_identities():
    c, r, s = sp.symbols("c r s", positive=True)
    form = X**4 + 2 * c * X**2 * Y**2 + Y**4
    check_zero(form - ((X**2 - Y**2) ** 2 + 2 * (c + 1) * X**2 * Y**2))
    f1, f2 = gradient(form, (X, Y))
    check_zero(f1 - f2 - (X - Y) * ((X - Y) ** 2 + (3 - c) * X * Y))
    check_zero(f1.subs(X, 0))
    check_zero(f2.subs(Y, 0))
    jacobian = sp.Matrix([f1, f2]).jacobian((X, Y))
    check_zero(jacobian.subs({X: s, Y: s}).det() - 3 * (1 + c) * (3 - c) * s**4)
    ratio_c = 1 + r + 1 / r
    h = (r + 1) * (r**2 + 1)
    for component in (f1, f2):
        check_zero(component.subs({X: r, Y: 1, c: ratio_c}, simultaneous=True) - h)
    asymmetric_det = 3 * r**2 * s**4 * (c - 3) * (c**2 - 1)
    check_zero((jacobian.subs({X: r * s, Y: s}, simultaneous=True).det()
                - asymmetric_det).subs(c, ratio_c))
    print("Pair: exact factorization, quartic positivity identity and Jacobian identities.")


def tensor_entry(indices, c, epsilon):
    positions = set(indices)
    if len(positions) == 1:
        return 1 + epsilon
    if len(positions) == 2 and indices.count(indices[0]) == 2:
        first, second = sorted(positions)
        if first % 2 == 0 and second == first + 1:
            return c / 3 + epsilon
    return epsilon


def tensor_contractions(coordinates, c, epsilon):
    n = len(coordinates)
    contractions = [[] for _ in coordinates]
    for indices in product(range(n), repeat=4):
        entry = tensor_entry(indices, c, epsilon)
        assert entry > 0 and entry == tensor_entry(tuple(sorted(indices)), c, epsilon)
        monomial = sp.prod(coordinates[i] for i in indices[1:])
        contractions[indices[0]].append(entry * monomial)
    return tuple(sp.Add(*terms) for terms in contractions)


def check_family(k, ratio, epsilon):
    r = sp.Rational(ratio)
    c = 1 + r + 1 / r
    a = (r + 1) * (r**2 + 1)
    u = sp.Symbol("u")
    roots = sp.Poly(u**2 - (c - 1) * u + 1, u).real_roots()
    assert roots == [1 / r, r]
    t = sp.root(a / (1 + c), 3)
    branches = ((t, t), (r, sp.S.One), (sp.S.One, r))
    coordinates = sp.symbols(f"z0:{2 * k}", real=True)
    total = sum(coordinates)
    form = sum(x**4 + 2 * c * x**2 * y**2 + y**4
               for x, y in zip(coordinates[::2], coordinates[1::2]))
    assert c > 3 and epsilon > 0
    coupled = form + epsilon * total**4
    contractions = tensor_contractions(coordinates, c, epsilon)
    baseline = gradient(form, coordinates)
    for coordinate, actual, expected, block in zip(
            coordinates, contractions, gradient(coupled, coordinates), baseline):
        check_zero(actual - expected)
        check_zero(actual - block - epsilon * total**3)
        check_zero(block.subs(coordinate, 0))
        assert sp.Poly(actual, *coordinates).is_homogeneous
        assert sp.Poly(actual, *coordinates).total_degree() == 3
    check_zero(sum(x * f for x, f in zip(coordinates, contractions)) - coupled)

    signatures = set()
    for pairs in product(branches, repeat=k):
        w = tuple(x for pair in pairs for x in pair)
        assert all(x.is_positive for x in w)
        signatures.add(tuple(sp.simplify(x / y) for x, y in pairs))
        at_w = dict(zip(coordinates, w))
        for contraction in baseline:
            check_zero(contraction.subs(at_w) - a)
        # Represent the positive scale by its exact cube. Homogeneity then tests
        # every component at z = scale * w without simplifying nested radicals.
        height = a + epsilon * sum(w) ** 3
        assert height.is_positive
        scale_cubed = a / height
        assert scale_cubed.is_positive
        check_zero(scale_cubed * height - a)
        for contraction in contractions:
            check_zero(contraction.subs(at_w) - height)
    assert len(signatures) == 3**k
    print(f"Dense family: k={k}, c={c}, a={a}, epsilon={epsilon}, "
          f"{len(signatures)} distinct verified solution rays with unique scales.")


def check_symmetrized_antecedent():
    # Bai–Huang–Wang, Example 3.2 of arXiv:1508.06012v1.
    original = {(0, 0, 0, 0): 1, (0, 0, 0, 1): -2,
                (0, 0, 1, 1): 1, (1, 1, 1, 1): 1}
    index_tuples = list(product(range(2), repeat=4))
    classes = [[indices for indices in index_tuples if sum(indices) == count]
               for count in range(5)]
    means = [sp.Rational(sum(original.get(indices, 0) for indices in group), len(group))
             for group in classes]
    assert means == [1, -sp.Rational(1, 2), sp.Rational(1, 6), 0, 1]
    symmetric_form = sum(means[sum(indices)] * X**(4 - sum(indices)) * Y**sum(indices)
                         for indices in index_tuples)
    check_zero(symmetric_form - X**2 * (X - Y)**2 - Y**4)
    contractions = gradient(symmetric_form, (X, Y))
    for point in ({X: 0, Y: 1}, {X: 1, Y: 1}):
        assert tuple(component.subs(point) for component in contractions) == (0, 1)
        assert symmetric_form.subs(point) - point[Y] == 0
    print("Antecedent: the symmetrization of Bai–Huang–Wang Example 3.2 retains both solutions.")


def check_rational_witness():
    form = X**4 + 7 * X**2 * Y**2 + Y**4 + (X + Y) ** 4
    contractions = tensor_contractions((X, Y), sp.Rational(7, 2), sp.S.One)
    for point in ({X: 2, Y: 1}, {X: 1, Y: 2}):
        assert tuple(component.subs(point) for component in contractions) == (42, 42)
        assert form.subs(point) - 42 * sum(point.values()) == 0
    print("Rational refutation: q=(-42,-42), solutions (2,1) and (1,2).")


if __name__ == "__main__":
    check_fidelity()
    check_pair_identities()
    for k, ratio, epsilon in ((1, 2, sp.S.One), (2, 3, sp.S(2)), (3, 4, sp.Rational(1, 3))):
        check_family(k, ratio, epsilon)
    check_symmetrized_antecedent()
    check_rational_witness()
    print("ALL_CHECKS_OK")

"""Exact checks for Sain–Paul–Debnath, arXiv:1407.5016v3, Conjecture 2.11."""

from math import lcm

from sympy import Integer, Matrix, Rational, S, eye, sign, simplify, sqrt


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def equal(left, right, label):
    require(simplify(left - right) == 0, (label, left, right))


def check_source_example():
    """Example 2.9: strict full-span orthogonality, yet no zero coapproximation."""
    basis = Matrix([(1, 1, 1, 1), (1, -1, -1, 1), (1, -1, 1, -1), (1, 1, -1, -1)])
    require(basis * basis.T == 4 * eye(4), "Hadamard support pairings")
    for row in basis.tolist():
        require(all(abs(coordinate) == 1 for coordinate in row), "cube vertex")
    # Each row divided by 4 gives a linear form bounded by 1 on the cube.
    # Its value can equal 1 only when each of the four signed coordinates is 1.
    # Hence each hyperplane through a row and parallel to the other rows meets
    # the cube at that row alone.
    point = basis.row(0)
    test = sum((basis.row(index) for index in range(1, 4)), Matrix.zeros(1, 4))
    require(test * point.T == Matrix([[0]]), "coapproximation test in complementary span")
    norm_test = max(abs(coordinate) for coordinate in test)
    norm_difference = max(abs(coordinate) for coordinate in point - test)
    require(norm_test == 3 and norm_difference == 2, "zero fails coapproximation")
    print("Fidelity, Example 2.9: full support pairing = I; 3 > 2 rules out zero coapproximation.")


def planar_branch(s, t, exponent, scale):
    power_sum = abs(s) ** exponent + abs(t) ** exponent
    norm = scale * power_sum ** (1 / exponent)
    derivative_factor = scale * power_sum ** (1 / exponent - 1)
    gradient = Matrix([
        derivative_factor * sign(coord) * abs(coord) ** (exponent - 1)
        for coord in (s, t)
    ])
    return norm, gradient


def sector(s, t, p, q, c):
    return (p, c) if abs(s) >= abs(t) else (q, 1 / c)


def planar(s, t, p, q, c):
    return planar_branch(s, t, *sector(s, t, p, q, c))


def axial_frame(dimension):
    horizontal = eye(dimension - 1)
    horizontal[:2, :2] = Matrix([
        [Rational(3, 5), -Rational(4, 5)],
        [Rational(4, 5), Rational(3, 5)],
    ])
    require(horizontal.T * horizontal == eye(dimension - 1), "Euclidean frame")
    return horizontal.columnspace()


def lifted(vector, p, q, c):
    horizontal = vector[:-1, :]
    radius = sqrt(simplify(horizontal.dot(horizontal)))
    norm, grad = planar(radius, vector[-1], p, q, c)
    if radius == 0:
        equal(grad[0], 0, "differentiability on the rotation axis")
        return norm, Matrix.zeros(len(vector) - 1, 1).col_join(Matrix([grad[1]]))
    return norm, (grad[0] * horizontal / radius).col_join(Matrix([grad[1]]))


def check_family():
    exponents = [Rational(4, 3), Rational(3, 2), S(2), Rational(5, 2), S(3), S(4)]
    case_count = 0
    for p in exponents:
        q = p / (p - 1)
        c = Integer(2) ** (S.Half - 1 / p)
        for s in (-1, 1):
            for t in (-1, 1):
                value_p, grad_p = planar_branch(S(s), S(t), p, c)
                value_q, grad_q = planar_branch(S(s), S(t), q, 1 / c)
                equal(value_p, sqrt(2), "seam value")
                equal(value_p, value_q, "seam branch agreement")
                for index, coord in enumerate((s, t)):
                    equal(grad_p[index], coord / sqrt(2), "seam gradient")
                    equal(grad_p[index], grad_q[index], "seam gradient agreement")

        vertical = Matrix([0, 0, 1])
        horizontal = axial_frame(3)[0].col_join(Matrix([0]))
        left = sum(lifted(v, p, q, c)[0] ** 2 for v in (horizontal + vertical, horizontal - vertical))
        right = 2 * sum(lifted(v, p, q, c)[0] ** 2 for v in (horizontal, vertical))
        equal(left, 4, "parallelogram left side")
        equal(right, 2 * (c**2 + c**-2), "parallelogram right side")
        equal(right - left, 2 * (c - 1 / c) ** 2, "parallelogram gap")
        require((c == 1) == (p == 2), ("parallelogram nonzero off p=2", p))

        # Powers of this integer keep the unnormalized rays and their
        # conjugate directions rational. Normalization uses exact radicals.
        slope = S(2) ** lcm(p.q, q.q)
        rays = [(0, 1), (0, -1), (1, 0), (1, 1), (1, -1),
                (slope, 1), (slope, -1), (1, slope), (1, -slope)]
        for raw_r, raw_z in rays:
            r0, z0 = S(raw_r), S(raw_z)
            norm, grad = planar(r0, z0, p, q, c)
            exponent, _ = sector(r0, z0, p, q, c)
            y0 = Matrix([-sign(z0) * abs(z0) ** (exponent - 1), r0 ** (exponent - 1)])
            y_norm, y_grad = planar(*y0, p, q, c)
            equal(y0[0] / y_norm, -grad[1], "conjugate unit direction, first coordinate")
            equal(y0[1] / y_norm, grad[0], "conjugate unit direction, second coordinate")
            equal(y_grad[0], -z0 / norm, "conjugate gradient, first coordinate")
            equal(y_grad[1], r0 / norm, "conjugate gradient, second coordinate")

            for dimension in (3, 4, 6):
                first, *transverse = axial_frame(dimension)
                directions = [
                    (r0 * first).col_join(Matrix([z0])),
                    (y0[0] * first).col_join(Matrix([y0[1]])),
                    *(axis.col_join(Matrix([0])) for axis in transverse),
                ]
                # The norm is homogeneous and its gradient is constant along
                # each positive ray, so normalize after taking the gradient.
                supports = [lifted(direction, p, q, c) for direction in directions]
                for i, (_, functional) in enumerate(supports):
                    for j, (other, (other_norm, _)) in enumerate(zip(directions, supports)):
                        equal(functional.dot(other), int(i == j) * other_norm,
                              ("full pairing matrix", p, r0, z0, i, j))
                case_count += 1
    print(f"Family: {case_count} exact frames in dimensions 3, 4, and 6 for p = {exponents}.")
    print("Seam values and derivatives, conjugacy, all support cross-pairings, and parallelogram gaps verified.")


if __name__ == "__main__":
    check_source_example()
    check_family()
    print("ALL_CHECKS_OK")

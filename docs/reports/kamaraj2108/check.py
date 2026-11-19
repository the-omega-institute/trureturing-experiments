"""Exact checks for Filomat 37(3), §4, and arXiv:2108.05873v1."""

from sympy import I, Matrix, Rational, diag, expand, simplify, symbols, zeros


def adjoint(matrix, target_weight, source_weight):
    return source_weight.inv() * matrix.H * target_weight


def equal(left, right):
    assert left.shape == right.shape
    assert (left - right).applyfunc(simplify) == zeros(*left.shape)


def check_penrose(matrix, inverse, target_weight, source_weight):
    range_projection = matrix * inverse
    domain_projection = inverse * matrix
    equal(range_projection * matrix, matrix)
    equal(domain_projection * inverse, inverse)
    equal(adjoint(range_projection, target_weight, target_weight), range_projection)
    equal(adjoint(domain_projection, source_weight, source_weight), domain_projection)


def check_source_example():
    weight = diag(1, -1)
    a = Matrix([[1, 2], [0, 0]])
    b = Matrix([[2, 1], [0, 0]])
    a_adjoint = adjoint(a, weight, weight)
    b_adjoint = adjoint(b, weight, weight)
    equal(a_adjoint, Matrix([[1, 0], [-2, 0]]))
    equal(b_adjoint, Matrix([[2, 0], [-1, 0]]))
    a_inverse = -Rational(1, 3) * a_adjoint
    b_inverse = Rational(1, 3) * b_adjoint
    ab = a * b
    equal(ab, b)
    check_penrose(a, a_inverse, weight, weight)
    check_penrose(b, b_inverse, weight, weight)
    check_penrose(ab, b_inverse, weight, weight)
    equal(b_inverse * a_inverse, -Rational(1, 9) * b_adjoint)
    assert b_inverse != b_inverse * a_inverse
    print("Source Example 3.2: all three inverses satisfy the four equations; reverse order fails.")


def check_family():
    t = symbols("t", real=True, nonzero=True)
    weight = diag(1, 1, -1)
    e = Matrix([1, 0, 0])
    v = Matrix([1, t, -t])
    p = e * e.H * weight
    q = v * v.H * weight
    c = p * q
    x = q * p
    check_penrose(p, p, weight, weight)
    check_penrose(q, q, weight, weight)
    check_penrose(c, x, weight, weight)
    equal(c * x, p)
    equal(x * c, q)

    for matrix in (p, q, c, x):
        assert matrix[0, 0] == 1
        equal(matrix, matrix[:, 0] * matrix[0, :])

    p_adjoint = adjoint(p, weight, weight)
    q_adjoint = adjoint(q, weight, weight)
    first_range = p_adjoint * p * q
    second_range = q * q_adjoint * p_adjoint
    equal(first_range[:, 0], e)
    equal(second_range[:, 0], v)
    equal(q, v * q[0, :])
    equal(p_adjoint, e * p_adjoint[0, :])
    assert Matrix.hstack(e, v)[:2, :].det() == t

    block = Matrix.vstack(
        Matrix.hstack(q_adjoint * p_adjoint, q_adjoint * q),
        Matrix.hstack(p * p_adjoint, c),
    )
    equal(block, Matrix.vstack(v, e) * Matrix.hstack(e.H, v.H * weight))
    assert block[0, 0] == 1
    side_by_side = Matrix.hstack(p_adjoint, q)
    assert side_by_side.extract([0, 1], [0, 3]).det() == t
    equal(side_by_side, Matrix.hstack(e, v) * Matrix.vstack(
        Matrix.hstack(p_adjoint[0, :], zeros(1, 3)),
        Matrix.hstack(zeros(1, 3), q[0, :]),
    ))
    print("Rank-one family for every real t != 0: reverse order holds and both range inclusions fail.")
    print("Theorem 3.16 ranks for this family: block = 1, concatenation = 2.")


def check_gram_identities():
    a, b, c, d, e, f, g, h, s1, s2 = symbols("a b c d e f g h s1 s2", real=True)
    lines = Matrix([[a + I * b, c + I * d], [e + I * f, g + I * h]])
    gram = lines.H * diag(s1, s2) * lines
    determinant = lines.det()
    assert expand(gram.det() - s1 * s2 * determinant.conjugate() * determinant) == 0
    print("Dimension-two Gram determinant formula checked symbolically over complex coordinates.")


if __name__ == "__main__":
    check_source_example()
    check_family()
    check_gram_identities()
    print("ALL_CHECKS_OK")

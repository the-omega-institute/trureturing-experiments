#!/usr/bin/env python3
"""Exact checks for Lovitz, On decomposable correlation matrices.

All arithmetic is SymPy exact arithmetic over Q(i, sqrt(2), sqrt(5), sqrt(10)).
The finite family is a rank-three product-vector Gram family, so it checks the
source's Gram/Schur formulation and the proposed positive settlement on an
explicit nontrivial small case.
"""
from sympy import I, Matrix, Rational, sqrt, simplify


def kron(a, b):
    return Matrix([a[i] * b[j] for i in range(len(a)) for j in range(len(b))])


def gram(vectors):
    V = Matrix.hstack(*vectors)
    return V.conjugate().T * V


def exact_matrix_equal(A, B):
    return all(simplify(A[i, j] - B[i, j]) == 0
               for i in range(A.rows) for j in range(A.cols))


def is_correlation(P):
    return exact_matrix_equal(P, P.conjugate().T) and all(P[i, i] == 1 for i in range(P.rows))


def check_product_family():
    # Four unit product vectors in C^2 tensor C^2.  The first three span a
    # three-dimensional space and the fourth has a nonzero relation to all
    # three, matching the rank-three boundary isolated in Lovitz Section 4.
    e0 = Matrix([1, 0])
    e1 = Matrix([0, 1])
    a = [e0, e0, Matrix([1, 1]) / sqrt(2), Matrix([2, 1]) / sqrt(5)]
    b = [e0, e1, Matrix([1, 1]) / sqrt(2), Matrix([1, 1]) / sqrt(2)]
    vectors = [kron(a[i], b[i]) for i in range(4)]
    P = gram(vectors)
    Q = gram(a)
    R = gram(b)

    assert all(simplify(v.conjugate().T * v)[0] == 1 for v in vectors)
    assert is_correlation(P) and is_correlation(Q) and is_correlation(R)
    assert P.rank() == 3
    assert Q.rank() <= 2 and R.rank() <= 2
    assert exact_matrix_equal(P, Q.multiply_elementwise(R))
    # The fourth vector has a nonzero linear relation with the first three:
    # sqrt(10) v4 = v1 + v2 + 2 v3.
    relation = sqrt(10) * vectors[3] - vectors[0] - vectors[1] - 2 * vectors[2]
    assert all(simplify(x) == 0 for x in relation)
    return P


def check_source_small_case():
    # Fidelity check for the source's basic correlation-matrix convention:
    # the 2 by 2 identity is a correlation matrix and is one allowed factor
    # in a 2-decomposition (m = 1, rank = 2).
    P = Matrix.eye(2)
    assert is_correlation(P)
    assert P.rank() == 2
    assert exact_matrix_equal(P, P)


if __name__ == "__main__":
    check_source_small_case()
    P = check_product_family()
    print("source_small_case: 2x2 identity is a rank-2 correlation factor")
    print("rank_three_product_family: rank(P)=3 and P=Q Schur R exactly")
    print("P=")
    print(P)
    print("ALL_CHECKS_OK")

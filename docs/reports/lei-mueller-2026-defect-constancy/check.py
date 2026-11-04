#!/usr/bin/env python3
"""Exact checks for the proposed refutation of Lei–Müller Conjecture 7.11.

Adjacency entry (i,j) counts arrows j -> i, as in the source. No numerical
eigenvalue approximations are used. This is an independent computational check,
not a Lean proof. Requires Python 3 and SymPy 1.14.0.
"""

import math
import sys

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


x, z = sp.symbols("x z")


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def root_order(poly, root=1):
    poly = sp.Poly(poly, x, domain=sp.ZZ)
    check(not poly.is_zero, "root order of the zero polynomial is undefined")
    divisor = sp.Poly(x - root, x, domain=sp.ZZ)
    count = 0
    while poly.eval(root) == 0:
        poly, remainder = poly.div(divisor)
        check(remainder.is_zero, "nonexact polynomial division")
        count += 1
    return count


def integer_rank(matrix):
    """Integer matrix rank equals exact Q rank; free cokernel rank is n-rank.

    Cross-check with Smith form on matrices of order at most 12. Torsion
    invariant factors are not needed for the conjecture.
    """
    rank = matrix.to_DM().convert_to(sp.QQ).rank()
    if matrix.rows <= 12:
        diagonal = smith_normal_form(matrix, domain=sp.ZZ)
        smith_rank = sum(diagonal[i, i] != 0 for i in range(min(diagonal.shape)))
        check(rank == smith_rank, "Smith/Q rank mismatch")
    return rank


def invariants(matrix):
    size = matrix.rows
    check(matrix.cols == size, "matrix not square")
    charpoly = matrix.charpoly(x).as_expr()
    a = root_order(charpoly)
    rank = integer_rank(sp.eye(size) - matrix)
    b = size - rank
    return a, b, a - b, rank


def cycle(size):
    check(size >= 1, "empty cyclic group")
    return sp.Matrix(size, size, lambda i, j: int(i == (j + 1) % size))


def strongly_connected(matrix, sheets=1):
    """Two exact reachability searches on the actual constant lift."""
    size = matrix.rows
    outgoing = [[i for i in range(size) if matrix[i, j] > 0]
                for j in range(size)]
    incoming = [[j for j in range(size) if matrix[i, j] > 0]
                for i in range(size)]
    for neighbors, increment in [(outgoing, 1), (incoming, -1)]:
        reached = {(0, 0)}
        queue = [(0, 0)]
        for v, phase in queue:
            for w in neighbors[v]:
                node = (w, (phase + increment) % sheets)
                if node not in reached:
                    reached.add(node)
                    queue.append(node)
        check(len(reached) == size * sheets, "disconnected constant lift")


def family(q, c=3):
    u, v = sp.Matrix([1, -1, 0]), sp.Matrix([1, 1, -2])
    n = u * v.T
    check(n * n == sp.zeros(3), "N squared")
    check(n.rank() == 1, "N rank")
    check(n * sp.ones(3, 1) == sp.zeros(3, 1), "N kills constants")
    check(sp.ones(1, 3) * n == sp.zeros(1, 3), "N has zero column sums")
    check(n * (v / 6) == u and n * u == sp.zeros(3, 1), "Jordan chain")
    return c * sp.ones(3 * q) + sp.kronecker_product(sp.eye(3) + n, cycle(q))


def simple_replacement(matrix):
    """Six-vertex classes: internal shifts 1,2,3; cross shifts 0,...,a_ij-1."""
    size = matrix.rows
    result = sp.zeros(6 * size)
    for i in range(size):
        for j in range(size):
            count = int(matrix[i, j])
            check(1 <= count <= 5, "entry outside replacement range")
            shifts = [1, 2, 3] if i == j else range(count)
            if i == j:
                check(count == 3, "diagonal not three")
            for t in range(6):
                for shift in shifts:
                    result[6 * i + (t + shift) % 6, 6 * j + t] += 1
    check(set(result) <= {0, 1}, "parallel edge in replacement")
    check(all(result[i, i] == 0 for i in range(result.rows)), "loop in replacement")
    projection = sp.kronecker_product(sp.eye(size), sp.ones(6, 1))
    check(result * projection == projection * matrix, "constant-class quotient")
    check(projection.T * result == matrix * projection.T, "sum-zero complement")
    # At every vertex: a 2-cycle through a different class, and a 3-cycle
    # through two other distinct classes, all at the same within-class index.
    for i in range(size):
        j, k = (i + 1) % size, (i + 2) % size
        for t in range(6):
            a, b, c = 6 * i + t, 6 * j + t, 6 * k + t
            check(result[b, a] == result[a, b] == 1, "missing 2-cycle")
            check(result[b, a] == result[c, b] == result[a, c] == 1,
                  "missing 3-cycle")
    return result


def main():
    print("Python", sys.version.split()[0], "SymPy", sp.__version__)
    # Audit the printed matrix in Example 7.2 without adopting its erroneous
    # claim of a double eigenvalue 1. Its actual defect is zero.
    source = sp.Matrix([[1, 0, 1], [1, 2, 1], [1, 1, 2]])
    check(invariants(source) == (1, 1, 0, 2), "printed Example 7.2 matrix")
    denominator = (sp.eye(3) - z * source).det()
    check(sp.expand(denominator - (1 - z) * (1 - 4 * z + 2 * z ** 2)) == 0,
          "printed matrix zeta denominator")
    print("SOURCE_EXAMPLE_7_2_DISCREPANCY: charpoly="
          + str(sp.factor(source.charpoly(x).as_expr())))
    print("SOURCE_EXAMPLE_7_2_DISCREPANCY: a=1 b=1 delta=0; source states a=2 delta=1")
    # Fidelity: Example 7.10, in its connected case r=3, p=2. The lift is a
    # directed cycle on 3*2^m vertices, a=b=1, defect zero.
    for m in range(4):
        sheets = 2 ** m
        lift = sp.kronecker_product(cycle(3), cycle(sheets))
        size = 3 * sheets
        check(invariants(lift) == (1, 1, 0, size - 1), "source Example 7.10")
        check(sp.Poly(lift.charpoly(x).as_expr() - (x ** size - 1), x).is_zero,
              "directed-cycle characteristic polynomial")
        strongly_connected(cycle(3), sheets)
        print("FIDELITY Example 7.10 p=2 r=3 m=%s vertices=%s a=1 b=1 delta=0 SC=yes"
              % (m, size))
    print("FAMILY c=3: p h q m sheets vertices a b delta rank(I-lift) SC")
    for p in [2, 3]:
        for h in [1, 2]:
            q, c = p ** h, 3
            matrix = family(q, c)
            size, r = 3 * q, 9 * q + 1
            check(min(matrix) >= 1, "nonpositive entry")
            check(all(matrix[i, i] == c for i in range(size)), "diagonal")
            check(matrix * sp.ones(size, 1) == r * sp.ones(size, 1), "row sums")
            check(sp.ones(1, size) * matrix == r * sp.ones(1, size), "column sums")
            expected_cp = (x - r) * (x - 1) ** 2 * sp.div(x ** q - 1, x - 1)[0] ** 3
            check(sp.Poly(matrix.charpoly(x).as_expr() - expected_cp, x).is_zero,
                  "base characteristic polynomial")
            previous = None
            for m in range(h + 2):
                sheets = p ** m
                power = matrix ** sheets
                # det(x I - A tensor C_s) = det(x^s I - A^s).
                # x^s-1 has a simple zero at 1 in characteristic zero.
                cp = power.charpoly(x).as_expr()
                a = root_order(cp)
                rank_small = integer_rank(sp.eye(size) - power)
                b = size - rank_small
                delta = a - b
                d = math.gcd(q, sheets)
                check((a, b, delta) == (3 * d - 1, 2 * d - 1, d), "family law")
                strongly_connected(matrix, sheets)
                # Direct lift cross-check for each required (p,h) at m=1,
                # including its full characteristic polynomial and exact rank.
                if m == 1:
                    lift = sp.kronecker_product(matrix, cycle(sheets))
                    check(invariants(lift) == (a, b, delta, size * sheets - b),
                          "direct lift invariants")
                    check(sp.Poly(lift.charpoly(x).as_expr()
                                  - cp.subs(x, x ** sheets), x).is_zero,
                          "direct determinant identity")
                if previous is not None:
                    check(delta >= previous, "source Theorem 7.6 monotonicity")
                previous = delta
                print(p, h, q, m, sheets, size * sheets, a, b, delta,
                      size * sheets - b, "yes", flush=True)
            simple = simple_replacement(matrix)
            for m in range(h + 2):
                strongly_connected(simple, p ** m)
            print("SIMPLE_STRUCTURE p=%s h=%s vertices=%s 0/1=yes diagonal=0 "
                  "quotient=A complement=invariant SC_all_tested_levels=yes"
                  % (p, h, simple.rows), flush=True)
            if (p, h) == (2, 1):
                a0, b0, delta0, _ = invariants(simple)
                a1, b1, delta1, _ = invariants(sp.kronecker_product(simple, cycle(p)))
                check(delta1 - delta0 >= p - 1, "simple defect increase")
                print("SIMPLE_EXACT p=2 h=1: base a=%s b=%s delta=%s; "
                      "level1 a=%s b=%s delta=%s; increase=%s >= 1"
                      % (a0, b0, delta0, a1, b1, delta1, delta1 - delta0), flush=True)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

"""Exact checks for a four-layer Hermitian portal with a rank-one cut."""

from itertools import permutations, product

import sympy as sp


def portal(sign):
    return sp.Matrix([
        [5, sign*sp.I, 1, 1],
        [-sign*sp.I, 6, 1, 1],
        [1, 1, 7, sp.I],
        [1, 1, -sp.I, 8],
    ])


def loop_weight(matrix, indices):
    successors = indices[1:] + indices[:1]
    return sp.prod(matrix[i, j] for i, j in zip(indices, successors))


def check_pencils_and_cycles():
    plus, minus = portal(1), portal(-1)
    mass = sp.diag(1, 2, 3, 4)
    s = sp.symbols("s", real=True)
    x = sp.symbols("x:4", real=True)
    # The x_i can be shifted independently, and can include the spectral parameter.
    expected = ((x[0]*x[1]-s**2)*(x[2]*x[3]-s**2)
                - s**2*(x[0]+x[1])*(x[2]+x[3]))
    for matrix in (plus, minus):
        assert matrix == matrix.H
        assert all(matrix[i, i] > sum(abs(matrix[i, j]) for j in range(4) if j != i)
                   for i in range(4))
        off_diagonal = matrix - sp.diag(*matrix.diagonal())
        assert sp.expand((sp.diag(*x)-s*off_diagonal).det(method="berkowitz") - expected) == 0

    for order in (1, 2, 3, 4):
        for indices in product(range(4), repeat=order):
            if len(set(indices)) < 4:
                assert sp.re(loop_weight(plus, indices)-loop_weight(minus, indices)) == 0
    assert [loop_weight(matrix, (0, 1, 2, 3)) for matrix in (plus, minus)] == [-1, 1]
    traces = [((mass*matrix-matrix*mass)**4).trace() for matrix in (plus, minus)]
    assert traces == [848, 736]
    print("Commutator fourth traces (+, -):", traces)


def check_separated_source():
    # Position labels are their coordinates in units of the common spacing.
    cycles = ((0, 1, 3, 7), (0, 1, 7, 3), (0, 3, 1, 7))
    sorted_lengths = []
    for cycle in cycles:
        lengths = tuple(abs(a-b) for a, b in zip(cycle, cycle[1:]+cycle[:1]))
        sorted_lengths.append(tuple(sorted(lengths)))
        # Strict rearrangement: largest mass on shortest edge, then proceed in order.
        mass_order = tuple(3-sorted(lengths).index(length) for length in lengths)
        assert len(set(lengths)) == 4
        if cycle == cycles[0]:
            assert mass_order == (3, 2, 1, 0)
            assert [loop_weight(portal(k), mass_order) for k in (1, -1)] == [-1, 1]
    assert sorted_lengths == [(1, 2, 4, 7), (1, 3, 4, 6), (2, 3, 6, 7)]
    masses = (2, sp.sqrt(3), sp.sqrt(2), 1)
    rates = [sum(length*m for length, m in zip(lengths, masses)) for lengths in sorted_lengths]
    assert rates[1]-rates[0] == sp.sqrt(3)-1
    assert rates[2]-rates[0] == 2+sp.sqrt(3)+2*sp.sqrt(2)
    # Check the six oriented spacetime cycles split into the three reversal pairs.
    assert {(0,)+p for p in permutations((1, 3, 7))} == {
        order for cycle in cycles for order in (cycle, (0,)+cycle[:0:-1])
    }
    print("Separated four-source minimum exponents:", rates)


if __name__ == "__main__":
    check_pencils_and_cycles()
    check_separated_source()

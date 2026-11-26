"""Finite algebra checks of the Higgs frame and Wronskian section."""
from fractions import Fraction
from itertools import combinations
from math import comb, prod
import sympy as sp


def color_frame(columns):
    basis = []
    for column in range(columns.cols):
        column = columns[:, column]
        orthogonal = column.copy()
        for unit in basis:
            orthogonal -= unit*(unit.H*column)[0]
        basis.append(orthogonal / sp.sqrt((orthogonal.H*orthogonal)[0]))
    first = sp.Matrix.hstack(*basis)
    assert first.shape == (5, 4)
    last = sp.Matrix(
        [(-1)**i * first.extract([j for j in range(5) if j != i], range(4)).det().conjugate()
         for i in range(5)])
    return first.row_join(last)


def verify_frame():
    base = sp.eye(5)
    vectors = sp.Matrix.hstack(base[:, 1], base[:, 2], base[:, 1]+base[:, 4], sp.I*base[:, 3])
    gauge = sp.eye(5)
    gauge[:2, :2] = sp.Matrix([[3*sp.I, 4], [-4, -3*sp.I]])/5
    assert gauge.H*gauge == base and gauge.det() == 1
    frame = color_frame(vectors)
    assert frame.H*frame == base and frame.det() == 1
    assert color_frame(gauge*vectors) == gauge*frame


def verify_wronskians():
    # For monomial rows and divided derivatives, every minor is a Vandermonde.
    for size in range(1, 8):
        for rank in range(1, size+1):
            maximum = 0
            for powers in combinations(range(size), rank):
                coefficients = sp.Matrix([[comb(i, j) for j in range(rank)] for i in powers])
                numerator = prod(b-a for a, b in combinations(powers, 2))
                assert coefficients.det() * prod(sp.factorial(j) for j in range(rank)) == numerator
                maximum = max(maximum, sum(powers) - rank*(rank-1)//2)
            assert maximum == rank*(size-rank)
    # At a complex root the projected derivatives lose rank.
    z = sp.symbols("z")
    v = sp.Matrix([1, z, z**2, z**3])
    rows = sp.Matrix([[1, 0, 1, 0], [0, 1, 0, 1]])
    w = (rows*sp.Matrix.hstack(v, v.diff(z))).det().factor()
    assert sp.simplify(w-(z*z+1)**2) == 0
    assert (rows*sp.Matrix.hstack(v, v.diff(z))).subs(z, sp.I).rank() == 1


def probability_budget():
    sites = 625
    dims = [4*sites*d for d in (5, 10)]
    roots = sum(n*n//4 for n in dims)
    assert roots == 195312500 and roots < 2**28
    # 4*(1-|z|^2)^3 is relative to normalized disk area.
    assert sum(Fraction(4*(-1)**j * comb(3, j), j+1) for j in range(4)) == 1
    kappa = 2**115
    conditional_bound = 8 * (4*(1+kappa))**5 * 2**(23*14)
    assert roots*conditional_bound * Fraction(1, 2**960) < Fraction(1, 2**17)
    print(f"Two Wronskians have at most {roots} exceptional complex values together.")


if __name__ == "__main__":
    verify_frame()
    verify_wronskians()
    probability_budget()

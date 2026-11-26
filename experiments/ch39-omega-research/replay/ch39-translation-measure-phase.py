"""Translation action on interpolation jets and the finite two-coordinate probability bound."""
from fractions import Fraction
from itertools import combinations
from math import comb
import sympy as sp


def interpolation_jet(labels, argument, rank):
    variable = sp.Symbol("v")
    rows = []
    for i in range(len(labels)):
        polynomial = sp.prod(variable-value for j, value in enumerate(labels) if j != i)
        rows.append([sp.diff(polynomial, variable, j).subs(variable, argument) / sp.factorial(j)
                     for j in range(rank)])
    return sp.Matrix(rows)


def volume_components(columns):
    return sp.Matrix([columns.extract(indices, range(columns.cols)).det()
                      for indices in combinations(range(columns.rows), columns.cols)])


def certify():
    sites = 3
    internal = 2
    labels = [sp.Rational(x+1, 7)+3*a for x in range(sites) for a in range(internal)]
    s = sp.I/5
    order = [(i+internal) % len(labels) for i in range(len(labels))]
    shift = sp.eye(len(labels)).extract(order, range(len(labels)))
    columns = sp.Matrix([[1, 2], [0, 1], [1, sp.I], [sp.I, 0], [-2, 1], [1, -1]])
    proj = columns * (columns.H*columns).inv() * columns.H
    jet = interpolation_jet(labels, s, columns.cols)
    shifted = interpolation_jet([labels[i] for i in order], s, columns.cols)
    assert shifted == shift*jet
    original = (proj*jet).applyfunc(sp.cancel)
    moved = ((shift*proj*shift.T)*shifted).applyfunc(sp.cancel)
    assert moved == shift*original and volume_components(original) != sp.zeros(15, 1)
    for i, point in enumerate(labels):
        evaluation = interpolation_jet(labels, point, 1)
        assert evaluation[i] != 0 and sum(v != 0 for v in evaluation) == 1
    # Pair density is 12*(1-s-t)^2 for squared radii s,t in the triangle.
    u, w = sp.symbols("u w")
    assert sp.integrate(12*(1-u-w)**2, (w, 0, 1-u), (u, 0, 1)) == 1
    n = 625
    root_count = sum((4*n*d)**2//4 for d in (5, 10))
    density_bound = 24*(4*(1+2**115))**5*2**(23*14)
    failure = density_bound*(comb(n, 2)+n*n*root_count)*Fraction(1, 2**992)
    assert failure < Fraction(1, 2**28)
    print("Translation covariance, interpolation rank, disk density and finite separation verified.")


if __name__ == "__main__":
    certify()

#!/usr/bin/env python3
"""Exact helper checks for Beshkov Conjecture 2.16.

All values are fractions.  The checks enumerate rational input grids and
rational boundary points and verify the stated simplex-boundary image formula
for d=1,2,3,4.  They are helper checks for the planned general argument.
"""
from fractions import Fraction as F
from itertools import product


def relu(z):
    return z if z > 0 else F(0)


def hidden(x):
    d = len(x)
    return tuple([relu(v) for v in x] + [relu(F(1) - sum(x, F(0)))])


def output_from_hidden(a):
    t = F(1) - sum(a, F(0))
    return tuple(relu(v + 2 * t) for v in a)


def network(x):
    return output_from_hidden(hidden(x))


def in_boundary(y):
    return all(v >= 0 for v in y) and sum(y, F(0)) <= 1 and (
        sum(y, F(0)) == 1 or any(v == 0 for v in y)
    )


def preimage_of_boundary(y):
    """Exact inverse used in the source's surjectivity calculation."""
    d = len(y) - 1
    q = sum(y, F(0))
    if q == 1:
        a = y
    else:
        zero = next(i for i, v in enumerate(y) if v == 0)
        t = (F(1) - q) / (2 * d - 1)
        a = tuple(v + (F(2) * t if i != zero else F(0)) for i, v in enumerate(y))
    assert sum(a, F(0)) >= 1
    assert sum(a, F(0)) == 1 or any(v == 0 for v in a)
    if q == 1 or a[-1] == 0:
        x = a[:-1]
    else:
        zero = next(i for i in range(d) if a[i] == 0)
        t = sum(a, F(0)) - 1
        x = tuple(-t if i == zero else a[i] for i in range(d))
    return x


def grid_points(d, denominator):
    vals = [F(i, denominator) for i in range(-denominator, 2 * denominator + 1)]
    return product(vals, repeat=d)


def boundary_grid(m, denominator):
    vals = [F(i, denominator) for i in range(denominator + 1)]
    for y in product(vals, repeat=m):
        if in_boundary(y):
            yield y


def run_dimension(d):
    forward = 0
    for x in grid_points(d, 2):
        y = network(tuple(x))
        assert in_boundary(y), (d, x, y)
        forward += 1
    inverse = 0
    for denominator in (2, 3):
        for y in boundary_grid(d + 1, denominator):
            x = preimage_of_boundary(y)
            assert all(F(-1) <= v <= F(2) for v in x), (d, y, x)
            assert network(x) == y, (d, y, x, network(x))
            inverse += 1
    return forward, inverse


def main():
    totals = []
    for d in range(1, 5):
        counts = run_dimension(d)
        totals.append(counts)
        print(f"d={d}: exact forward-grid points={counts[0]}, boundary-inverse points={counts[1]}")
    # The source's displayed indexing has d inputs but a sum through d+1;
    # this check uses the dimensionally consistent sum through the d inputs.
    assert network((F(0),)) == (F(1), F(0))
    assert network((F(1),)) == (F(0), F(1))
    print("checked d=1..4: P([-1,2]^d) grid outputs lie in the simplex boundary,")
    print("and every enumerated rational boundary point has an exact inverse in [-1,2]^d.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

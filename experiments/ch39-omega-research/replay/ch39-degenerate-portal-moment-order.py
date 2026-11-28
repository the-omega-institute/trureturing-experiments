"""Exact witnesses and determinant checks for single-mass portal inversion."""

from itertools import combinations_with_replacement
from math import prod
from sympy import Matrix


def moments(values, k):
    return sum(x**k for x in values)


def sign_changes(x, y):
    points = sorted(set(x + y))
    signs = []
    for a, b in zip(points, points[1:]):
        twice_t = a + b
        h = sum(2*z > twice_t for z in x) - sum(2*z > twice_t for z in y)
        if h:
            signs.append(1 if h > 0 else -1)
    return sum(s != t for s, t in zip(signs, signs[1:]))


def main():
    x, y = (1, 12), (9, 10)
    assert moments(x, 3) == moments(y, 3) == 1729
    assert (moments(x, 4), moments(y, 4)) == (20737, 16561)
    for n in range(1, 7):
        points = range(1, n+1)
        jacobian = Matrix([[k*z**(k-1) for z in points] for k in range(3,n+3)])
        determinant = prod(range(3,n+3)) * prod(z*z for z in points)
        determinant *= prod(b-a for a in points for b in points if a < b)
        assert jacobian.det() == determinant
    for n in range(1, 5):
        lists = list(combinations_with_replacement(range(-2, 3), n))
        assert all(sign_changes(x, y) <= n-1 for i, x in enumerate(lists) for y in lists[i:])
        assert len({tuple(moments(x,k) for k in range(3,n+3)) for x in lists}) == len(lists)
    print("Cube collision, fourth power sums, Vandermonde determinants and signed grid verified.")


if __name__ == "__main__":
    main()

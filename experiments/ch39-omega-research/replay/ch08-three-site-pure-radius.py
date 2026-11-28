"""Check the five three-site cyclic difference patterns and exact pure witnesses."""

from collections import Counter
from itertools import combinations
from math import gcd
from sympy import I, Rational, S, conjugate, simplify, sqrt


def geometric_type(sites, dimension):
    for middle in sites:
        left, right = (x for x in sites if x != middle)
        if (left + right - 2 * middle) % dimension == 0:
            order = dimension // gcd(middle - left, dimension)
            return {3: "triangle", 4: "square"}.get(order, "progression")
    has_diameter = any(2 * (left - right) % dimension == 0
                       for left, right in combinations(sites, 2))
    return "diameter" if has_diameter else "distinct"


def verify_geometries():
    # Direct ordered-difference multiplicities, independent of the midpoint test.
    signatures = {
        "distinct": (1, 1, 1, 1, 1, 1),
        "diameter": (1, 1, 1, 1, 2),
        "progression": (1, 1, 2, 2),
        "square": (2, 2, 2),
        "triangle": (3, 3),
    }
    counts = Counter()
    for dimension in range(3, 37):
        for sites in combinations(range(dimension), 3):
            differences = Counter((x-y) % dimension for x in sites for y in sites if x != y)
            name = geometric_type(sites, dimension)
            assert tuple(sorted(differences.values())) == signatures[name]
            counts[name] += 1
    assert set(counts) == set(signatures)
    return counts


def check_witness(dimension, sites, amplitudes, expected):
    norm = simplify(sum(z * conjugate(z) for z in amplitudes))
    correlations = [S.Zero] * dimension
    for x, a in zip(sites, amplitudes):
        for y, b in zip(sites, amplitudes):
            correlations[(x-y) % dimension] += a * conjugate(b) / norm
    assert simplify(correlations[0]) == 1
    energy = sum(z * conjugate(z) for z in correlations[1:])
    assert simplify(energy - expected) == 0


def leakage_window(left, right):
    # left and right are the original two positive, unnormalized amplitudes.
    norm = left**2 + right**2
    product_sq = (left * right)**2
    return tuple(Rational(product_sq, product_sq + norm * width**2)
                 for width in (left+right, left-right))


if __name__ == "__main__":
    counts = verify_geometries()
    check_witness(7, (0, 1, 3), (1, 1, 1), Rational(2, 3))
    check_witness(8, (0, 1, 4), (1, 1, I), Rational(4, 9))
    check_witness(5, (0, 1, 2), (3, 2, -1), Rational(25, 98))
    square_phase = (1-I) / sqrt(2)
    check_witness(4, (0, 1, 2), (1, square_phase, I), 0)
    check_witness(12, (1, 4, 7), (3, 2*square_phase, I), Rational(8, 49))
    check_witness(3, (0, 1, 2), (3, 2, -1), Rational(1, 98))
    check_witness(9, (1, 4, 7), (1, (-1 + I*sqrt(3))/2, 1), 0)

    lower, upper = leakage_window(2, 1)
    assert (lower, upper) == (Rational(4, 49), Rational(4, 9))
    check_witness(3, (0, 1, 2), (6, 3, -2), 0)
    check_witness(3, (0, 1, 2), (2, -1, 2), 0)
    assert leakage_window(1, 1) == (Rational(1, 9), 1)
    print(f"Checked {sum(counts.values())} triples. Pattern counts: {dict(counts)}")
    print("All seven pure witnesses and leakage endpoints verified.")

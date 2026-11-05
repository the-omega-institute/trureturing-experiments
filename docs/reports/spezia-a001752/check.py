#!/usr/bin/env python3
"""Exact finite checks for Spezia's OEIS A001752 COMMENTS conjecture.

Source: https://oeis.org/A001752/internal, revision 119.
This program is numerical evidence on the stated finite ranges, not a proof.
"""
from math import comb
import platform
import sympy as sp

SOURCE_DATA = (1,4,11,24,46,80,130,200,295,420,581,784,1036,1344,1716,2160,2685,3300,4015,4840,5786,6864,8086,9464,11011,12740,14665,16800,19160,21760,24616,27744,31161,34884,38931,43320,48070,53200,58730,64680,71071,77924,85261)


def series_coefficients(limit):
    """Invert (1+t)(1-t)^5 by exact coefficient recursion."""
    t = sp.Symbol("t")
    denominator = sp.Poly((1 + t) * (1 - t)**5, t)
    d = [int(denominator.nth(j)) for j in range(7)]
    assert d == [1, -4, 5, 0, -5, 4, -1]
    values = [1]
    for k in range(1, limit + 1):
        values.append(-sum(d[j] * values[k-j] for j in range(1, min(k, 6)+1)))
    return values


def source_closed_formula(k):
    u = (k + 3)**2
    return ((u - 1) * (u - 3)) // 48


def floor_matrix(n):
    # Python // by positive 2 is mathematical floor, also for negative numerators.
    return sp.Matrix(n, n, lambda i, j: (j - i + 1) // 2)


def distance_sum(n):
    return sum((n-d) * ((d*d)//4) for d in range(1, n))


def main():
    print("Python", platform.python_version(), "SymPy", sp.__version__)
    a = series_coefficients(200)
    assert a[:len(SOURCE_DATA)] == list(SOURCE_DATA)
    print(f"FIDELITY: all {len(SOURCE_DATA)} published terms, offset 0, agree with the generating function")
    for k in range(201):
        alternating_binomial = sum((-1)**(k-j) * comb(j+4, 4) for j in range(k+1))
        assert a[k] == source_closed_formula(k) == alternating_binomial
    print("FIDELITY: generating-function, floor closed formula, and alternating binomial sum agree for k=0..200")
    assert list(floor_matrix(4)[0, :]) == [0, 1, 1, 2]
    assert list(floor_matrix(4)[:, 0]) == [0, 0, -1, -1]
    print("FLOOR: M(4) first row [0,1,1,2], first column [0,0,-1,-1]")
    x = sp.Symbol("x")
    for n in range(2, 33):
        matrix = floor_matrix(n)
        c = int(matrix.charpoly(x).as_poly().nth(n-2))
        minors = sum(int(matrix.extract([i,j], [i,j]).det())
                     for i in range(n) for j in range(i+1,n))
        for i in range(n):
            for j in range(i+1, n):
                assert int(matrix.extract([i,j], [i,j]).det()) == ((j-i)**2)//4
        assert c == minors == distance_sum(n)
        assert c >= 0
        expected = a[n-3] if n >= 3 else 0
        assert abs(c) == expected
        print(f"n={n:2d} charpoly_coeff={c} minor_sum={minors} distance_sum={distance_sum(n)} expected={expected}")
    for n in range(3, 204):
        assert distance_sum(n) == a[n-3]
    print("SETTLEMENT CHECK: exact characteristic polynomials n=3..32; distance sums n=3..203")
    print("LIMIT: finite exact checks do not prove the universal conjecture")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

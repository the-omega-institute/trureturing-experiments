#!/usr/bin/env python3
"""Exact certificates for Sra, arXiv:2609.39979v1, Section 4(b,d).

Python 3 standard library only. Matrix data and source fidelity target:
Suvrit Sra, A Koteljanskii inequality for permanents, Theorem 1.2.
All arithmetic below is integral or fractions.Fraction; no numerical spectra.
"""
from fractions import Fraction as F
from itertools import combinations, permutations
from math import prod


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def multiply(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def inverse(a):
    n = len(a)
    rows = [list(map(F, row)) + identity(n)[i] for i, row in enumerate(a)]
    for k in range(n):
        pivot = next(i for i in range(k, n) if rows[i][k])
        rows[k], rows[pivot] = rows[pivot], rows[k]
        t = rows[k][k]
        rows[k] = [x / t for x in rows[k]]
        for i in range(n):
            if i != k:
                t = rows[i][k]
                rows[i] = [x - t * y for x, y in zip(rows[i], rows[k])]
    return [row[n:] for row in rows]


def determinant(a):
    if not a:
        return F(1)
    return sum((-1)**j * a[0][j] * determinant(
        [row[:j] + row[j + 1:] for row in a[1:]]) for j in range(len(a)))


def cycle_count(p):
    seen, count = set(), 0
    for i in range(len(p)):
        if i not in seen:
            count += 1
            j = i
            while j not in seen:
                seen.add(j)
                j = p[j]
    return count


def cycle_coefficients(a, subset):
    """Coefficient of alpha^c, including fixed points among the cycles."""
    subset = tuple(sorted(subset))
    result = [0] * (len(subset) + 1)
    for p in permutations(range(len(subset))):
        result[cycle_count(p)] += prod(a[subset[i]][subset[p[i]]]
                                       for i in range(len(subset)))
    return result


def subsets(n):
    return [frozenset(s) for k in range(n + 1)
            for s in combinations(range(n), k)]


def polynomial_multiply(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def polynomial_subtract(a, b):
    return [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
            for i in range(max(len(a), len(b)))]


def evaluate(a, x):
    return sum(c * x**i for i, c in enumerate(a))


def check_inverse_m(a, m):
    n = len(a)
    assert all(x > 0 for row in a for x in row)
    assert multiply(a, m) == multiply(m, a) == identity(n)
    assert inverse(a) == m
    assert all(m[i][j] <= 0 for i in range(n) for j in range(n) if i != j)


def derivative(poly, variable):
    result = {}
    for exponent, coefficient in poly.items():
        if exponent[variable]:
            reduced = list(exponent)
            reduced[variable] -= 1
            result[tuple(reduced)] = coefficient * exponent[variable]
    return result


def main():
    # Published Theorem 1.2: a fidelity check distinct from the new Hessian reduction.
    a = [[13, 1, 7, 7, 7], [1, 13, 7, 7, 7], [1, 1, 7, 7, 7],
         [1, 1, 1, 7, 7], [1, 1, 1, 1, 7]]
    m = [[F(x, 84) for x in row] for row in
         [[7, 0, -7, 0, 0], [0, 7, -7, 0, 0], [0, 0, 14, -14, 0],
          [0, 0, 0, 14, -14], [-1, -1, 0, 0, 14]]]
    check_inverse_m(a, m)
    assert m[0][2] != 0 and m[2][0] == 0  # excludes diagonal symmetrizability
    r, s, t, u = map(frozenset, [(2, 3, 4), (0, 2, 3, 4), (1, 2, 3, 4), range(5)])
    coeff = {z: cycle_coefficients(a, z) for z in subsets(5)}
    permanents = {z: sum(c) for z, c in coeff.items()}
    assert coeff[r] == [0, 56, 147, 343]
    assert coeff[s] == coeff[t] == [0, 546, 2051, 2940, 4459]
    assert coeff[u] == [0, 3444, 21308, 48125, 51940, 57967]
    assert [permanents[z] for z in [r, s, t, u]] == [546, 9996, 9996, 182784]
    assert permanents[r] * permanents[u] == 99800064
    assert permanents[s] * permanents[t] == 99920016
    assert F(permanents[r] * permanents[u], permanents[s] * permanents[t]) == F(832, 833)
    assert all(p > 0 for p in permanents.values())
    print('FIDELITY Theorem 1.2: A*M=M*A=I; all A entries positive; inverse off-diagonal <=0')
    print('FIDELITY cycle sums:', [coeff[z][1:] for z in [r, s, t, u]])
    print('FIDELITY permanents [R,S,T,U]:', [permanents[z] for z in [r, s, t, u]])
    print('FIDELITY product ratio: 832/833; difference: -119952')

    # Build the full polynomial from all 32 principal permanents, then differentiate.
    poly = {(5 - len(z),) + tuple(int(i in z) for i in range(5)): F(1, p)
            for z, p in permanents.items()}
    assert len(poly) == 32 and all(sum(e) == 5 for e in poly)
    q = poly
    for variable in [3, 4, 5]:
        q = derivative(q, variable)
    expected = {(2, 0, 0, 0, 0, 0): F(1, 546),
                (1, 1, 0, 0, 0, 0): F(1, 9996),
                (1, 0, 1, 0, 0, 0): F(1, 9996),
                (0, 1, 1, 0, 0, 0): F(1, 182784)}
    assert q == expected
    h = [[derivative(derivative(q, i), j).get((0,) * 6, F(0))
          for j in range(6)] for i in range(6)]
    h3 = [row[:3] for row in h[:3]]
    assert determinant(h3) == -F(1, 7597732565090304)
    assert sum(h3[i][i] for i in range(3)) == F(1, 273)
    basis = [[1, 0], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0]]
    transpose = [list(row) for row in zip(*basis)]
    restriction = multiply(multiply(transpose, h), basis)
    assert restriction == [[F(1, 273), F(1, 4998)], [F(1, 4998), F(1, 91392)]]
    assert restriction[0][0] > 0 and determinant(restriction) == F(1, 20783363328)
    print('ITEM (d) differentiated polynomial:', {e: str(c) for e, c in q.items()})
    print('ITEM (d) Hessian 3x3:', [[str(x) for x in row] for row in h3])
    print('ITEM (d) det(H3)=-1/7597732565090304; trace(H3)=1/273')
    print('ITEM (d) positive two-plane: first pivot=1/273; determinant=1/20783363328')

    # New symmetric positive witness: B=6*(I-C_cycle/4)^(-1).
    b = [[7, 2, 1, 2], [2, 7, 2, 1], [1, 2, 7, 2], [2, 1, 2, 7]]
    mb = [[F(1, 6) if i == j else -F(1, 24) if (i-j) % 4 in (1, 3) else F(0)
           for j in range(4)] for i in range(4)]
    check_inverse_m(b, mb)
    assert b == [list(row) for row in zip(*b)]  # positive diagonal similarity is I
    rr, ss, tt, uu = map(frozenset, [(0, 2), (0, 1, 2), (0, 2, 3), range(4)])
    bc = {z: cycle_coefficients(b, z) for z in subsets(4)}
    assert bc[rr] == [0, 1, 49]
    assert bc[ss] == bc[tt] == [0, 8, 63, 343]
    assert bc[uu] == [0, 48, 257, 882, 2401]
    assert bc[frozenset()] == [1]
    assert all(all(c >= 0 for c in cs) and any(c > 0 for c in cs) for cs in bc.values())
    difference = polynomial_subtract(polynomial_multiply(bc[uu], bc[rr]),
                                     polynomial_multiply(bc[ss], bc[tt]))
    factorized = polynomial_multiply([0, 0, 1, 1], [-16, 1617, 2401])
    assert difference == factorized + [0]
    assert factorized == [0, 0, -16, 1601, 4018, 2401]
    endpoint = F(1, 200)
    assert endpoint > 0 and endpoint not in (F(1, 2), F(1))
    upper_bound = evaluate([-16, 1617, 2401], endpoint)
    assert upper_bound == -F(314199, 40000) < 0
    # The nonconstant coefficients are positive, so this upper bound holds
    # for every 0<alpha<=endpoint, without sampling the interval.
    assert all(c > 0 for c in [1617, 2401])
    assert evaluate(difference, endpoint) == -F(63153999, 320000000000)
    # Known positive orders on every pair of subsets for this witness.
    for alpha in [F(1, 2), F(1)]:
        for v in bc:
            for w in bc:
                assert evaluate(bc[v | w], alpha) * evaluate(bc[v & w], alpha) >= \
                       evaluate(bc[v], alpha) * evaluate(bc[w], alpha)
    print('ITEM (b) B=', b)
    print('ITEM (b) B inverse=', [[str(x) for x in row] for row in mb])
    print('ITEM (b) symmetric; diagonal symmetrizer=I; all 16 entries strictly positive')
    print('ITEM (b) alpha-permanent coefficients [R,S,T,U]:', [bc[z] for z in [rr, ss, tt, uu]])
    print('ITEM (b) difference=alpha^2*(alpha+1)*(2401*alpha^2+1617*alpha-16)')
    print('ITEM (b) bracket <= -314199/40000 for 0<alpha<=1/200')
    print('ITEM (b) difference at alpha=1/200: -63153999/320000000000')
    print('ITEM (b) all 256 subset pairs satisfy supermodularity at alpha=1/2 and alpha=1')
    print('INTERPRETATION (d): positive definite restriction has dimension two, excluding Lorentzianity')
    print('INTERPRETATION (b): universal affirmative extension is false; remaining alpha ranges unclassified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

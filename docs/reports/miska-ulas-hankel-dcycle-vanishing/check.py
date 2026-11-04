#!/usr/bin/env python3
"""Exact finite checks for Miska–Ulas, arXiv:1904.03395v1, §8 and Conj. 9.5(3).

These checks do not prove the universal conjecture. Only integer and rational
arithmetic is used. The source's polynomial recurrence is compared with its
factorial expansion and with direct weighted enumeration of permutations.
"""
from fractions import Fraction
from itertools import permutations
from math import comb, factorial, prod


def exact_quotient(a, b):
    assert b != 0
    q, r = divmod(a, b)
    assert r == 0, (a, b, r)
    return q


def det(a):
    """Fraction-free Bareiss elimination, including row swaps and empty matrices."""
    n = len(a)
    if n == 0:
        return 1
    b = [list(row) for row in a]
    assert all(len(row) == n for row in b)
    sign, previous = 1, 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if b[i][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            b[k], b[pivot_row] = b[pivot_row], b[k]
            sign = -sign
        pivot = b[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                b[i][j] = exact_quotient(b[i][j] * pivot - b[i][k] * b[k][j], previous)
        for i in range(k + 1, n):
            b[i][k] = 0
        previous = pivot
    return sign * b[-1][-1]


def det_leibniz(a):
    n = len(a)
    return sum((-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
               * prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))


def source_polynomials(d, maximum):
    """Coefficient vectors, constant coefficient first, from the literal recurrence."""
    out = []
    for n in range(maximum + 1):
        if n < d:
            out.append([0] * n + [1])
        else:
            a = [0] + out[n - 1]
            factor = prod(range(n - d + 1, n))
            for j, c in enumerate(out[n - d]):
                a[j] += factor * c
            out.append(a)
    return out


def closed_polynomial(d, n):
    a = [0] * (n + 1)
    for j in range(n // d + 1):
        a[n - d * j] = exact_quotient(factorial(n), factorial(n - d * j) * factorial(j) * d ** j)
    return a


def evaluate(a, x):
    y = 0
    for c in reversed(a):
        y = y * x + c
    return y


def cycles(p):
    used, lengths = set(), []
    for i in range(len(p)):
        if i in used:
            continue
        j, size = i, 0
        while j not in used:
            used.add(j)
            size += 1
            j = p[j]
        lengths.append(size)
    return lengths


def moment(d, j):
    return exact_quotient(factorial(d * j), d ** j * factorial(j))


def rising(a, n):
    return prod((a + j for j in range(n)), start=Fraction(1))


def main():
    # Cross-check the determinant implementation against its definition, including
    # zero pivots, negative entries, singular matrices and the empty determinant.
    samples = [[], [[0]], [[3]], [[0, 1], [1, 0]], [[1, 2], [2, 4]]]
    samples += [[[((i + 2) * (j + 5) + i * i - 3 * j) % 11 - 5
                  for j in range(n)] for i in range(n)] for n in range(2, 7)]
    for a in samples:
        assert det(a) == det_leibniz(a)
    print('Determinant implementation: 10 exact Leibniz comparisons OK')

    fidelity = 0
    for d in range(2, 13):
        polys = source_polynomials(d, 48)
        for n, a in enumerate(polys):
            assert a == closed_polynomial(d, n), (d, n)
            fidelity += 1
    print(f'Fidelity: section 8 factorial expansion agrees with literal recurrence in {fidelity} cases (2<=d<=12, 0<=n<=48)')

    permutation_cases = 0
    for n in range(9):
        counts = {d: [0] * (n + 1) for d in range(2, 9)}
        for p in permutations(range(n)):
            lengths = cycles(p)
            fixed = lengths.count(1)
            nontrivial = {k for k in lengths if k != 1}
            for d in range(2, 9):
                if nontrivial <= {d}:
                    counts[d][fixed] += 1
        for d, a in counts.items():
            assert a == source_polynomials(d, n)[n], (d, n, a)
            permutation_cases += 1
    print(f'Fidelity: Theorem 8.2 weighted permutation formula checked by direct enumeration in {permutation_cases} cases (2<=d<=8, 0<=n<=8)')

    comparisons, block_checks, zero_cases = 0, 0, 0
    for d in range(3, 13):
        polys = source_polynomials(d, 46)
        for n in range(1, 25):
            k, s = divmod(n, d)
            base = None
            for x in (0, 1, 2, 3):
                values = [evaluate(a, x) for a in polys]
                value = det([[values[i + j] for j in range(n)] for i in range(n)])
                assert (value == 0) == (2 <= s <= d - 1), (d, n, x, value)
                if base is None:
                    base = value
                assert value == base, (d, n, x, value, base)
                comparisons += 1
            if s in (0, 1):
                b0 = det([[moment(d, i + j) for j in range(k + (s == 1))]
                          for i in range(k + (s == 1))])
                b1 = det([[moment(d, i + j + 1) for j in range(k)] for i in range(k)])
                assert b0 > 0 and b1 > 0, (d, n)
                assert base == (-1) ** (k * ((d - 1) // 2)) * b0 * b1 ** (d - 1), (d, n)
                block_checks += 1
            else:
                zero_cases += 1
    print(f'Conjecture 9.5(3): {comparisons} exact determinant checks OK (3<=d<=12, 1<=n<=24, x=0,1,2,3)')
    print(f'Residue factorization: {block_checks} nonzero cases with positive factors; {zero_cases} zero cases OK')

    # For fixed d,n the polynomial determinant has degree at most n(n-1):
    # every Leibniz product has that bound. Evaluation at degree+1 distinct
    # integers therefore certifies constancy for each of these finite pairs.
    interpolated_pairs, interpolation_values = 0, 0
    for d in range(3, 9):
        polys = source_polynomials(d, 14)
        for n in range(1, 9):
            expected = det([[polys[i + j][0] for j in range(n)] for i in range(n)])
            for x in range(n * (n - 1) + 1):
                values = [evaluate(a, x) for a in polys]
                assert det([[values[i + j] for j in range(n)] for i in range(n)]) == expected
                interpolation_values += 1
            interpolated_pairs += 1
    print(f'Exact polynomial constancy: {interpolated_pairs} finite (d,n) pairs certified by {interpolation_values} integer evaluations (3<=d<=8, 1<=n<=8)')

    shifted = 0
    for d in range(3, 13):
        for e in (0, 1):
            for j in range(21):
                rhs = moment(d, e) * d ** ((d - 1) * j) * prod(
                    (rising(Fraction(e) + Fraction(r, d), j) for r in range(1, d)),
                    start=Fraction(1))
                assert rhs == moment(d, j + e), (d, e, j)
                shifted += 1
    print(f'Shifted moment/Pochhammer identity: {shifted} exact rational checks OK')

    cholesky = 0
    for a in (Fraction(1, 7), Fraction(1, 3), Fraction(1, 2), Fraction(1), Fraction(3, 2)):
        for i in range(9):
            for j in range(9):
                rhs = sum((Fraction(comb(i, q) * comb(j, q) * factorial(q))
                           * rising(a, i) * rising(a, j) / rising(a, q)
                           for q in range(min(i, j) + 1)), start=Fraction(0))
                assert rising(a, i + j) == rhs, (a, i, j)
                cholesky += 1
    print(f'Finite positive LDL identity: {cholesky} exact rational entry checks OK')
    print('Representative determinants for d=3, x=0, n=1..7:',
          [det([[source_polynomials(3, 12)[i + j][0] for j in range(n)] for i in range(n)]) for n in range(1, 8)])
    print('Finite checks only; no universal or Lean proof is claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

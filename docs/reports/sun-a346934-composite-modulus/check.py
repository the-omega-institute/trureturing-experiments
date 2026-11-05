#!/usr/bin/env python3
"""Exact finite checks for OEIS A346934 Conjecture 1's composite clause.

Sources: https://oeis.org/A346934 (revision 29) and Sun,
Arithmetic properties of some permanents, arXiv:2108.07723v7.
These finite checks do not prove the universally quantified conjecture.
"""

from itertools import permutations, product
from math import factorial, isqrt, prod
from random import Random


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def elementary(values):
    """Coefficients of product(1 + value*t), using integer arithmetic."""
    coefficients = [1]
    for value in values:
        coefficients.append(0)
        for r in range(len(coefficients) - 1, 0, -1):
            coefficients[r] += value * coefficients[r - 1]
    return coefficients


def expansion(x, y):
    assert len(x) == len(y)
    n = len(x)
    ex, ey = elementary(x), elementary(y)
    return sum(factorial(r) * factorial(n - r) * ex[r] * ey[n - r]
               for r in range(n + 1))


def permutation_permanent(matrix):
    n = len(matrix)
    return sum(prod(matrix[i][sigma[i]] for i in range(n))
               for sigma in permutations(range(n)))


def subset_permanent(matrix):
    """Independent row-by-row injection recurrence; no rank assumption."""
    n = len(matrix)
    table = [0] * (1 << n)
    table[0] = 1
    for mask in range(1, 1 << n):
        row = mask.bit_count() - 1
        todo = mask
        while todo:
            bit = todo & -todo
            column = bit.bit_length() - 1
            table[mask] += matrix[row][column] * table[mask ^ bit]
            todo ^= bit
    return table[-1]


def additive_matrix(x, y):
    return [[u + v for v in y] for u in x]


def difference_matrix(n):
    return [[i - j for j in range(1, n + 1)] for i in range(1, n + 1)]


def main():
    # Fidelity to the exact OEIS source, including its displayed a(1) example.
    source_values = [
        -1, 52, -18660, 24446016, -85000104000, 647188836814080,
        -9486416237249952000, 244072502056661870592000,
        -10282514440038927957603532800,
        671904022157076034864609763328000,
        -65203712913305114275839483698454528000,
    ]
    assert permutation_permanent([[0, -1], [1, 0]]) == -1
    for n, expected in enumerate(source_values, 1):
        x = list(range(1, 2 * n + 1))
        assert expansion(x, [-v for v in x]) == expected
    print('FIDELITY OEIS revision 29: a(1) example and all 11 displayed terms agree')

    # Independently evaluate actual matrices, including zero and odd orders.
    for size in range(19):
        matrix = difference_matrix(size)
        x = list(range(1, size + 1))
        actual = subset_permanent(matrix)
        assert actual == expansion(x, [-v for v in x])
        if size <= 8:
            assert actual == permutation_permanent(matrix)
        if size % 2:
            assert actual == 0
    print('INDEPENDENT: subset permanents agree with expansion at every order 0..18')
    print('INDEPENDENT: direct permutation sums agree at every order 0..8')

    # Fidelity to the source's proved prime-modulus theorem, not the target.
    primes = [p for p in range(3, 20, 2) if prime(p)]
    for p in primes:
        value = subset_permanent(difference_matrix(p - 1))
        assert value % p == 2
    print('FIDELITY prime-modulus theorem: primes', primes, 'all give residue 2')

    # Exhaustive signed and degenerate inputs test the general rank-two identity.
    exhaustive = 0
    for size in range(5):
        vectors = list(product((-1, 0, 1), repeat=size))
        for x in vectors:
            for y in vectors:
                assert permutation_permanent(additive_matrix(x, y)) == expansion(x, y)
                exhaustive += 1
    random = Random(346934)
    independent = 0
    for size in range(5, 11):
        for _ in range(12):
            x = [random.randrange(-7, 8) for _ in range(size)]
            y = [random.randrange(-7, 8) for _ in range(size)]
            assert subset_permanent(additive_matrix(x, y)) == expansion(x, y)
            independent += 1
    print('RANK_TWO exhaustive vectors in {-1,0,1}, orders 0..4:', exhaustive)
    print('RANK_TWO independent subset checks, orders 5..10, seed 346934:', independent)

    factorials = [factorial(r) for r in range(1000)]
    composites = [q for q in range(3, 1000, 2) if not prime(q)]
    split_checks = 0
    for q in composites:
        for r in range(q):
            assert (factorials[r] * factorials[q - 1 - r]) % q == 0
            split_checks += 1
    print('FACTORIAL_SPLIT odd composite q<=999:', len(composites),
          'moduli;', split_checks, 'splits; all residues zero')
    print('EXCEPTION q=9 splits r=3,4,5:',
          [factorial(r) * factorial(8 - r) for r in (3, 4, 5)])

    target_checks = 0
    for n in range(1, 101):
        q = 2 * n + 1
        x = list(range(1, 2 * n + 1))
        value = expansion(x, [-v for v in x])
        if not prime(q):
            assert value % q == 0
            target_checks += 1
    first = subset_permanent(difference_matrix(8))
    assert first == 9 * 2716224
    print('TARGET every composite q=2n+1 with 1<=n<=100:', target_checks,
          'cases; all residues zero')
    print('TARGET first composite q=9: a(4)=', first, '=9*2716224')

    general_checks = 0
    for q in (9, 15):
        for _ in range(8):
            x = [random.randrange(-9, 10) for _ in range(q - 1)]
            y = [random.randrange(-9, 10) for _ in range(q - 1)]
            value = subset_permanent(additive_matrix(x, y))
            assert value == expansion(x, y) and value % q == 0
            general_checks += 1
    print('GENERAL additive matrices: q=9,15, seed 346934:', general_checks,
          'independent integer-permanent checks')
    print('SCOPE: sign sentence and Conjecture 2 are not nominated or settled here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

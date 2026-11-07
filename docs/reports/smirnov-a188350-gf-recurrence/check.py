#!/usr/bin/env python3
"""Exact finite checks for the two conjectural formulas in OEIS A188350.

The sequence is computed from native open words, independently of either formula.
Only Python's standard library is required. No floating point is used.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path
from math import comb

LETTERS = (-1, 0, 1)
LIMIT = 200


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def convolution(p, q, limit):
    result = [0] * (limit + 1)
    for i, x in enumerate(p[:limit + 1]):
        if x:
            for j, y in enumerate(q[:limit + 1 - i]):
                if y:
                    result[i + j] += x * y
    return result


def binomial_series(q, exponent, limit):
    """Truncate (1+q)^exponent by its formal generalized-binomial expansion."""
    require(q[0] == 0, 'The substituted series must have zero constant term')
    power = [1] + [0] * limit
    result = [Fraction(0) for _ in range(limit + 1)]
    factor = Fraction(1)
    for k in range(limit + 1):
        for i, coefficient in enumerate(power):
            result[i] += factor * coefficient
        factor *= (exponent - k) / (k + 1)
        power = convolution(power, q, limit)
    return result


def native_distributions(limit):
    """dp[(last, energy)] counts actual words of n+1 letters."""
    dp = {(letter, 0): 1 for letter in LETTERS}
    rows = []
    for n in range(limit + 1):
        require(sum(dp.values()) == 3 ** (n + 1), f'word mass at n={n}')
        require(all(dp.get((-last, e), 0) == count
                    for (last, e), count in dp.items()), f'sign reversal at n={n}')
        rows.append(dict(dp))
        following = defaultdict(int)
        for (last, e), count in dp.items():
            for appended in LETTERS:
                following[appended, e + last * appended] += count
        dp = following
    return rows


def collapsed_distributions(limit):
    """Two endpoint states: u for last +1, v for last 0; -1 is sign reversal."""
    u, v = {0: 1}, {0: 1}
    rows = []
    for _ in range(limit + 1):
        rows.append((dict(u), dict(v)))
        next_u, next_v = defaultdict(int), defaultdict(int)
        for e, count in u.items():
            next_u[e + 1] += count
            next_u[e - 1] += count
            next_v[e] += 2 * count
        for e, count in v.items():
            next_u[e] += count
            next_v[e] += count
        u, v = next_u, next_v
    return rows


def main():
    require(all(isinstance(s, int) for s in LETTERS), 'Integer alphabet')
    native = native_distributions(LIMIT)
    counts = [sum(row.get((last, 0), 0) for last in LETTERS) for row in native]
    initial = [3, 5, 15, 29, 91, 197, 627, 1477, 4671, 11693,
               36487, 95405, 293695, 792525]
    require(counts[:len(initial)] == initial, 'Native initial values')
    for n in range(9):
        brute = defaultdict(int)
        for word in product(LETTERS, repeat=n + 1):
            energy = sum(word[i] * word[i + 1] for i in range(n))
            brute[word[-1], energy] += 1
        require(dict(brute) == native[n], f'Brute-force distribution at n={n}')
    print('Native open-word DP: edges 0..200; brute-force distributions agree at 0..8')
    print('Initial counts a_star(0)..a_star(13):', ','.join(map(str, initial)))

    source = {}
    for line in Path(__file__).with_name('b188350.txt').read_text().splitlines():
        fields = line.split()
        if fields and not fields[0].startswith('#'):
            n, value = map(int, fields)
            source[n] = value
    require(set(source) == set(range(1, LIMIT + 1)), 'Source b-file range')
    require(all(counts[n] == value for n, value in source.items()), 'OEIS fidelity')
    print('Source fidelity: all 200 entries of Hardin OEIS b188350.txt agree (n=1..200)')

    collapsed = collapsed_distributions(LIMIT)
    for n, (u, v) in enumerate(collapsed):
        for e in range(-n, n + 1):
            require(u.get(e, 0) == native[n].get((1, e), 0), f'u at {(n,e)}')
            require(v.get(e, 0) == native[n].get((0, e), 0), f'v at {(n,e)}')
    print('Two-state endpoint transfer agrees with every native energy coefficient at 0..200')

    # D=(1-3t)(1+t-4t^2)=1-2t-7t^2+12t^3.
    d = [1, -2, -7, 12]
    require(convolution([1, -3], [1, 1, -4], 3) == d, 'Discriminant factors')
    A, B = [1, -1, -2], [0, 1, -1]
    aa, bb = convolution(A, A, 4), convolution(B, B, 4)
    require([x - 4*y for x, y in zip(aa, bb)] == d + [0], 'Kernel discriminant')
    h = binomial_series([0, -2, -7, 12], Fraction(-1, 2), LIMIT)
    s = binomial_series([0, -2, -7, 12], Fraction(1, 2), LIMIT)
    identity = [1] + [0] * LIMIT
    require(convolution(s, s, LIMIT) == d + [0] * (LIMIT - 3), 'Square root branch')
    require(convolution(h, s, LIMIT) == identity, 'Reciprocal square root')
    require(convolution(d, convolution(h, h, LIMIT), LIMIT) == identity, 'D H^2=1')
    numerator = [2*x for x in h]
    numerator[0] -= 2
    numerator[1] += 3
    gf, total = [], Fraction(0)
    for coefficient in numerator:
        total += coefficient
        gf.append(total)
    native_gf = [0] + counts[1:]
    require(gf == native_gf, 'Published generating function')
    lhs = convolution(convolution([1, -1], s, LIMIT), native_gf, LIMIT)
    rhs = convolution([-2, 3], s, LIMIT)
    rhs[0] += 2
    require(lhs == rhs, 'Denominator-cleared binding GF equality')
    print('Formal GF and normalized root: exact generalized-binomial coefficients agree through degree 200')

    for n in range(5, LIMIT + 1):
        residual = (n*counts[n] - (3*n-1)*counts[n-1] - (5*n-6)*counts[n-2]
                    + (19*n-25)*counts[n-3] - 6*(2*n-3)*counts[n-4])
        require(residual == 0, f'Published recurrence at n={n}')
    print('Published order-four recurrence: residual zero for every n=5..200')
    require(all(comb(2*j, j) == sum(1 for word in product((-1, 1), repeat=2*j)
                                  if sum(word) == 0) for j in range(7)),
            'Central-binomial constant-term check')
    print('Central-binomial constant-term counts agree for j=0..6')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

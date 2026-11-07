#!/usr/bin/env python3
"""Exact finite checks for Provvedi's prime-period conjecture, OEIS A249938.

Only Python integer arithmetic is used. These finite checks are not a proof
of the quantified conjecture. Run with Python 3.8 or later.
"""
from math import comb, factorial, isqrt


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def primes_up_to(limit):
    return [p for p in range(2, limit + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def fubini_by_stirling(limit):
    """S(N,k)=k*S(N-1,k)+S(N-1,k-1); F_N=sum k! S(N,k)."""
    row = [1]
    facts = [factorial(k) for k in range(limit + 1)]
    result = [1]
    for n in range(1, limit + 1):
        row = [0] + [k * (row[k] if k < len(row) else 0) + row[k - 1]
                     for k in range(1, n + 1)]
        result.append(sum(facts[k] * row[k] for k in range(n + 1)))
    return result


def fubini_by_first_block(limit):
    """Independent ordered-partition recurrence, with F_0=1."""
    result = [1]
    for n in range(1, limit + 1):
        result.append(sum(comb(n, k) * result[k] for k in range(n)))
    return result


def multiplicative_order(base, p):
    x = 1
    for d in range(1, p):
        x = x * base % p
        if x == 1:
            return d
    raise AssertionError('No multiplicative order')


def least_cyclic_period(word):
    return next(t for t in range(1, len(word) + 1)
                if all(word[i] == word[(i + t) % len(word)]
                       for i in range(len(word))))


def main():
    primes = primes_up_to(199)
    max_m = (primes[-1] - 1) // 2
    max_n = 4 * max_m + 17
    f = fubini_by_stirling(2 * max_n)
    a = f[::2]
    require(f[:81] == fubini_by_first_block(80), 'Independent recurrence')
    require(a[:5] == [1, 3, 75, 4683, 545835], 'OEIS initial terms')
    example = [a[n] % 23 for n in range(1, 12)]
    require(example == [3, 6, 14, 22, 5, 3, 10, 2, 4, 5, 0],
            'Source p=23 example')
    print('FIDELITY a(0..4) =', a[:5])
    print('FIDELITY p=23, n=1..11 =', example)
    print('Independent ordered-partition recurrence: F_0 through F_80 agree')
    print('p expected least square_count surviving zero_weights ord_p(4)')
    for p in primes:
        expected = 1 if p == 2 else (p - 1) // 2
        squares = {j * j % p for j in range(1, p)}
        require(len(squares) == expected, ('square count', p))
        word = [a[n] % p for n in range(1, expected + 1)]
        least = least_cyclic_period(word)
        require(least == expected, ('least cyclic period', p, least))
        for n in range(1, max_n - expected + 1):
            require(a[n + expected] % p == a[n] % p, ('period check', p, n))
        for onset in [1, 2, 17, 2 * expected + 1]:
            for t in range(1, expected):
                require(any((a[n + t] - a[n]) % p != 0
                            for n in range(onset, onset + expected)),
                        ('shorter tail shift not excluded', p, onset, t))
        if p == 2:
            require(all(x % 2 == 1 for x in a[1:]), 'p=2 positive terms')
            print(p, expected, least, len(squares), '-', '-', '-')
            continue
        require(a[0] % p != a[expected] % p, ('onset zero exclusion', p))
        coeffs = [sum((-1) ** (k - j) * comb(k, j)
                      for k in range(j, p)) % p for j in range(p)]
        require(all(coeffs[j] == pow(2, -j, p) for j in range(1, p)),
                ('finite-sum coefficients', p))
        for n in range(1, 2 * max_n + 1):
            spectral = sum(pow(2, -j, p) * pow(j, n, p)
                           for j in range(1, p)) % p
            require(spectral == f[n] % p, ('full exponential formula', p, n))
        m = expected
        lambdas = [r * r % p for r in range(1, m + 1)]
        weights = [(pow(2, -r, p) + pow(2, r - 1, p)) % p
                   for r in range(1, m + 1)]
        require(set(lambdas) == squares and len(set(lambdas)) == m,
                ('square enumeration', p))
        require(all((w == 0) == (pow(4, r, p) == (-2) % p)
                    for r, w in enumerate(weights, 1)), ('zero weights', p))
        for n in range(1, max_n + 1):
            paired = sum(w * pow(lam, n, p)
                         for w, lam in zip(weights, lambdas)) % p
            require(paired == a[n] % p, ('paired formula', p, n))
        surviving = sum(w != 0 for w in weights)
        zeros = m - surviving
        d = multiplicative_order(4, p)
        require(m % d == 0, ('order divides square group', p))
        require(zeros in (0, m // d), ('zero-weight count', p))
        if p >= 7:
            require(d >= 3 and 3 * zeros <= m, ('support bound', p))
        if p >= 5:
            require(2 * surviving > m, ('strictly more than half survive', p))
        if p == 3:
            require(surviving == 0 and all(x % 3 == 0 for x in a[1:]),
                    'p=3 zero sequence')
        determinant = 1
        for lam in lambdas:
            determinant = determinant * lam % p
        for i, lam in enumerate(lambdas):
            for other in lambdas[i + 1:]:
                determinant = determinant * (other - lam) % p
        require(determinant != 0, ('Vandermonde determinant', p))
        print(p, expected, least, len(squares), surviving, zeros, d)
    print('Checked', len(primes), 'primes from 2 through 199; positive indices through', max_n)
    print('All smaller shifts excluded in four finite tail windows for each prime')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

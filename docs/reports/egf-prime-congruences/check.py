#!/usr/bin/env python3
"""Exact checks for both conjectures in Melika Tebni's OEIS A350759.

Source: https://oeis.org/search?q=id:A350759&fmt=text, revision 17.
DATA and the recurrence are attributed to the OEIS contributors; see README.md.
Only integer and fractions.Fraction arithmetic is used as evidence.
"""
from fractions import Fraction
from math import comb, factorial, isqrt

MAX_N = 220
EXPANSION_N = 60
SOURCE_DATA = (
    1, 0, -1, 1, 1, -4, 1, 29, -167, 1000, -7989, 75857,
    -794639, 9058180, -111944923, 1492748581, -21369667087,
    326932765840, -5323818187817, 91947960224097,
    -1678914212753599, 32317295442288844,
    -654084630476955479, 13886774070229667213,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def source_recurrence(limit):
    a = [1]
    for n in range(1, limit + 1):
        a.append(sum(factorial(n - 2 - k) * comb(n - 1, k)
                     * (-1) ** (n - 1 - k) * a[k]
                     for k in range(n - 1)))
    return a


def formal_exponential(limit):
    # G = X - (1+X) log(1+X): g_0 = g_1 = 0,
    # g_m = (-1)^(m+1)/(m(m-1)) for m >= 2.
    # Differentiate exp(G): n*f_n = sum_{k=1}^n k*g_k*f_(n-k).
    g = [Fraction(0), Fraction(0)] + [
        Fraction((-1) ** (m + 1), m * (m - 1))
        for m in range(2, limit + 1)
    ]
    f = [Fraction(1)]
    for n in range(1, limit + 1):
        f.append(sum((k * g[k] * f[n - k]
                      for k in range(1, n + 1)), Fraction(0)) / n)
    scaled = [factorial(n) * f[n] for n in range(limit + 1)]
    require(all(v.denominator == 1 for v in scaled), 'EGF integrality')
    return [v.numerator for v in scaled]


def stirling_transform(limit):
    # A345652 recurrence from that entry's FORMULA.
    b = [1]
    for n in range(1, limit + 1):
        b.append(sum((n - 1) * comb(n - 2, k)
                     * (-1) ** (n - 1 - k) * b[k]
                     for k in range(n - 1)))
    # Signed first-kind Stirling numbers: s(n+1,k)=s(n,k-1)-n*s(n,k).
    row = [1]
    a = [1]
    for n in range(1, limit + 1):
        prev = row
        row = [(prev[k - 1] if k else 0)
               - (n - 1) * (prev[k] if k < n else 0)
               for k in range(n + 1)]
        a.append(sum((-1) ** k * b[k] * row[k] for k in range(n + 1)))
    return a


def rising_expansion(limit):
    # P_k(X) = product_{r=1}^k (X+r), with coefficient list polys[k].
    polys = [[1]]
    for k in range(1, limit + 1):
        old = polys[-1]
        new = [0] * (len(old) + 1)
        for degree, value in enumerate(old):
            new[degree] += k * value
            new[degree + 1] += value
        polys.append(new)
    values = []
    count = 0
    for n in range(limit + 1):
        total = Fraction(0)
        for k in range(n + 1):
            require(polys[k][0] == factorial(k), f'P_{k}(0)')
            for ell in range(min(k, n - k) + 1):
                j = n - k - ell
                term = Fraction((-1) ** k * factorial(n) * polys[k][ell],
                                factorial(k) * factorial(j))
                require(term.denominator == 1, f'summand integrality {n,k,ell}')
                if ell:
                    ratio = Fraction(factorial(n), factorial(k) * factorial(j))
                    factorization = (n * factorial(ell - 1) * comb(n - 1, k)
                                     * comb(n - 1 - k, j))
                    require(ratio == factorization, f'factorial identity {n,k,ell}')
                if n and (k or ell):
                    require(term.numerator % n == 0, f'summand divisibility {n,k,ell}')
                    count += 1
                total += term
        require(total.denominator == 1, f'expansion integrality n={n}')
        values.append(total.numerator)
    return values, count


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def main():
    a = source_recurrence(MAX_N)
    require(tuple(a[:len(SOURCE_DATA)]) == SOURCE_DATA, 'OEIS DATA fidelity')
    require(a[9] == 1000, 'source example a(9)')
    print('FIDELITY: all 24 OEIS DATA terms matched; source example a(9)=1000')
    print('DATA:', ','.join(map(str, a[:24])))
    require(a == formal_exponential(MAX_N), 'published EGF vs recurrence')
    print(f'EGF: exact rational formal exponential matched recurrence for n=0..{MAX_N}')
    require(a[:EXPANSION_N + 1] == stirling_transform(EXPANSION_N),
            'signed Stirling-transform definition')
    print(f'NAME: signed first-kind Stirling transform matched for n=0..{EXPANSION_N}')
    expanded, count = rising_expansion(EXPANSION_N)
    require(a[:EXPANSION_N + 1] == expanded, 'rising-factorial expansion')
    print(f'EXPANSION: matched for n=0..{EXPANSION_N}; {count} nonexceptional summands divisible by n')
    for n in range(1, MAX_N + 1):
        require((a[n] - 1) % n == 0, f'all-index congruence n={n}')
    print(f'UNIFORM: a(n)=1 mod n for all n=1..{MAX_N}')
    primes = [p for p in range(2, 212) if is_prime(p)]
    for p in primes:
        require((a[p] + a[p + 1] - 2) % p == 0, f'first conjecture p={p}')
        require((a[p + 1] - 1) % (p * (p + 1)) == 0, f'second conjecture p={p}')
        require((a[p + 1] - 1) % p == 0, f'successor congruence p={p}')
    require(a[2] == -1 and a[3] == 1, 'p=2 boundary')
    print(f'TARGETS: both complete assertions checked on all {len(primes)} primes p<=211')
    print('PRIMES:', ','.join(map(str, primes)))
    print('Finite checks do not establish the universal claims.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

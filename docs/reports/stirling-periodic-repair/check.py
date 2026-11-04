#!/usr/bin/env python3
"""Exact finite checks for Miska--Ward, arXiv:2102.07561v2.

These computations test the proposed general argument. Truncated denominator
lcms are empirical lower divisors of Fail, not certified infinite repair factors.
Only Python 3 standard-library integers and Fraction arithmetic are used.
"""
from fractions import Fraction
from math import factorial, isqrt, lcm, prod

N = 3000
K = 40
SOURCE_TABLE = [
    1, 1, 2, 6, 12, 60, 30, 210, 840, 2520, 1260, 13860, 13860,
    180180, 90090, 30030, 240240, 4084080, 6126120, 116396280,
    58198140, 58198140, 29099070, 669278610, 892371480, 4461857400,
    2230928700, 20078358300, 20078358300, 582272390700, 291136195350,
    9025222055850, 144403552893600, 48134517631200, 24067258815600,
    24067258815600, 36100888223400, 1335732864265800, 333933216066450,
    333933216066450,
]


def primes_upto(n):
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, isqrt(p) + 1))]


def valuation(n, p):
    assert n > 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def floor_log(p, n):
    t, power = 0, 1
    while power * p <= n:
        t, power = t + 1, power * p
    return t


def mobius_table(n):
    mu = [1] * (n + 1)
    mu[0] = 0
    for p in primes_upto(n):
        for m in range(p, n + 1, p):
            mu[m] *= -1
        for m in range(p * p, n + 1, p * p):
            mu[m] = 0
    return mu


def shifted(k, nmax, modulus=None):
    row = [1] + [0] * k
    values = [0] * (nmax + 1)
    for r in range(1, nmax + k):
        for j in range(k, 0, -1):
            x = j * row[j] + row[j - 1]
            row[j] = x if modulus is None else x % modulus
        row[0] = 0
        if r >= k:
            values[r - k + 1] = row[k]
    return values


def mul_mod_f(a, b, f, modulus):
    """Multiply coefficient lists modulo a monic f and an integer modulus."""
    q = len(f) - 1
    c = [0] * (2 * q - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % modulus
    for i in range(2 * q - 2, q - 1, -1):
        x = c[i]
        for j in range(q):
            c[i - q + j] = (c[i - q + j] - x * f[j]) % modulus
    return c[:q]


def power_x_mod(e, f, modulus):
    q = len(f) - 1
    a, x = [1] + [0] * (q - 1), [0, 1] + [0] * (q - 2)
    while e:
        if e % 2:
            a = mul_mod_f(a, x, f, modulus)
        x = mul_mod_f(x, x, f, modulus)
        e //= 2
    return a


def root_polynomial(k):
    f = [1]
    for j in range(k + 1):
        g = [0] * (len(f) + 1)
        for r, c in enumerate(f):
            g[r] -= j * c
            g[r + 1] += c
        f = g
    return f


def main():
    mu = mobius_table(N)
    ps = primes_upto(K)
    sampled = []
    source_cases = 0
    for k in range(1, K + 1):
        a = shifted(k, N)
        if k == 3:
            assert a[1:5] == [1, 6, 25, 90]  # Source (2).
        sums = [0] * (N + 1)
        for d in range(1, N + 1):
            for h in range(1, N // d + 1):
                sums[d * h] += mu[h] * a[d]
        empirical = 1
        for n in range(1, N + 1):
            # Reproduce source Theorem 3 in this finite window.
            assert sums[n] >= 0
            denominator = Fraction(sums[n], n).denominator
            assert factorial(k - 1) % denominator == 0
            empirical = lcm(empirical, denominator)
            source_cases += 1
        assert empirical == SOURCE_TABLE[k - 1], (k, empirical)
        sampled.append(empirical)
        upper = lcm(*range(1, k)) if k >= 2 else 1
        radical = prod(p for p in ps if p < k)
        assert upper % empirical == 0
        assert empirical % radical == 0
        for p in ps:
            if p < k:
                assert 1 <= valuation(empirical, p) <= floor_log(p, k - 1)
                if k <= p * p:
                    assert valuation(empirical, p) == 1
            else:
                assert empirical % p != 0
    print(f'FIDELITY: source S_3 initial terms and all {K} Table 1 entries reproduced at n <= {N}')
    print(f'FIDELITY: nonnegative Mobius sums and factorial denominator bound: {source_cases} cases')
    print('EMPIRICAL_LCMS k=1..40: ' + ','.join(map(str, sampled)))
    print(f'TARGET: radical lower divisor and lcm upper divisor for the {K} truncated repair factors')

    sharp = 0
    for p in ps:
        j, power = 1, p
        while power + 1 <= K:
            k = power + 1
            length = p ** (j - 1) * (p - 1)
            a = shifted(k, N, p)
            assert all(a[n] == ((n - 1) // length + 1) % p
                       for n in range(1, N + 1))
            assert (a[power] - a[power // p]) % p == 1
            assert valuation(sampled[k - 1], p) == j
            sharp += 1
            j, power = j + 1, power * p
    print(f'SHARPNESS: {sharp} prime-power pairs with p^j+1 <= {K}; full residue formula at n <= {N}')

    polynomial_cases = 0
    for k in range(3, 19):
        f = root_polynomial(k)
        # Inclusion-exclusion annihilation is independently tested using exact Stirling recurrence.
        unshifted = [0] * (3 * k + 1)
        row = [1] + [0] * k
        for r in range(3 * k + 1):
            unshifted[r] = row[k]
            for z in range(k, 0, -1):
                row[z] = z * row[z] + row[z - 1]
            row[0] = 0
        for r in range(k):
            assert sum(c * unshifted[i + r] for i, c in enumerate(f)) == 0
        for p in ps:
            if p >= k:
                continue
            t = floor_log(p, k - 1)
            for a in range(1, 4):
                modulus = p ** a
                q = k // p + 1
                for n in range(1, 4):
                    low = k - 1 + n * p ** (a + t - 1)
                    high = k - 1 + n * p ** (a + t)
                    assert low >= a * q
                    assert power_x_mod(low, f, modulus) == power_x_mod(high, f, modulus)
                    polynomial_cases += 1
    print(f'POLYNOMIAL: {polynomial_cases} exact quotient congruences, 3 <= k <= 18, 1 <= a,n <= 3')
    print('SCOPE: finite arithmetic checks only; no universal claim or infinite lcm is certified here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for Paul D. Hanna's OEIS A156216 second comment.

Source: https://oeis.org/A156216/internal, revision #3.
This checks finite truncations; it does not prove the universal assertion.
Only Python's standard library is required.
"""
from fractions import Fraction
from itertools import product
from math import comb, factorial
from random import Random

SOURCE_TERMS = (
    1, 1, 5, 26, 634, 32928, 5704263, 2470113915, 2978904483553,
    9401949327631932, 79268874871208384494, 1762019469678472912173354,
    103537245443913551792800303420, 16030602885085486700462431379649694,
)


def mul(a, b, degree):
    return [sum(a[k] * b[n-k] for k in range(n+1))
            for n in range(degree+1)]


def divisors(n):
    return [d for d in range(1, n+1) if n % d == 0]


def mobius(n):
    sign, p = 1, 2
    while p*p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
        p += 1
    return -sign if n > 1 else sign


def ghosts(a):
    """Coefficients of x F'/F, using F(0)=1 and integer operations."""
    assert a[0] == 1
    c = [0] * len(a)
    for n in range(1, len(a)):
        c[n] = n*a[n] - sum(c[k]*a[n-k] for k in range(1, n))
    return c


def exp_recurrence(c):
    """Source recurrence: n a_n = sum_{k=1}^n c_k a_{n-k}."""
    a = [Fraction(1)]
    for n in range(1, len(c)):
        a.append(sum(c[k]*a[n-k] for k in range(1, n+1)) / n)
    return a


def exp_definition(c):
    """Independent finite formal exponential: sum_{j=0}^N H^j/j!."""
    degree = len(c)-1
    h = [Fraction(0)] + [Fraction(c[n], n) for n in range(1, degree+1)]
    term = [Fraction(1)] + [Fraction(0)]*degree
    ans = term.copy()
    for j in range(1, degree+1):
        term = mul(term, h, degree)
        ans = [a + b/factorial(j) for a, b in zip(ans, term)]
    return ans


def ledger(c):
    b = [0] * len(c)
    for n in range(1, len(c)):
        numerator = sum(mobius(n//d)*c[d] for d in divisors(n))
        b[n], remainder = divmod(numerator, n)
        assert remainder == 0, ('nonintegral ledger', n, numerator)
    return b


def euler_product(b):
    degree = len(b)-1
    a = [1] + [0]*degree
    for d in range(1, degree+1):
        exponent = b[d]
        if exponent == 0:
            continue
        factor = [0]*(degree+1)
        factor[0] = 1
        for j in range(1, degree//d+1):
            if exponent > 0:
                factor[d*j] = comb(exponent+j-1, j)
            elif j <= -exponent:
                factor[d*j] = (-1)**j * comb(-exponent, j)
        a = mul(a, factor, degree)
    return a


def check_local(c):
    for n in range(2, len(c)):
        for p in divisors(n):
            if p < 2 or any(p % k == 0 for k in range(2, p)):
                continue
            m, modulus = n, 1
            while m % p == 0:
                m //= p
                modulus *= p
            assert (c[n]-c[n//p]) % modulus == 0, (n, p, modulus)


def check(a, direct_exp=False):
    c = ghosts(a)
    assert exp_recurrence(c) == a
    assert euler_product(ledger(c)) == a
    check_local(c)
    d = [0] + [c[n]**n for n in range(1, len(c))]
    check_local(d)
    e = ledger(d)
    for n in range(1, len(c)):
        assert d[n] == sum(k*e[k] for k in divisors(n))
    g = exp_recurrence(d)
    assert all(x.denominator == 1 for x in g)
    assert g == euler_product(e)
    if direct_exp:
        assert exp_definition(c) == a
        assert exp_definition(d) == g
    return c, g


def main():
    degree = 40
    lucas = [2, 1]
    for n in range(2, degree+1):
        lucas.append(lucas[-1]+lucas[-2])
    fib_gf = [1, 1]
    for n in range(2, degree+1):
        fib_gf.append(fib_gf[-1]+fib_gf[-2])
    c, g = check(fib_gf)
    assert c[1:] == lucas[1:]
    assert tuple(g[:len(SOURCE_TERMS)]) == SOURCE_TERMS
    print('FIDELITY: all 14 OEIS A156216 displayed terms, degrees 0..13, agree')
    print('FIDELITY: Lucas ghosts of 1/(1-x-x^2), degrees 1..40, agree')
    print('A156216 prefix:', ','.join(str(int(x)) for x in g[:14]))

    count = 0
    for coefficients in product(range(-2, 3), repeat=5):
        check([1, *coefficients] + [0]*19)
        count += 1
    print(f'POLYNOMIAL BOX: {count} inputs, a_1..a_5 in -2..2, degrees 0..24: OK')

    count = 0
    for entries in product(range(-2, 3), repeat=4):
        b = [0, *entries] + [0]*36
        check(euler_product(b))
        count += 1
    print(f'SIGNED LEDGER BOX: {count} inputs, b_1..b_4 in -2..2, degrees 0..40: OK')

    rng = Random(156216)
    extra = [[1]+[0]*24, [1]+[-1]*24, [1]+[(-1)**n for n in range(1, 25)]]
    extra += [[1]+[rng.randrange(-20, 21) for _ in range(24)] for _ in range(9)]
    for a in extra:
        check(a, direct_exp=True)
    print('INDEPENDENT EXPONENTIAL: 12 signed inputs, exact sum H^j/j!, degrees 0..24: OK')

    for z in range(-30, 31):
        for p in (2, 3, 5, 7, 11):
            for r in range(1, 5):
                modulus = p**r
                assert (pow(z, p**r, modulus)-pow(z, p**(r-1), modulus)) % modulus == 0
    print('FERMAT LIFT: z=-30..30, p=2,3,5,7,11, r=1..4: OK')
    print('No positivity condition; zero/negative ghosts, negative Euler exponents and p=2 included')
    print('Finite exact checks only; no universal proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

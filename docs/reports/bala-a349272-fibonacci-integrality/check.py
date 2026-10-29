#!/usr/bin/env python3
"""Exact finite checks for Bala's OEIS A349272 polynomial conjecture.

Polynomials are tuples of integer coefficients, in increasing degree order.
The only divisions use fractions.Fraction or monic polynomial long division.
These checks are finite evidence; they do not prove the unbounded conjecture.
"""
from fractions import Fraction
from math import prod

ZERO = (0,)
ONE = (1,)
X = (0, 1)


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        c[i] += v
    for i, v in enumerate(b):
        c[i] += v
    return trim(c)


def sub(a, b):
    return add(a, tuple(-v for v in b))


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, v in enumerate(a):
        if v:
            for j, w in enumerate(b):
                c[i + j] += v * w
    return trim(c)


def power(a, n):
    r = ONE
    for _ in range(n):
        r = mul(r, a)
    return r


def product(items):
    r = ONE
    for a in items:
        r = mul(r, a)
    return r


def total(items):
    r = ZERO
    for a in items:
        r = add(r, a)
    return r


def evaluate(a, x):
    r = 0
    for v in reversed(a):
        r = r * x + v
    return r


def monic_divmod(a, b):
    assert b != ZERO and b[-1] == 1
    r = list(a)
    q = [0] * max(1, len(a) - len(b) + 1)
    while trim(r) != ZERO and len(r) >= len(b):
        shift = len(r) - len(b)
        lead = r[-1]
        q[shift] += lead
        for j, v in enumerate(b):
            r[shift + j] -= lead * v
        r = list(trim(r))
    return trim(q), trim(r)


def fibonacci(limit):
    f = [0, 1]
    for _ in range(2, limit + 1):
        f.append(f[-1] + f[-2])
    return f


def fibonacci_polynomials(limit):
    f = [ZERO, ONE]
    for _ in range(2, limit + 1):
        f.append(add(mul(X, f[-1]), f[-2]))
    return f


def determinant_one_sequence(t, limit):
    u = [ZERO, ONE]
    for _ in range(2, limit + 1):
        u.append(sub(mul(t, u[-1]), u[-2]))
    return u


def main():
    # OEIS A349272 revision #13, DATA: indices 0 through 9.
    expected = [
        1, 2, 315, 2471040, 918185538816, 16047302734562299200,
        13178031727820369629763174400,
        508406658175888466343652105865846784000,
        921456090985190879093613420564815806955580862464000,
        78458394721620642094151397745899367347021362840662985785265356800,
    ]
    f = fibonacci(500)
    values = []
    for n, wanted in enumerate(expected):
        terms = [f[2 * k] for k in range(1, 2 * n + 2)]
        ratio = Fraction(prod(terms), sum(terms))
        assert ratio.denominator == 1 and ratio.numerator == wanted
        # Independently reproduce the source's first FORMULA line.
        source_formula = Fraction(f[2 * n + 1], f[2 * n + 2]) * prod(terms[:-1])
        assert source_formula == ratio
        values.append(ratio.numerator)
    print('FIDELITY OEIS A349272 DATA n=0..9:', values)
    print('FIDELITY source m=2 FORMULA: 10/10 exact matches')

    u = determinant_one_sequence(X, 27)
    additions = 0
    for a in range(13):
        for b in range(13):
            assert u[a + b + 1] == sub(mul(u[a + 1], u[b + 1]), mul(u[a], u[b]))
            additions += 1
    for n in range(13):
        assert total(u[1:2 * n + 2]) == mul(u[n + 1], add(u[n + 1], u[n]))
    print(f'IDENTITIES in Z[T]: {additions} additions; 13 odd partial sums')

    fp = fibonacci_polynomials(132)
    cases = 0
    reductions = 0
    max_numerator_degree = 0
    for m in range(2, 13, 2):
        c = fp[m]
        t = add(fp[m + 1], fp[m - 1])
        ut = determinant_one_sequence(t, 11)
        for k in range(12):
            assert fp[m * k] == mul(c, ut[k])
            reductions += 1
        for n in range(6):
            terms = [fp[m * k] for k in range(1, 2 * n + 2)]
            numerator = product(terms)
            denominator = total(terms)
            assert denominator != ZERO
            assert len(denominator) - 1 == m * (2 * n + 1) - 1
            assert denominator[-1] == 1
            # Polynomial nonzero does not mean nonzero at every x.
            assert evaluate(denominator, 0) == 0
            if n == 0:
                h = ONE
            else:
                w = product(ut[k] for k in range(1, 2 * n + 1) if k != n + 1)
                h = mul(mul(power(c, 2 * n), sub(ut[n + 1], ut[n])), w)
            assert mul(denominator, h) == numerator
            q, r = monic_divmod(numerator, denominator)
            assert r == ZERO and q == h
            integer_terms = [f[m * k] for k in range(1, 2 * n + 2)]
            assert Fraction(prod(integer_terms), sum(integer_terms)) == evaluate(h, 1)
            max_numerator_degree = max(max_numerator_degree, len(numerator) - 1)
            cases += 1
    print(f'EVEN-STEP REDUCTION: {reductions} polynomial identities')
    print(f'POLYNOMIAL m=2,4,6,8,10,12; n=0..5: {cases} identities D*H=P')
    print(f'INDEPENDENT monic long division: {cases} integral quotients, zero remainders')
    print(f'DENOMINATORS: {cases} monic nonzero polynomials; all vanish at x=0')
    print(f'MAXIMUM numerator degree checked: {max_numerator_degree}')

    integer_cases = 0
    for m in range(2, 21, 2):
        for n in range(13):
            terms = [f[m * k] for k in range(1, 2 * n + 2)]
            ratio = Fraction(prod(terms), sum(terms))
            assert ratio.denominator == 1
            integer_cases += 1
    print(f'INTEGER x=1; m=2,4,...,20; n=0..12: {integer_cases} integral quotients')
    assert Fraction(f[1] * f[2] * f[3], f[1] + f[2] + f[3]) == Fraction(1, 2)
    print('BOUNDARY excluded odd step m=1,n=1,x=1 gives 1/2')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

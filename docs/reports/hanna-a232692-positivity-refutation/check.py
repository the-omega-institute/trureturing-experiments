#!/usr/bin/env python3
"""Exact checks for the Hanna A232692 positivity conjecture.

Only integer and Fraction arithmetic is used. This is experimental evidence,
not a Lean proof. Source: https://oeis.org/A232692/internal, revision 19;
threshold: https://oeis.org/A047399/internal, published offset 1.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path

DEGREE = 430
PUBLISHED = [1, 1, 3, 24, 213, 3096, 46071, 967608, 20251809,
             555747048, 15004870731, 508165972056, 16810393586733,
             677183788645704, 26523956467895103, 1238567261126084856,
             56056407696184372281, 2976966230117448265128,
             152872356339113679491859, 9098430770913969095416728]


def require(condition, description):
    if not condition:
        raise AssertionError(description)


def derivatives(n, k, degree):
    """G^(r)=exp((1+(k-n)r)y) P_r(y)/(1+ny)^(2r-1), y=log G."""
    a, p = [1, 1], [1]
    for r in range(1, degree):
        nxt = [0] * (len(p) + 1)
        for j, value in enumerate(p):
            nxt[j] += (1 + n + (k - 3*n)*r + n*j) * value
            nxt[j+1] += n * (1 + (k-n)*r) * value
            if j:
                nxt[j-1] += j * value
        p = nxt
        a.append(p[0])
    return a[:degree+1]


def egf_mul(a, b, degree):
    return [sum(comb(r, j)*a[j]*b[r-j] for j in range(r+1))
            for r in range(degree+1)]


def ode_coefficients(degree):
    """Independently solve A^7 A'' - 6 A^6 (A')^2 + 3 (A')^3=0.

    Multiplication is binomial convolution of factorial-normalized coefficients.
    At order r the coefficient of a_(r+2) is exactly 1.
    """
    a = [1, 1]
    powers = [[1] for _ in range(8)]
    d2, d3 = [], []
    for r in range(degree-1):
        if r:
            powers[1].append(a[r])
            for j in range(2, 8):
                powers[j].append(sum(comb(r,i)*powers[j-1][i]*a[r-i]
                                     for i in range(r+1)))
        d2.append(sum(comb(r,i)*a[i+1]*a[r-i+1] for i in range(r+1)))
        d3.append(sum(comb(r,i)*d2[i]*a[r-i+1] for i in range(r+1)))
        rhs = 6*sum(comb(r,i)*powers[6][i]*d2[r-i] for i in range(r+1))-3*d3[r]
        known = sum(comb(r,i)*powers[7][i]*a[r-i+2] for i in range(1,r+1))
        a.append(rhs-known)
    return a


def mul(a, b, degree):
    return [sum((a[j]*b[r-j] for j in range(r+1)), Q(0))
            for r in range(degree+1)]


def power(a, exponent, degree):
    result = [Q(1)] + [Q(0)]*degree
    while exponent:
        if exponent % 2:
            result = mul(result, a, degree)
        exponent //= 2
        if exponent:
            a = mul(a, a, degree)
    return result


def inverse(a, degree):
    require(a[0] != 0, 'inverse constant is nonzero')
    b = [1/a[0]]
    for r in range(1,degree+1):
        b.append(-sum((a[i]*b[r-i] for i in range(1,r+1)),Q(0))/a[0])
    return b


def exponential(f, degree):
    require(f[0] == 0, 'formal exponential has zero-constant argument')
    b = [Q(1)]
    for r in range(1,degree+1):
        b.append(sum((j*f[j]*b[r-j] for j in range(1,r+1)),Q(0))/r)
    return b


def integral(f, degree):
    return [Q(0)] + [f[r-1]/r for r in range(1,degree+1)]


def fixed_point(n, k, degree):
    """Literal source iteration; degree j is stable after j iterations."""
    g = [Q(1)]+[Q(0)]*degree
    for _ in range(degree):
        g = exponential(mul(power(inverse(g,degree),n,degree),
                            integral(power(g,k,degree),degree),degree),degree)
    rhs = exponential(mul(power(inverse(g,degree),n,degree),
                          integral(power(g,k,degree),degree),degree),degree)
    require(g == rhs, f'literal defining equation, n={n}, k={k}')
    return g


def threshold(n):
    require(n >= 1, 'one-based threshold domain')
    return 8*((n-1)//3)+3*((n-1)%3)


def main():
    a = derivatives(3,8,DEGREE)
    require(a[:len(PUBLISHED)] == PUBLISHED, 'published DATA terms 0..19')
    fixture = Path(__file__).with_name('b232692.txt')
    rows = [tuple(map(int,line.split())) for line in fixture.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    require([r for r,_ in rows] == list(range(301)), 'source table covers exactly 0..300')
    require(all(a[r] == value for r,value in rows), 'all published b-file terms')
    print('FIDELITY: DATA 0..19 and all 301 b-file coefficients 0..300 agree')
    print('BFILE_SHA256:',sha256(fixture.read_bytes()).hexdigest())
    independent = ode_coefficients(DEGREE)
    require(a == independent, 'independent polynomial-ODE recurrence at every degree 0..430')
    print('INDEPENDENT_ODE: all 431 coefficients 0..430 agree')
    require(all(v > 0 for v in a[:430]), 'all coefficients before degree 430 strictly positive')
    require(a[430] < 0, 'degree-430 coefficient strictly negative')
    print('FIRST_NEGATIVE_DEGREE:',430)
    print('A430_FACTORIAL_NORMALIZED:',a[430])
    print('A430_SHA256:',sha256(str(a[430]).encode('ascii')).hexdigest())
    print('ORDINARY_COEFFICIENT: A430_FACTORIAL_NORMALIZED / 430! < 0')
    expected = [r for r in range(80) if r%8 in (0,3,6)]
    require([threshold(i) for i in range(1,len(expected)+1)] == expected,
            'threshold matches increasing residue enumeration')
    # Literal formal integral/exp fixed point is independent of the P recurrence.
    small = 14
    pairs = [(1,3),(2,6),(3,7),(3,8),(6,16),(9,24),(12,32)]
    for n,k in pairs:
        g = fixed_point(n,k,small)
        require([factorial(r)*g[r] for r in range(small+1)] == derivatives(n,k,small),
                f'literal-source fidelity vs derivative recurrence for n={n},k={k}')
    print('LITERAL_SOURCE: exact Fraction fixed-point checks through degree 14 for',pairs)
    print('SOURCE_EXAMPLES: A232690 and A232691 first six EGF coefficients:',
          derivatives(1,3,5),derivatives(2,6,5))
    require(derivatives(1,3,5)==[1,1,2,7,33,202], 'A232690 example')
    require(derivatives(2,6,5)==[1,1,3,19,161,1857], 'A232691 example')
    for m in range(1,7):
        gm = derivatives(3*m,8*m,DEGREE)
        product = [1]+[0]*DEGREE
        for _ in range(m):
            product = egf_mul(product,gm,DEGREE)
        require(product == [m**r*a[r] for r in range(DEGREE+1)],
                f'power/rescaling identity through degree 430, m={m}')
        negatives = [r for r,v in enumerate(gm) if v < 0]
        require(negatives and min(negatives)<=430, f'negative coefficient for m={m}')
        require(threshold(3*m)==8*m-2<=8*m, 'published one-based threshold')
        require(3*(3*m)-(3*m)//3==8*m, 'shifted zero-based threshold')
        print(f'FAMILY: m={m}; thresholds={8*m-2},{8*m}; k={8*m}; '
              f'first_negative={min(negatives)}; power_identity_0..430=OK')
    print('FAMILY_SCOPE: m=1..6 checked; arbitrary m>=1 requires the symbolic uniqueness proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

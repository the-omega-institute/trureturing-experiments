#!/usr/bin/env python3
"""Exact finite checks for OEIS A138738; no unbounded proof is claimed."""
from fractions import Fraction as Q
from math import factorial, comb

# OEIS A138738 revision 6, DATA (indices 0 through 21).
SOURCE = [
    1, 1, 1, 7, 25, 121, 1561, 9871, 101137, 1293265,
    15765841, 226501111, 3355388521, 55348183177, 981791212585,
    19611613003711, 387083895223201, 8241127257050401,
    193814428658623777, 4572332597344204135,
    116627907470591924281, 3142752742361513926681,
]
LIMIT = 40


def multiply(a, b):
    n = len(a) - 1
    return [sum((a[j] * b[k-j] for j in range(k+1)), Q(0))
            for k in range(n+1)]


def exp_differential(g):
    """Compute exp(g) from f'=g'f, requiring g(0)=0."""
    assert g[0] == 0
    f = [Q(1)] + [Q(0)] * (len(g)-1)
    for k in range(1, len(g)):
        f[k] = sum((j*g[j]*f[k-j] for j in range(1, k+1)), Q(0))/k
    return f


def exp_powers(g):
    """Independent direct sum of powers g^j/j! modulo x^(n+1)."""
    assert g[0] == 0
    n = len(g)-1
    term = [Q(1)] + [Q(0)]*n
    total = term[:]
    for j in range(1, n+1):
        term = [c/j for c in multiply(term, g)]
        total = [a+b for a, b in zip(total, term)]
    return total


def tower(depth, n, first=1, exp=exp_differential):
    f = [Q(1)] + [Q(0)]*n
    for d in range(depth, first-1, -1):
        g = [Q(0)]*(n+1)
        for k in range(d, n+1):
            g[k] = f[k-d]
        f = exp(g)
    return f


def multiplicities(n):
    """Enumerate all weighted partitions of n-1 with parts 2,...,n."""
    def visit(d, remainder, previous, factor, total_count, denominator):
        if d > n:
            if remainder == 0:
                yield factor, total_count, denominator
            return
        for m in range(remainder//d+1):
            yield from visit(d+1, remainder-d*m, m,
                             factor*Q(previous**m, factorial(m)),
                             total_count+m, denominator*factorial(m))
    return visit(2, n-1, 1, Q(1), 0, 1)


def main():
    f = tower(LIMIT, LIMIT)
    a = [factorial(n)*f[n] for n in range(LIMIT+1)]
    assert a[:len(SOURCE)] == SOURCE
    print('FIDELITY: all 22 DATA terms (n=0..21) match OEIS revision 6')
    assert f[:6] == [Q(1), Q(1), Q(1, 2), Q(7, 6), Q(25, 24), Q(121, 120)]
    print('FIDELITY: displayed EGF expansion through degree 5 matches')
    assert all(x.denominator == 1 and (x.numerator-1) % 6 == 0 for x in a)
    print('SETTLEMENT CHECK: integral a(n) == 1 (mod 6) for n=0..40')
    for n in range(LIMIT+1):
        assert tower(n, n)[n] == f[n]
    for depth in [LIMIT+1, LIMIT+5]:
        assert tower(depth, LIMIT) == f
    print('STABILIZATION CHECK: depth n agrees at coefficient n for n=0..40; depths 40,41,45 agree through degree 40')
    assert tower(24, 24, exp=exp_powers) == f[:25]
    print('INDEPENDENT EXP CHECK: direct power sum agrees through degree 24')
    f2 = tower(LIMIT, LIMIT, first=2)
    c = [Q(0)] + [factorial(n)*f2[n-1] for n in range(1, LIMIT+1)]
    count = 0
    for n in range(1, 21):
        total = Q(0)
        for factor, m, denominator in multiplicities(n):
            term = factorial(n)*factor
            total += term
            if n >= 2:
                assert m >= 1 and n >= 2*m+1
                multinomial = Q(factorial(m), denominator)
                gap = factorial(n)//factorial(m)
                assert multinomial.denominator == 1 and gap % 6 == 0
                assert term.denominator == 1 and term.numerator % 6 == 0
            count += 1
        assert total == c[n]
    print(f'MULTIPLICITY CHECK: c(n) agrees for n=1..20; {count} weighted-partition terms checked')
    assert c[0] == 0 and c[1] == 1 and c[2] == 0
    assert all(x.denominator == 1 and x.numerator % 6 == 0 for x in c[2:])
    for n in range(LIMIT):
        assert a[n+1] == sum(comb(n,j)*c[j+1]*a[n-j] for j in range(n+1))
    print('RECURRENCE CHECK: a(n+1)=sum binom(n,j)c(j+1)a(n-j) for n=0..39')
    assert a[3]-1 == 6 and a[6] == 1561
    print('SHARPNESS CHECK: a(3)-1=6; distinction from A138292: a(6)=1561')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

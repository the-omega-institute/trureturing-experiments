#!/usr/bin/env python3
"""Exact fidelity and gap-certificate checks for arXiv:2606.12484v1, C1."""
from fractions import Fraction as Q
from math import gcd, prod
import sympy as sp


def factor(n):
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def source_functions(n):
    factors = factor(n)
    phi = prod(p ** (e - 1) * (p - 1) for p, e in factors)
    psi = prod(p ** (e - 1) * (p + 1) for p, e in factors)
    sigma = prod((p ** (e + 1) - 1) // (p - 1) for p, e in factors)
    return phi, psi, sigma


def lhs(a, b, c):
    return (a + b + c) * (Q(1, a) + Q(1, b) + Q(1, c))


def rhs(n, k):
    return Q(5 * (n * n - 1) ** k + 2 * (n + 1) ** (2 * k)
             + 2 * (n - 1) ** (2 * k), (n * n - 1) ** k)


def gap_terms(x, y, r):
    return ((x - r) * (r * x - 1) / (r * x),
            (y - r) * (r * y - 1) / (r * y),
            (x - y) ** 2 / (x * y))


def main():
    # An exact symbolic identity, checked independently of arithmetic-function data.
    x, y, r = sp.symbols('x y r', nonzero=True)
    identity_lhs = (1 + x + y) * (1 + 1 / x + 1 / y) - 5 - 2 * (r + 1 / r)
    assert sp.cancel(identity_lhs - sum(gap_terms(x, y, r))) == 0
    print('Symbolic gap identity: exact residual 0')

    assert source_functions(1) == (1, 1, 1)
    fidelity = prime_count = composite_count = case_count = 0
    examples = []
    for n in range(2, 301):
        phi, psi, sigma = source_functions(n)
        # Independent counting definitions and the existing D5 psi representation.
        assert phi == sum(gcd(j, n) == 1 for j in range(1, n + 1))
        assert sigma == sum(d for d in range(1, n + 1) if n % d == 0)
        radical = prod(p for p, _ in factor(n))
        assert radical > 0 and n % radical == 0
        assert psi == (n // radical) * prod(p + 1 for p, _ in factor(n))
        assert 0 < phi <= n - 1 and psi >= n + 1 and sigma >= n + 1
        # Source fidelity: Theorem 1, equation (1), exactly as printed.
        assert lhs(phi, psi, sigma) >= Q(9 * n * n - 1, n * n - 1)
        assert rhs(n, 1) == Q(9 * n * n - 1, n * n - 1)
        fidelity += 1
        is_prime = factor(n) == [(n, 1)]
        prime_count += is_prime
        composite_count += not is_prime
        for k in range(1, 13):
            a, b, c = phi ** k, psi ** k, sigma ** k
            a0, b0 = (n - 1) ** k, (n + 1) ** k
            xx, yy, rr = Q(b, a), Q(c, a), Q(b0, a0)
            assert xx >= rr > 1 and yy >= rr
            terms = gap_terms(xx, yy, rr)
            gap = lhs(a, b, c) - rhs(n, k)
            assert rhs(n, k) == 5 + 2 * (rr + 1 / rr)
            assert gap == sum(terms) and all(t >= 0 for t in terms)
            assert gap >= 0
            if is_prime:
                assert gap == 0
            else:
                assert gap > 0
            case_count += 1
            if n in (2, 4, 6) and k == 4:
                examples.append((n, k, str(gap)))
    rational_cases = 0
    for rr in (Q(1), Q(3, 2), Q(2), Q(5)):
        for i in range(13):
            for j in range(13):
                xx, yy = rr + Q(i, 3), rr + Q(j, 5)
                terms = gap_terms(xx, yy, rr)
                assert ((1 + xx + yy) * (1 + 1 / xx + 1 / yy)
                        - 5 - 2 * (rr + 1 / rr)) == sum(terms)
                assert all(t >= 0 for t in terms)
                rational_cases += 1
    print('Definitions: n=1 normalizations; independent counts and exact psi division for 2<=n<=300')
    print(f'Source fidelity: Theorem 1 for {fidelity} integers; k=1 RHS equivalence')
    print(f'Conjecture and gap identity: {case_count} pairs, 2<=n<=300, 1<=k<=12')
    print(f'Prime equality: {prime_count} primes x 12 exponents; strict gaps for {composite_count} composites x 12 exponents')
    print(f'Independent rational boundary grid: {rational_cases} triples, including R=1 and x=y=R')
    for n, k, gap in examples:
        print(f'Example n={n}, k={k}: exact gap={gap}')
    print('Finite checks and symbolic algebra do not certify the universal Lean claim.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

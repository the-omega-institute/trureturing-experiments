#!/usr/bin/env python3
"""Exact finite checks for Aliabadi et al., arXiv:2607.24847v2, Question 2.17.

Finite checks support fidelity and the proposed limit argument; they do not
establish an infinite limit. All evidence uses integers or Fraction.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, prod


def primes_to(bound):
    ps = []
    for n in range(2, bound + 1):
        if all(n % p for p in ps if p * p <= n):
            ps.append(n)
    return ps


def totient_sieve(bound):
    phi = list(range(bound + 1))
    for p in range(2, bound + 1):
        if phi[p] == p:
            for n in range(p, bound + 1, p):
                phi[n] -= phi[n] // p
    return phi


def minfac_sieve(bound):
    spf = [0] * (bound + 1)
    for p in range(2, bound + 1):
        if spf[p] == 0:
            for n in range(p, bound + 1, p):
                if spf[n] == 0:
                    spf[n] = p
    return spf


def direct_cyclic_chowla(n):
    # Compute every element's additive order directly, without the criterion
    # or the divisor/totient formula. Sorted orders give the largest feasible k.
    orders = sorted((n // gcd(a, n) for a in range(n)), reverse=True)
    return max(k for k in range(1, n) if orders[k - 1] > k)


def exhaustive_subsets(n):
    # Independent literal-definition oracle for the smallest groups.
    return max(len(s) for k in range(1, n)
               for s in combinations(range(n), k)
               if all(n // gcd(a, n) > len(s) for a in s))


def phi_probability(p, cutoff, threshold):
    # Finite Bernoulli prime-pattern law, compared using rational products.
    states = [(Q(1), Q(1))]
    for q in primes_to(cutoff):
        if q <= p:
            continue
        states = [(x, w * Q(q - 1, q)) for x, w in states] + [
            (x * Q(q, q - 1), w * Q(1, q)) for x, w in states]
    assert sum((w for x, w in states), Q(0)) == 1
    assert len({x for x, w in states}) == len(states)
    return sum((w for x, w in states if x <= threshold), Q(0))


def mixture(cutoff):
    residual = Q(1)
    total = Q(0)
    for p in primes_to(cutoff):
        dp = residual / p
        total += dp * phi_probability(p, cutoff, Q(p - 1))
        residual *= Q(p - 1, p)
    return total, residual


def rational_certificate(m):
    residual = Q(1)
    lower, upper = Q(0), Q(0)
    for p in primes_to(m):
        dp = residual / p
        lower += dp * phi_probability(p, m ** 3, Q((p - 1) * m, m + 1))
        upper += dp * phi_probability(p, m ** 3, Q(p - 1))
        residual *= Q(p - 1, p)
    return lower - Q(m + 1, m ** 3), upper + residual


def main():
    bound = 4000
    phi = totient_sieve(bound)
    spf = minfac_sieve(200000)
    for n in range(2, 13):
        assert direct_cyclic_chowla(n) == exhaustive_subsets(n)
    print('Literal subset oracle agrees with sorted element orders for 2 <= n <= 12')
    for n in range(2, bound + 1):
        c = direct_cyclic_chowla(n)
        divisors = [d for d in range(1, n + 1) if n % d == 0]
        div_formula = max(m for m in range(1, n)
                          if sum(phi[d] for d in divisors if d > m) >= m)
        # Source Theorem 2.4 and Proposition 2.6, using independent computations.
        assert c == div_formula
        assert (c == phi[n]) == (phi[n] + 1 >= n // spf[n])
        assert (c == phi[n]) == (Q(phi[n] + 1, n) >= Q(1, spf[n]))
    print(f'Source cyclic divisor formula and iff criterion: all {bound - 1} moduli 2..{bound}')
    for cutoff in (7, 11, 13):
        ps = primes_to(cutoff)
        modulus = prod(ps)
        count = 0
        for r in range(modulus):
            dividing = [p for p in ps if r % p == 0]
            if dividing and prod((Q(p, p - 1) for p in dividing[1:]), start=Q(1)) <= dividing[0] - 1:
                count += 1
        finite, residual = mixture(cutoff)
        assert Q(count, modulus) == finite
        assert residual == Q(sum(gcd(r, modulus) == 1 for r in range(modulus)), modulus)
        print(f'CRT finite-prime mixture Y={cutoff}: {count}/{modulus} = {finite}; residual={residual}')
    # Exact finite counts of the actual source criterion (not a float estimate).
    limit = 200000
    phis = totient_sieve(limit)
    for nmax in (4000, 20000, 200000):
        count = sum(phis[n] + 1 >= n // spf[n] for n in range(2, nmax + 1))
        print(f'Empirical source density N={nmax}: {count}/{nmax} = {Q(count, nmax)}')
    for m in (2, 3, 4):
        lo, hi = rational_certificate(m)
        assert lo <= hi
        # Finite observations are compared, without asserting an error bound for N.
        count = sum(phis[n] + 1 >= n // spf[n] for n in range(2, limit + 1))
        empirical = Q(count, limit)
        assert lo <= empirical <= hi
        print(f'Proposed limit enclosure m={m}: [{lo}, {hi}]; empirical N={limit} lies inside')
    witness_primes = (3, 5, 7, 11, 13, 17, 19, 23)
    P, phiP = prod(witness_primes), prod(p - 1 for p in witness_primes)
    assert (P, phiP) == (111546435, 36495360)
    assert 3 * (phiP + 1) < P
    print(f'Follow-up upper-bound witness: P={P}, phi(P)={phiP}, 3*(phi(P)+1)={3 * (phiP + 1)} < P')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

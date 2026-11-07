#!/usr/bin/env python3
"""Exact recurrence checks for Sun's arXiv:2608.13192v1, Conjecture 4.1.

Finite checks support the proposed all-prime mechanism; they do not prove it.
Only Python's standard library, integers, and Fraction arithmetic are used.
"""

from fractions import Fraction as Q
from math import factorial, gcd, prod
import json


def primes_through(bound):
    return [n for n in range(2, bound + 1)
            if all(n % d for d in range(2, n) if d * d <= n)]


def sequence(x, length):
    """The literal initial values and three-term recurrence, over Q."""
    x = Q(x)
    a = [Q(1), x]
    for n in range(1, length - 1):
        a.append(((2*n + 1) * (n*(n + 1) + x) * a[n]
                  - n**3 * a[n-1]) / (n + 1)**3)
    return a[:length]


def residue(x, modulus):
    x = Q(x)
    assert gcd(x.denominator, modulus) == 1
    return x.numerator * pow(x.denominator, -1, modulus) % modulus


def moment(a, exponent):
    return sum((Q((2*n + 1)**exponent) * v*v
                for n, v in enumerate(a)), Q(0))


def product_formula(n, x):
    """Source Theorem 3.1, v_n(x), with empty product equal to one."""
    return sum((prod((Q((2*n+1)**2 - (2*r-1)**2)
                      * (2*x - 1 - (2*r-1)**2)
                      for r in range(1, k+1)), start=Q(1))
                / (16**k * factorial(k)**4)
                for k in range(n+1)), Q(0))


def add(a, b, modulus):
    return [((a[i] if i < len(a) else 0)
             + (b[i] if i < len(b) else 0)) % modulus
            for i in range(max(len(a), len(b)))]


def mul(a, b, modulus):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] = (c[i+j] + x*y) % modulus
    return c


def scale(a, c, modulus):
    return [x*c % modulus for x in a]


def eval_poly(a, x, modulus):
    value = 0
    for c in reversed(a):
        value = (value*x + c) % modulus
    return value


def reduced_F(p):
    """Independent polynomial arithmetic for the route's truncation lemma."""
    total, basis, b = [0], [1], 1
    for k in range(p):
        if k:
            root = (2*k - 1)**2
            basis = mul(basis, [-root, 1], p)
            b = b * (4-root) * pow(16*k**4, -1, p) % p
        total = add(total, scale(basis, b, p), p)
    while total and total[-1] == 0:
        total.pop()
    return total


def main():
    primes = primes_through(97)
    fidelity_product = 0
    for x in [Q(5, 2), Q(8), Q(-7, 3), Q(13, 8)]:
        a = sequence(x, 16)
        for n in range(16):
            assert a[n] == product_formula(n, x)
            fidelity_product += 1
    print("FIDELITY Theorem 3.1: exact product formula,", fidelity_product, "checks")

    fidelity_corollary = 0
    for p in primes:
        if p > 3:
            linear = moment(sequence(Q(5, 2), p), 1)
            assert residue(linear, p**4) == 0
            fidelity_corollary += 1
    print("FIDELITY Corollary 4.3, first congruence: linear squared moment",
          "at x=5/2 is 0 mod p^4;", fidelity_corollary, "primes through 97")

    fidelity_theorem = 0
    for p in primes:
        if not 3 < p <= 43:
            continue
        for r in range(p):
            if pow((2*r-1) % p, (p-1)//2, p) != 1:
                continue
            x0 = next(j for j in range(1, (p-1)//2 + 1)
                      if ((2*j-1)**2 - (2*r-1)) % p == 0)
            for shift in [Q(0), Q(1, 3)]:
                x = Q(r) + p*shift
                x1 = (2*x - 1 - (2*x0-1)**2) / p
                rhs = ((x1+4*x0-2) / (2*(2*x-1)) * p**2
                       + x1*(x1+8*x0-4) / (8*(2*x-1)**2) * p**3)
                assert residue(moment(sequence(x, p), 1) - rhs, p**4) == 0
                fidelity_theorem += 1
    print("FIDELITY Theorem 4.4: exact linear congruence,",
          fidelity_theorem, "rational inputs at primes 5 through 43")

    shifts = [Q(0), Q(-2), Q(-1), Q(1), Q(2), Q(1, 2), Q(-1, 3), Q(2, 3)]
    family_count = 0
    summary = []
    for p in primes:
        if p <= 7:
            continue
        expected = residue(Q(p, 64), p*p)
        assert p > 4*2 and expected != 0
        f = reduced_F(p)
        h = (p-3)//2
        assert len(f) - 1 == h
        assert f[-1] == residue(Q(-1, 8), p)
        top = f[-1]**2 % p
        assert 5+4*h == 2*p-1 and top == residue(Q(1, 64), p)
        # Odd-polynomial reflection, including the boundary and a term
        # whose coefficient vanishes modulo p but whose degree is larger.
        for powers in [{1: 1}, {p: 3}, {2*p-1: 7},
                       {1: 2, p: 3, 2*p-1: 5, 2*p+1: p}]:
            value = sum(sum(c*pow(2*n+1, d, p*p) for d, c in powers.items())
                        for n in range(p)) % (p*p)
            assert value == p * powers.get(2*p-1, 0) % (p*p)
        for t in shifts:
            x = Q(5, 2) + p*t
            assert residue(x-Q(1, 2), p) != 0
            assert residue(x+Q(1, 2), p) != 0
            assert residue(2*x-1, p) == 4
            assert pow(residue(2*x-1, p), (p-1)//2, p) == 1
            a = sequence(x, p)
            assert all(gcd(v.denominator, p) == 1 for v in a)
            actual = residue(moment(a, 5), p*p)
            assert actual == expected
            family_count += 1
        summary.append([p, expected, h, f[-1]])
    print("FAMILY m=2, x=5/2+p*t:", family_count,
          "exact rational checks; p ranges over all primes 7<p<=97")
    print("FAMILY rows [p, quintic residue mod p^2, degree(F mod p), leading coefficient]:")
    print(json.dumps(summary))

    for p in [11, 13]:
        for x in [Q(5, 2), Q(p+5, 2)]:
            a = sequence(x, p)
            s = moment(a, 5)
            record = {"p": p, "x": str(x), "m": 2,
                      "sequence": list(map(str, a)),
                      "residues_mod_p2": [residue(v, p*p) for v in a],
                      "weighted_terms_mod_p2": [residue((2*n+1)**5*v*v, p*p)
                                                 for n, v in enumerate(a)],
                      "exact_sum": str(s), "sum_mod_p2": residue(s, p*p)}
            if p == 11 and x == Q(5, 2):
                assert s == Q(653267525096522474214586592269, 295147905179352825856)
                assert record["residues_mod_p2"] == [1,63,23,87,12,2,100,21,23,96,56]
                assert record["weighted_terms_mod_p2"] == [1,97,23,43,23,0,109,56,120,68,120]
            print("WITNESS", json.dumps(record, sort_keys=True))

    # Independent recurrence in (Z/121Z)[X], followed by coefficient comparison.
    modulus = 121
    a = [[1], [0, 1]]
    for n in range(1, 10):
        numerator = add(scale(mul([n*(n+1), 1], a[n], modulus), 2*n+1, modulus),
                        scale(a[n-1], -n**3, modulus), modulus)
        a.append(scale(numerator, pow((n+1)**3, -1, modulus), modulus))
    s_poly = [0]
    for n, v in enumerate(a):
        s_poly = add(s_poly, scale(mul(v, v, modulus), (2*n+1)**5, modulus), modulus)
    factorization = [-44 % modulus]
    for factor in [[-5, 1], [-3, 1], [-2, 1], [-1, 1], [4, 1], [2, -2, 1]]:
        factorization = mul(factorization, factor, modulus)
    assert add(s_poly, scale(factorization, -1, modulus), modulus) == [0]*len(s_poly)
    for x in range(modulus):
        assert eval_poly(s_poly, x, modulus) == residue(moment(sequence(Q(x), 11), 5), modulus)
    assert all(eval_poly(s_poly, 8+11*t, modulus) == 55 for t in range(11))
    print("POLYNOMIAL p=11: S(X)=-44(X-5)(X-3)(X-2)(X-1)(X+4)(X^2-2X+2) mod 121")
    print("POLYNOMIAL: coefficient identity and all 121 residues checked; x=8+11*t gives 55")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()

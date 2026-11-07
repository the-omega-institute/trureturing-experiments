#!/usr/bin/env python3
"""Exact checks for the positive symmetric unit-step weighting existence question.

Python 3 and SymPy. All evidence uses integers, Fraction, symbolic algebra,
or Q[x]/Phi_N(x); no numerical approximation is used.
"""
from fractions import Fraction
from functools import reduce
from math import gcd, lcm
import sys
import sympy as sp


class Cyclotomic:
    """Integer coordinates in Z[zeta_N], with the standard polynomial basis."""
    def __init__(self, n):
        self.n = n
        x = sp.Symbol('x')
        phi = sp.Poly(sp.cyclotomic_poly(n, x), x)
        self.d = phi.degree()
        self.rel = tuple(int(phi.nth(i)) for i in range(self.d))
        self.zero = (0,) * self.d
        self.one = (1,) + (0,) * (self.d - 1)
        self.powers = [self.one]
        for _ in range(1, n):
            v = self.powers[-1]
            shifted = [0] + list(v[:-1])
            self.powers.append(tuple(a - v[-1] * b for a, b in zip(shifted, self.rel)))
        assert self.shift(self.powers[-1], 1) == self.one

    def scale(self, k, v):
        return tuple(k * a for a in v)

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def shift(self, v, e):
        out = [0] * self.d
        for i, c in enumerate(v):
            if c:
                for j, a in enumerate(self.powers[(i + e) % self.n]):
                    out[j] += c * a
        return tuple(out)

    def conj(self, v):
        out = self.zero
        for i, c in enumerate(v):
            if c:
                out = self.add(out, self.scale(c, self.powers[-i % self.n]))
        return out

    def mul(self, a, b):
        out = self.zero
        for i, c in enumerate(a):
            if c:
                out = self.add(out, self.scale(c, self.shift(b, i)))
        return out


def ramanujan(n, k):
    d = n // gcd(n, k)
    numerator = int(sp.mobius(d)) * int(sp.totient(n))
    denominator = int(sp.totient(d))
    assert numerator % denominator == 0
    return numerator // denominator


def inverse_numerators(spectrum):
    n = len(spectrum)
    ring = Cyclotomic(n)
    rows = []
    for j in range(n):
        value = ring.zero
        for k, c in enumerate(spectrum):
            if c:
                value = ring.add(value, ring.scale(c, ring.powers[-k * j % n]))
        rows.append(value)
    return ring, rows


def check_weights(n, spectrum, base, perturbation):
    """Inverse weights are rows/n; the norm estimate certifies strict positivity."""
    ring, rows = inverse_numerators(spectrum)
    assert len(spectrum) == n
    assert spectrum == [spectrum[-k % n] for k in range(n)]
    assert all(spectrum[k] == base * ramanujan(n, k) + perturbation[k]
               for k in range(n))
    perturbation_bound = Fraction(sum(abs(c) for c in perturbation), n)
    assert base > perturbation_bound
    for j, row in enumerate(rows):
        assert row == rows[-j % n] == ring.conj(row)
        if gcd(j, n) != 1:
            assert row == ring.zero
        # On a unit, w_j = base + delta_j and |delta_j| <= perturbation_bound.
        # The strict rational inequality above therefore certifies w_j > 0.
    # Recompute the entire forward transform in the quotient ring.
    for k in range(n):
        value = ring.zero
        for j, row in enumerate(rows):
            value = ring.add(value, ring.shift(row, k * j))
        assert value == ring.scale(n * spectrum[k], ring.one)
    # Exact prime-coset vanishing, for every index and prime divisor.
    for q in sp.factorint(n):
        q = int(q)
        for k in range(n):
            assert sum(spectrum[(k + ell * (n // q)) % n] for ell in range(q)) == 0
    return perturbation_bound


def check_revival(n, spectrum, time_multiple):
    """At t=2*pi*time_multiple compute the full column of exp(-it A)."""
    den, num = time_multiple.denominator, time_multiple.numerator
    order = lcm(n, den)
    ring = Cyclotomic(order)
    phases = [ring.powers[-num * value * (order // den) % order]
              for value in spectrum]
    column = []  # Numerators, with common denominator n.
    for j in range(n):
        value = ring.zero
        for k, phase in enumerate(phases):
            value = ring.add(value, ring.shift(phase, -k * j * (order // n)))
        column.append(value)
    assert all(value == ring.zero for j, value in enumerate(column)
               if j not in (0, n // 2))
    alpha, beta = column[0], column[n // 2]
    assert beta != ring.zero
    assert ring.add(ring.mul(alpha, ring.conj(alpha)),
                    ring.mul(beta, ring.conj(beta))) == ring.scale(n * n, ring.one)
    return alpha == ring.zero


def source_fidelity():
    # The source's unweighted n=8 instance has no nontrivial QFR.
    s8 = [ramanujan(8, k) for k in range(8)]
    assert s8 == [4, 0, 0, 0, -4, 0, 0, 0]
    # Each parity contains a zero eigenvalue, hence both revival phases would be 1.
    assert s8[1] == s8[2] == 0
    # Also reproduce the source's n=6 probability |beta|^2=3/4 at 2*pi/3.
    s6 = [ramanujan(6, k) for k in range(6)]
    assert s6 == [2, 1, -1, -2, -1, 1]
    assert not check_revival(6, s6, Fraction(1, 3))
    ring = Cyclotomic(6)
    beta_num = ring.zero
    for k, eigenvalue in enumerate(s6):
        beta_num = ring.add(beta_num,
                            ring.powers[(-2 * eigenvalue - 3 * k) % 6])
    beta_sq_num = ring.mul(beta_num, ring.conj(beta_num))
    assert beta_sq_num == ring.scale(27, ring.one)  # 27/6^2 = 3/4.
    print('SOURCE_FIDELITY: unweighted n=8 spectrum and no-QFR obstruction; n=6 |beta|^2=3/4')


def check_supplied_eight():
    r = sp.sqrt(2)
    weights = [0, 1 + r / 4, 0, 1 - r / 4, 0, 1 - r / 4, 0, 1 + r / 4]
    spectrum = [sp.simplify(sum(weights[j] *
                 (sp.cos(sp.pi * k * j / 4) + sp.I * sp.sin(sp.pi * k * j / 4))
                 for j in range(8))) for k in range(8)]
    assert spectrum == [4, 1, 0, -1, -4, -1, 0, 1]
    assert r.is_positive and (4 - r).is_positive
    assert check_revival(8, list(map(int, spectrum)), Fraction(1, 2))
    print('SUPPLIED_n8: exact weights, spectrum (4,1,0,-1,-4,-1,0,1), exp(-i*pi*A)=P_4')


def check_excluded_twelve():
    # Every symmetric positive unit weighting is (a,b) on +/-1 and +/-5.
    a, b = sp.symbols('a b', positive=True)
    spectrum = [sp.expand(2 * a * sp.cos(sp.pi * k / 6) +
                         2 * b * sp.cos(5 * sp.pi * k / 6)) for k in range(12)]
    assert spectrum[0] == 2 * spectrum[2]
    assert spectrum[2] == a + b
    assert spectrum[3] == 0
    # Equal even-index phases z imply z=z^2; z is invertible.
    # The zero odd eigenvalue forces the odd phase to be 1.
    z, zi = sp.symbols('z zi')
    basis = sp.groebner([z**2 - z, z * zi - 1], z, zi)
    assert list(basis) == [z - 1, zi - 1]
    print('EXCLUDED_n12: all positive real a,b and all times; exact phase ideal gives z_even=z_odd=1')


def power_two(h):
    n = 2**h
    if h in (1, 2):
        spectrum = [ramanujan(n, k) for k in range(n)]
        check_weights(n, spectrum, 1, [0] * n)
        assert check_revival(n, spectrum, Fraction(1, 4))
        return
    f = [0 if k % 2 == 0 else (1 if k < n // 4 or k > 3 * n // 4 else -1)
         for k in range(n)]
    spectrum = [ramanujan(n, k) + f[k] for k in range(n)]
    assert check_weights(n, spectrum, 1, f) == Fraction(1, 2)
    assert check_revival(n, spectrum, Fraction(1, 2))


def twice_odd_prime_power(p, a):
    m, q = p**a, p**(a - 1)
    n = 2 * m
    reps = {r % m for r in range(-(q // 2), q // 2 + 1)}
    assert len(reps) == q and {-r % m for r in reps} == reps
    assert all(sum((k + ell * q) % m in reps for ell in range(p)) == 1
               for k in range(m))
    f = [p - 1 if k in reps else -1 for k in range(m)]
    mu = [p * ramanujan(m, k) + f[k] for k in range(m)]
    inv_two = pow(2, -1, m)
    perturbation = [(-1)**k * f[k * inv_two % m] for k in range(n)]
    spectrum = [(-1)**k * mu[k * inv_two % m] for k in range(n)]
    assert all(spectrum[k] % p == (-1 if k % 2 == 0 else 1) % p for k in range(n))
    assert check_weights(n, spectrum, p, perturbation) == Fraction(2 * (p - 1), p)
    assert not check_revival(n, spectrum, Fraction(1, p))


def all_small_cases(limit):
    allowed, excluded = [], []
    for n in range(2, limit + 1):
        factors = sp.factorint(n)
        if set(factors) == {2}:
            power_two(int(factors[2]))
            allowed.append(n)
        elif factors.get(2) == 1 and len(factors) == 2:
            p = next(int(p) for p in factors if p != 2)
            twice_odd_prime_power(p, int(factors[p]))
            allowed.append(n)
        else:
            excluded.append(n)
            if n % 2 == 0:
                odd_primes = [int(p) for p in factors if p != 2]
                # Necessary phase orders, derived from bipartition and coset sums.
                constraints = ([2] if n % 4 == 0 else []) + odd_primes
                assert reduce(gcd, constraints) == 1
            # For odd n, inversion fixing 0 leaves no distinct involutive target.
    print('SMALL_CONSTRUCTIONS n<=%d: %s' % (limit, allowed))
    print('SMALL_OBSTRUCTIONS n<=%d: %s' % (limit, excluded))
    print('Every construction: exact inverse/forward DFT, reality, symmetry, zero nonunit support,')
    print('strict-positive norm certificate, every prime-coset identity, full revival column and normalization.')


def main():
    source_fidelity()
    check_supplied_eight()
    check_excluded_twelve()
    all_small_cases(64)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

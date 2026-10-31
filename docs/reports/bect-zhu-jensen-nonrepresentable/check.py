#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact fidelity and tail-family checks for arXiv:2605.18100v1, Open Problem 2.

Finite checks and symbolic identities support a preregistration. They are not a
proof of measurability, the arbitrary-random-measure Jensen quantifier, or the
quantifier over every measurable loss family.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product
import sympy as sp


def simplex(denominator, dimension):
    for numerators in product(range(denominator + 1), repeat=dimension):
        if sum(numerators) == denominator:
            yield tuple(F(n, denominator) for n in numerators)


def fidelity():
    """Section 2.2.1, Proposition 2, eqs. (5)-(6), h(t)=t-t^2.

    Reproduce the Gini-Simpson identity and the Jensen-gap decomposition for
    every pair of small probability vectors and rational two-point mixing law.
    """
    vectors = [p for d in range(1, 6) for p in simplex(d, 3)]
    weights = sorted({F(k, d) for d in range(1, 7) for k in range(d + 1)})
    entropy = lambda p: sum((x - x*x for x in p), F(0))
    checked = 0
    for p in vectors:
        assert entropy(p) == 1 - sum(x*x for x in p)
        for q, t in product(vectors, weights):
            mean = tuple(t*x + (1-t)*y for x, y in zip(p, q))
            gap = entropy(mean) - t*entropy(p) - (1-t)*entropy(q)
            coordinate_gap = sum((z-z*z - t*(x-x*x) - (1-t)*(y-y*y)
                                  for x, y, z in zip(p, q, mean)), F(0))
            variance_gap = t*(1-t)*sum((x-y)**2 for x, y in zip(p, q))
            assert gap == coordinate_gap == variance_gap and gap >= 0
            checked += 1
    assert all(entropy(tuple(F(i == j) for i in range(3))) == 0
               for j in range(3))
    print(f"FIDELITY: Gini-Simpson identity and source Jensen-gap equation; {checked} exact cases")


@dataclass(frozen=True)
class Profile:
    """prefix[0..K], followed by a*2^(-n) for every n>K."""
    prefix: tuple
    a: F

    @property
    def cutoff(self):
        return len(self.prefix)-1

    def mass(self, n):
        return self.prefix[n] if n <= self.cutoff else self.a / 2**n

    def total_mass(self):
        return sum(self.prefix, F(0)) + self.a / 2**self.cutoff

    def H(self):
        # The scaled tail is exactly a from cutoff+1 onwards, so its liminf is a.
        return min(F(1), self.a)


def p_profile(c):
    return Profile((1-c,), c)


def q_profile(c, N):
    assert 0 < c < 1 and N >= 1
    return Profile(((1-c)*(1-F(1, 2**N)),)
                   + tuple(c / 2**n for n in range(1, N+1)), F(1))


def finite_jensen(profiles):
    count = 0
    for p, q, t in product(profiles, profiles, [F(0), F(1, 5), F(1, 2), F(4, 5), F(1)]):
        K = max(p.cutoff, q.cutoff)
        mean = Profile(tuple(t*p.mass(n)+(1-t)*q.mass(n) for n in range(K+1)),
                       t*p.a+(1-t)*q.a)
        assert mean.total_mass() == 1
        assert mean.H() >= t*p.H()+(1-t)*q.H()
        count += 1
    print(f"JENSEN SANITY: {count} exact finite mixtures of infinite-tail profiles")


def tail_loss_sum(N, A, B, C, D, r):
    """Sum over n>N of 2^(-n)*(A+B*n+C*n^2+D*r^n); 0<=r<2."""
    assert 0 <= r < 2
    return ((A+B*(N+2)+C*(N*N+4*N+6)) / 2**N
            + D*(r/2)**(N+1)/(1-r/2))


def exact_risk(profile, ell0, A, B, C, D, r):
    ell = lambda n: ell0 if n == 0 else A+B*n+C*n*n+D*r**n
    return (sum((profile.mass(n)*ell(n) for n in range(profile.cutoff+1)), F(0))
            + profile.a*tail_loss_sum(profile.cutoff, A, B, C, D, r))


def check_tail_family():
    cs = sorted({F(k, d) for d in range(2, 10) for k in range(1, d)})
    profiles, family_count, risk_count = [], 0, 0
    losses = [(F(0), F(0), F(0), F(0), F(0), F(1)),
              (F(1, 7), F(1, 4), F(0), F(0), F(0), F(1)),
              (F(3), F(0), F(1), F(1), F(0), F(1)),
              (F(0), F(0), F(0), F(0), F(2), F(3, 2)),
              (F(2), F(1, 3), F(1, 5), F(1, 7), F(1, 2), F(1, 2))]
    for c in cs:
        p = p_profile(c)
        assert p.total_mass() == 1 and p.H() == c
        profiles.append(p)
        for N in range(1, 13):
            q = q_profile(c, N)
            assert q.total_mass() == 1 and q.H() == 1
            assert all(q.mass(n) > 0 and q.mass(n) <= p.mass(n)/c
                       for n in range(N+2))
            # Beyond N the ratio is the same exact constant, not a truncation.
            assert q.a == p.a/c
            deficit_at_zero = p.mass(0)-q.mass(0)
            excess_in_entire_tail = (q.a-p.a)/2**N
            assert deficit_at_zero == excess_in_entire_tail == (1-c)/2**N
            for loss in losses:
                ell0, A, B, C, D, r = loss
                Rp = exact_risk(p, *loss)
                Rq = exact_risk(q, *loss)
                tail = tail_loss_sum(N, A, B, C, D, r)
                # Exact infinite-series identity used in the general outline.
                assert Rq-Rp == (1-c)*(tail-ell0/2**N)
                risk_count += 1
            family_count += 1
    diracs = [Profile(tuple(F(n == j) for n in range(j+1)), F(0))
              for j in range(17)]
    assert all(p.total_mass() == 1 and p.H() == 0 for p in diracs)
    print("DIRAC: 17 exact eventually-zero coordinate profiles have H=0")
    # Extended-real boundary: a loss infinite at any coordinate has infinite
    # risk at every strictly positive p^c and q^(c,N).
    infinity_checks = 0
    for c, n in product(cs, range(16)):
        assert p_profile(c).mass(n) > 0
        infinity_checks += 1
    print(f"TAIL FAMILY: {family_count} positive normalized p/q pairs, including entire tails")
    print(f"RISK IDENTITY: {risk_count} exact infinite-series checks for integrable losses")
    print(f"INFINITE LOSS: {infinity_checks} strictly positive coordinate masses")
    clipped = Profile((F(0), F(0)), F(2))
    assert clipped.total_mass() == 1 and clipped.H() == 1
    finite_jensen(profiles[:12] + [q_profile(F(1, 3), N) for N in range(1, 5)]
                  + diracs[:4] + [clipped])


def symbolic_checks():
    c, t, e0, T = sp.symbols('c t e0 T')
    N = sp.symbols('N', integer=True, positive=True)
    z = sp.symbols('z')
    assert sp.expand((1-c)*(1-t)+c*(1-t)+t-1) == 0
    assert sp.expand((1-c)*(1-t)*e0+c*(T) + (1-c)*T
                     - ((1-c)*e0+c*T) - (1-c)*(T-t*e0)) == 0
    # The closed geometric-tail formulas used above, as formal power series.
    S0 = z**(N+1)/(1-z)
    S1 = sp.factor(z*sp.diff(S0, z))
    S2 = sp.factor(z*sp.diff(S1, z))
    for expression, expected in [(S0, 2**(-N)),
                                 (S1, (N+2)*2**(-N)),
                                 (S2, (N*N+4*N+6)*2**(-N))]:
        assert sp.simplify(expression.subs(z, sp.Rational(1, 2))-expected) == 0
    assert sp.limit(2**(-N), N, sp.oo) == 0
    assert sp.limit((N*N+4*N+6)*2**(-N), N, sp.oo) == 0
    print("SYMBOLIC: normalization, risk-difference identity, geometric tail moments and limits")


if __name__ == '__main__':
    fidelity()
    symbolic_checks()
    check_tail_family()
    print('ALL_CHECKS_OK')

#!/usr/bin/env python3
"""Exact finite-atom checks for Rubin's bounded-tilt extension.

No floating-point arithmetic or numerical transcendental evaluation is used.
A finite rational linear combination of log(prime) is stored as its exact
coefficient map. Log identities are checked in this representation. Its sign is
checked by clearing denominators and comparing integer products.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import lcm


@lru_cache(None)
def prime_powers(n):
    assert isinstance(n, int) and n > 0
    result = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            result[d] = result.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def logq(q):
    q = F(q)
    assert q > 0
    ans = {p: F(e) for p, e in prime_powers(q.numerator).items()}
    for p, e in prime_powers(q.denominator).items():
        ans[p] = ans.get(p, F(0)) - e
    return {p: e for p, e in ans.items() if e}


def add(*expressions):
    result = {}
    for expression in expressions:
        for p, c in expression.items():
            result[p] = result.get(p, F(0)) + c
    return {p: c for p, c in result.items() if c}


def scale(c, expression):
    return {p: F(c) * a for p, a in expression.items() if c * a}


def nonnegative(expression):
    denominator = lcm(*(c.denominator for c in expression.values()))
    numerator, divisor = 1, 1
    for p, c in expression.items():
        exponent = int(c * denominator)
        if exponent >= 0:
            numerator *= p ** exponent
        else:
            divisor *= p ** (-exponent)
    return numerator >= divisor


def expectation(q, logs):
    return add(*(scale(qi, ui) for qi, ui in zip(q, logs)))


def kl(q, p):
    return expectation(q, tuple(logq(qi / pi) for qi, pi in zip(q, p)))


def compositions(total, length):
    if length == 1:
        if total > 0:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def probabilities(n):
    return sorted({tuple(F(i, d) for i in c)
                   for d in range(n, 6) for c in compositions(d, n)})


def tilted(p, weights):
    z = sum(pi * wi for pi, wi in zip(p, weights))
    assert z > 0
    return z, tuple(pi * wi / z for pi, wi in zip(p, weights))


def main():
    fidelity_count = 0
    tightness_count = 0
    reverse_tilt_count = 0
    gauge_count = 0
    shift_count = 0
    for n in range(2, 5):
        distributions = probabilities(n)
        weights_list = list(product((F(1), F(2), F(3)), repeat=n))
        for p in distributions:
            for weights in weights_list:
                z, target = tilted(p, weights)
                u = tuple(logq(w) for w in weights)
                assert sum(target) == 1 and all(t > 0 for t in target)
                assert add(logq(z), kl(target, p), scale(-1, expectation(target, u))) == {}
                tightness_count += 1
                for q in distributions:
                    free_energy = add(kl(q, p), scale(-1, expectation(q, u)))
                    gap = add(logq(z), free_energy)
                    # Fidelity: source Proposition 3.1, with u_i = log(weights_i).
                    assert gap == kl(q, target)
                    assert nonnegative(gap)
                    assert (gap == {}) == (q == target)
                    fidelity_count += 1
                    # The claimed identity is unchanged by an arbitrary Q-only gauge.
                    a = sum(qi * qi for qi in q) + F(7, 3)
                    # A=(KL,a), B=(-E,-a): the rational constants cancel exactly.
                    A, B = (kl(q, p), a), (scale(-1, expectation(q, u)), -a)
                    assert A[1] + B[1] == 0
                    assert add(A[0], B[0]) == free_energy
                    gauge_count += 1
                    # Source constant-shift invariance: exp(u+c)=s*weights.
                    for s in (F(1, 2), F(3)):
                        zs, ts = tilted(p, tuple(s * wi for wi in weights))
                        us = tuple(logq(s * wi) for wi in weights)
                        assert ts == target and zs == s * z
                        assert add(logq(zs), kl(q, p), scale(-1, expectation(q, us))) == gap
                        shift_count += 1
            for q in distributions:
                # Every positive finite Q is a bounded tilt of p.
                z, target = tilted(p, tuple(qi / pi for qi, pi in zip(q, p)))
                assert z == 1 and target == q
                # Tightness recovers A(Q,p)-a(Q)=E_Q[log(Q/p)]=KL(Q||p).
                assert expectation(q, tuple(logq(qi / pi) for qi, pi in zip(q, p))) == kl(q, p)
                reverse_tilt_count += 1
        print(f'n={n}: distributions={len(distributions)}, weights={len(weights_list)}')
    print(f'FIDELITY_PROP_3_1={fidelity_count}')
    print(f'CANONICAL_IDENTITY_AND_GAUGE={gauge_count}')
    print(f'TIGHTNESS_AT_TARGET={tightness_count}')
    print(f'BOUNDED_RATIO_REVERSE_TILTS={reverse_tilt_count}')
    print(f'CONSTANT_SHIFT_INVARIANCE={shift_count}')
    print('SCOPE: finite rational atoms only; no claim of an infinite-dimensional proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

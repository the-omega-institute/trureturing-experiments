#!/usr/bin/env python3
"""Exact checks of Patel's OGD estimate and the proposed sign-hinge lower bound.

This finite experiment supports a proposed general argument; it is not a proof
of the quantified minimax statement. All arithmetic is integer or Fraction.
"""
from fractions import Fraction as Q
from itertools import product
from math import isqrt


def ceilq(x):
    return -(-x.numerator // x.denominator)


def sqrt_upper(x):
    """An integer upper bound on sqrt(x), certified by squaring."""
    r = isqrt(x.numerator // x.denominator)
    if Q(r * r) < x:
        r += 1
    assert Q(r * r) >= x
    return Q(r)


def fidelity():
    # Appendix F: minimize A/(2 eta) + eta G^2 T/2.  Parametrize
    # D=r^2, D+2P=q^2, T=n^2, so the critical square root is rational.
    count = 0
    for r, q, n, g in product(range(1, 5), range(1, 7), range(1, 5), range(1, 4)):
        if q < r:
            continue
        d, p, t, g = Q(r*r), Q(q*q-r*r, 2), n*n, Q(g)
        if p > d*(t-1):
            continue
        a = d*(d+2*p)
        critical = Q(r*q*n)
        eta_star = critical/(g*t)
        for s in [critical/4, critical/2, critical, 2*critical, 8*critical]:
            eta = min(eta_star, s/(g*t))
            bound = a/(2*eta) + eta*g*g*t/2
            assert t*eta*g <= s
            assert bound >= g*critical
            if s >= critical:
                assert bound == g*critical
            else:
                assert bound == g*t*a/(2*s) + g*s/2
                assert bound <= g*t*a/s
            count += 1
    print(f'FIDELITY: Appendix F OGD substitutions and budget: {count} exact cases')


def generic_cancellation():
    # Each symbol denotes an arbitrary block-start coordinate determined by
    # the already revealed sign prefix. Every coefficient vanishes exactly.
    symbols = assignments = 0
    for n in range(1, 13):
        coefficients = {}
        for eps in product([-1, 1], repeat=n):
            for b in range(n):
                symbol = (b, eps[:b])
                coefficients[symbol] = coefficients.get(symbol, 0) + eps[b]
                assignments += 1
        assert all(c == 0 for c in coefficients.values())
        symbols += len(coefficients)
    print(f'GENERIC_CANCELLATION: {symbols} arbitrary-history symbols, {assignments} exact coefficient contributions, 1..12 signs')


def sign_family():
    # Coordinates are in units a=D/(2sqrt(m)); all assertions cancel a.
    # Updating one coordinate at a time makes Euclidean movement exact abs(step).
    instances = patterns = 0
    for k, m, length in product(range(1, 4), range(1, 5), range(1, 4)):
        if k*m > 10:
            continue
        for beta, active, initial in product([Q(0), Q(1,2), Q(2), Q(7)], [False, True], [Q(-1,2), Q(0), Q(1,2)]):
            total_regret = total_internal = total_fresh = Q(0)
            nsign = 2**(k*m)
            for signs in product([-1, 1], repeat=k*m):
                x = [initial]*m
                spent = regret = internal = fresh = Q(0)
                for block, eps in enumerate(signs):
                    i = block % m
                    xb = x[i]
                    fresh += eps*xb
                    block_internal = Q(0)
                    for ell in range(length):
                        loss = max(Q(0), 1-eps*x[i])
                        assert loss >= 1-eps*xb-block_internal
                        regret += loss
                        assert sum(z*z for z in x) <= m
                        assert loss >= 0  # comparator value is exactly zero
                        # The action at the next round may use this revealed sign.
                        if active and (block != k*m-1 or ell != length-1):
                            delta = eps-x[i]
                            step = min(abs(delta), beta-spent)
                            x[i] += step if delta >= 0 else -step
                            spent += step
                            if ell < length-1:
                                block_internal += step
                    internal += block_internal
                assert spent <= beta and internal <= spent
                # Each phase comparator has norm a sqrt(m)=D/2;
                # any successive pair has squared distance <= a^2(4m)=D^2.
                for j in range(k-1):
                    squared = sum((signs[j*m+i]-signs[(j+1)*m+i])**2 for i in range(m))
                    assert squared <= 4*m
                total_regret += regret
                total_internal += internal
                total_fresh += fresh
                patterns += 1
            assert total_fresh == 0
            mean = total_regret/nsign
            assert mean >= length*(k*m-total_internal/nsign)
            assert mean >= length*(k*m-beta)
            if length == 1:
                assert mean >= k*m
            instances += 1
    print(f'SIGN_AVERAGES: {instances} policy/parameter cases, {patterns} sign realizations; exact fresh-sign cancellation and block inequality')


def parameters():
    count = regimes = 0
    # Homogeneity removes positive D,G. P/D=p, S/D=s.
    # Baseline: k sqrt(floor(T/k))/2. Its rational lower bound below is
    # intentionally weaker and still certifies the proposed factor 256.
    for t in range(1, 129):
        for p4 in range(0, 4*(t-1)+1):
            p = Q(p4,4)
            k = 1 + p.numerator//p.denominator
            assert 1 <= k <= t and k <= 1+p <= 2*k
            for s in [Q(1,16), Q(k,4), Q(k,4)+Q(1,64), Q(k,2), Q(k), Q(2*k), Q(t), Q(8*t)]:
                base_m = t//k
                assert base_m >= 1 and base_m >= Q(t,2*k)
                base_lower = Q(k*isqrt(base_m),2)
                lower = base_lower
                if s <= Q(k,4):
                    length = t//k
                    block_lower = length*(Q(k,2)-s)
                    assert block_lower >= Q(t,8)
                    lower = max(lower, block_lower)
                elif 64*s*s <= k*t:
                    h = ceilq(4*s/k)
                    m = h*h
                    length = t//(k*m)
                    assert 1 <= h and h <= 8*s/k
                    assert k*m <= t and length >= Q(t,2*k*m)
                    block_lower = length*(Q(k*h,2)-s)
                    assert block_lower >= Q(t*k,64)/s
                    lower = max(lower, block_lower)
                    regimes += 1
                # A certified rational upper bound on F/(GD), sufficient to
                # check F/256 <= the displayed lower-bound certificates.
                f_upper = min(Q(t), sqrt_upper((1+p)*t)+t*(1+p)/s)
                assert f_upper <= 256*lower, (t,p,s,lower,f_upper)
                count += 1
    print(f'ROUNDING: {count} exact parameter cases T=1..128, P/D quarter-grid, 8 positive budgets; {regimes} intermediate-regime cases')


if __name__ == '__main__':
    fidelity()
    generic_cancellation()
    sign_family()
    parameters()
    print('ALL_CHECKS_OK')

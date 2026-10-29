#!/usr/bin/env python3
"""Exact checks for the critical halo sum in arXiv:2610.00133v1.

All evidence uses integers and Fraction. No floating point, approximate
special-function evaluator, or early-small-term truncation is used.
The finite enclosures test the proposed asymptotic; they do not prove a limit.
"""
from fractions import Fraction as Q
from math import factorial, isqrt, prod

GRID = 10**45
PRINT_GRID = 10**12


def down(x, scale=GRID):
    return Q(x.numerator * scale // x.denominator, scale)


def up(x, scale=GRID):
    return -down(-x, scale)


def log_bounds(y, n=100):
    """For y>1: log y = 2 sum z^(2r+1)/(2r+1), z=(y-1)/(y+1).

    Bound the positive remainder by replacing every denominator by 2n+1
    and summing the resulting geometric series. All operations are exact.
    """
    assert y > 1
    z = (y - 1) / (y + 1)
    lower = 2 * sum((z**(2*r+1) / (2*r+1) for r in range(n)), Q(0))
    upper = lower + 2*z**(2*n+1) / ((2*n+1)*(1-z*z))
    return lower, upper


def root_bounds():
    lo, hi = Q(3, 8), Q(2, 5)
    assert log_bounds(1+3*lo)[0] > 2*lo
    assert log_bounds(1+3*hi)[1] < 2*hi
    # h'(p)=(1-6p)/(1+3p)<0 throughout this bracket.
    for _ in range(100):
        mid = (lo+hi)/2
        a, b = log_bounds(1+3*mid)
        if a > 2*mid:
            lo = mid
        else:
            assert b < 2*mid
            hi = mid
    assert Q(381344280425169490, 10**18) < lo < hi < Q(381344280425169492, 10**18)
    return lo, hi


def sqrt_bounds(x):
    n = isqrt(x.numerator * GRID**2 // x.denominator)
    lo, hi = Q(n, GRID), Q(n+1, GRID)
    assert lo*lo <= x < hi*hi
    return lo, hi


def lambda_squared(p):
    # At the critical root, exp(2p)=1+3p, so this equals exp(6p)/(27p).
    return (1+3*p)**3 / (27*p)


def constant_squared(p):
    return 3*(1+3*p)**2 / (6*p-1)


def ratio(K, q, j):
    return (q*K*K*(3*j+1)*(3*j+2)*(3*j+3) /
            ((j+1)**2*(K+3*j)*(K+3*j+1)*(K+3*j+2)))


def coefficient(K, q, j):
    rising = prod((K+r for r in range(3*j)), start=Q(1))
    return q**j * K**(2*j) * factorial(3*j) / (factorial(j)**2 * rising)


def fidelity():
    count = 0
    for K in (Q(3, 2), Q(2), Q(7, 2), Q(4), Q(8), Q(16)):
        for q in (Q(1, 4), Q(1), Q(4)):
            for j in range(1, 21):
                u, v = coefficient(K, q, j), coefficient(K, q, j+1)
                assert v/u == ratio(K, q, j)
                assert ratio(K, q, j) <= 27*q*(j+1)/K
                assert ratio(K, q, j) <= 2*q*K*K/(j*j)
                count += 1
    for K in (Q(3, 2), Q(2), Q(7, 2), Q(4), Q(8), Q(16)):
        for m in range(31):
            w = K**m / prod((K+r for r in range(m)), start=Q(1))
            product_form = prod((K/(K+r) for r in range(m)), start=Q(1))
            assert w == product_form
    print(f'FIDELITY: source weight product and Appendix D.1 ratio; {count} exact ratio cases OK')


def prefactor_identity():
    # a(p)^2*2*pi/A = C(p)^2, with pi cancelled exactly.
    for p in (Q(1, 4), Q(3, 8), Q(2, 5), Q(1, 2)):
        A = (6*p-1)/(p*(1+3*p))
        a_squared_times_two_pi = 3*(1+3*p)/p
        assert a_squared_times_two_pi/A == constant_squared(p)
    print('PREFACTOR: a(p)^2 * (2*pi/A) = 3*(1+3*p)^2/(6*p-1), exact rational specializations OK')


def sum_bounds(K, qlo, qhi):
    """Enclose the entire sum, including both humps and its infinite tail.

    Outward rounding at each positive recurrence step preserves inclusion.
    Beyond J, ratio <= 2*qhi*K^2/J^2 = rho < 1; the tail after J is at most
    u_J*rho/(1-rho). No decision to stop depends on a term getting small.
    """
    K = Q(K)
    J = 8*int(K)+100
    rho = 2*qhi*K*K/(J*J)
    assert 0 < rho < Q(1, 2)
    ulo = down(6*qlo*K/((K+1)*(K+2)))
    uhi = up(6*qhi*K/((K+1)*(K+2)))
    slo, shi = ulo, uhi
    for j in range(1, J):
        ulo = down(ulo*ratio(K, qlo, j))
        uhi = up(uhi*ratio(K, qhi, j))
        slo += ulo
        shi += uhi
    shi += uhi*rho/(1-rho)
    assert slo <= shi
    return slo, shi, J


def interval(label, lo, hi):
    a, b = down(lo, PRINT_GRID), up(hi, PRINT_GRID)
    print(f'{label}: [{a}, {b}]')


def main():
    fidelity()
    prefactor_identity()
    plo, phi = root_bounds()
    # R(p) is increasing for p>1/6. C(p)^2 is decreasing for 1/6<p<2/3.
    qlo, qhi = down(lambda_squared(plo)), up(lambda_squared(phi))
    clo = sqrt_bounds(constant_squared(phi))[0]
    chi = sqrt_bounds(constant_squared(plo))[1]
    assert Q(3272073540231, 10**12) < clo < chi < Q(3272073540232, 10**12)
    interval('critical root p', plo, phi)
    interval('critical lambda squared', qlo, qhi)
    interval('proposed limit C', clo, chi)
    for K in (Q(3, 2), 2, 3, Q(7, 2), 4, 8, 16, 32, 64, 128, 256, 512, Q(1025, 2), 1024):
        lo, hi, J = sum_bounds(K, qlo, qhi)
        assert hi-lo < Q(1, 10**20)
        interval(f'S_K K={K} cutoff={J}', lo, hi)
        error_lower, error_upper = lo-chi, hi-clo
        if K == 2:
            assert error_upper < 0
        if K == 1024:
            assert max(abs(error_lower), abs(error_upper)) < Q(1, 20)
    print('EXPECTED SETTLEMENT CHECK: exact finite enclosures; K=2 is below C and K=1024 is within 1/20 of C')
    print('LIMIT SCOPE: finite checks do not prove convergence for real K tending to infinity')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

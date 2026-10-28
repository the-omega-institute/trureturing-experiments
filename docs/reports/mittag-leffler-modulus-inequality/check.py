#!/usr/bin/env python3
"""Exact certificates for the GGPS Conjecture 25 endpoint obstruction.

Only integers and fractions are used. The analytic family proof is in README.md;
this program certifies actual series values, independently of asymptotic input.
"""
from fractions import Fraction as Q
from math import factorial, isqrt


def atan_bounds(q, terms=160):
    """Alternating-series enclosure for atan(1/q), q a positive integer."""
    s = sum((Q((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
             for k in range(terms)), Q(0))
    other = s + Q((-1) ** terms, (2 * terms + 1) * q ** (2 * terms + 1))
    return min(s, other), max(s, other)


def rho_bounds():
    """Machin's identity and integer square roots enclose 1/sqrt(pi)."""
    a, b = atan_bounds(5), atan_bounds(239)
    pi_lo, pi_hi = 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]
    assert Q(3) < pi_lo < pi_hi < Q(22, 7)
    scale = 10 ** 180
    lo = Q(isqrt(scale * scale * pi_hi.denominator // pi_hi.numerator), scale)
    hi = Q(isqrt(scale * scale * pi_lo.denominator // pi_lo.numerator) + 1, scale)
    assert 0 < lo < hi < 1
    assert lo * lo * pi_hi <= 1 <= hi * hi * pi_lo
    return lo, hi


RHO = rho_bounds()


def rgamma_half_integer(s):
    """Return (a,b) with 1/Gamma(s)=a+b/sqrt(pi), for positive half-integer s."""
    assert s > 0 and (2 * s).denominator == 1
    if s.denominator == 1:
        return Q(1, factorial(int(s) - 1)), Q(0)
    k = int(s - Q(1, 2))
    return Q(0), Q(4 ** k * factorial(k), factorial(2 * k))


def cmul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def affine_rho(a, b):
    x, y = a + b * RHO[0], a + b * RHO[1]
    return min(x, y), max(x, y)


def interval_ml_endpoint(z, cutoff=1200):
    """Rectangular enclosure for E_(3/2,1/2)(z), including the infinite tail.

    Terms two places apart have quotient z^2/[s(s+1)(s+2)],
    s=(3*n+1)/2. Each parity tail is bounded by a geometric series,
    using |z| <= |Re z|+|Im z| and 1/sqrt(pi) < 1.
    """
    z = tuple(map(Q, z))
    rational, half = [Q(0), Q(0)], [Q(0), Q(0)]
    power = (Q(1), Q(0))
    for n in range(cutoff):
        a, b = rgamma_half_integer(Q(3 * n + 1, 2))
        for j in range(2):
            rational[j] += a * power[j]
            half[j] += b * power[j]
        power = cmul(power, z)
    radius = abs(z[0]) + abs(z[1])
    error = Q(0)
    for n in [cutoff, cutoff + 1]:
        s = Q(3 * n + 1, 2)
        a, b = rgamma_half_integer(s)
        ratio = radius * radius / (s * (s + 1) * (s + 2))
        assert 0 <= ratio < 1
        error += (a + b) * radius ** n / (1 - ratio)
    out = []
    for j in range(2):
        lo, hi = affine_rho(rational[j], half[j])
        out.append((lo - error, hi + error))
    return out, error


def outer_grid(interval, scale=10 ** 12):
    """Coarsen an interval outwards to short rational endpoints."""
    lo, hi = interval
    return Q(lo.numerator * scale // lo.denominator, scale), Q(
        -((-hi.numerator * scale) // hi.denominator), scale)


def max_abs(interval):
    return max(abs(interval[0]), abs(interval[1]))


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def main():
    # Fidelity: the source states E_(1,1)=exp and E_(2,2)=sinh(sqrt(z))/sqrt(z).
    for n in range(81):
        assert rgamma_half_integer(Q(n + 1)) == (Q(1, factorial(n)), 0)
        assert rgamma_half_integer(Q(2 * n + 2)) == (Q(1, factorial(2 * n + 1)), 0)
    print('FIDELITY: source exponential and hyperbolic-sine coefficients agree for n=0..80')

    # Independent algebraic checks of the family mechanism, with u=alpha-1>0.
    assert poly_mul(poly_mul([0, 1], [1, 1]), [2, 1]) == [0, 2, 3, 1]
    assert all(a > 0 for a in [2, 3, 1])
    # Thus alpha(alpha-1)(alpha+1)=u(u+1)(u+2)>0.
    # At alpha=3/2, c=15/(8*sqrt(pi)); the recurrence product is exact.
    alpha = Q(3, 2)
    assert alpha * (alpha - 1) * (alpha + 1) == Q(15, 8)
    w = (Q(-1), Q(1, 2))
    assert w[0] ** 2 + w[1] ** 2 == Q(5, 4)
    assert 1 / (w[0] ** 2 + w[1] ** 2) == Q(4, 5) < 1
    # 1-(1+t^2)^(-2) has positive numerator 2t^2+t^4 for every t>0.
    assert poly_mul([1, 0, 1], [1, 0, 1]) == [1, 0, 2, 0, 1]
    print('FAMILY ALGEBRA: positive normalized inverse-square coefficient; ray ratio 4/5')

    for r in [1024, 1536, 2048]:
        real_box, err_r = interval_ml_endpoint((-r, 0))
        complex_box, err_c = interval_ml_endpoint((-r, Q(r, 2)))
        re_real = outer_grid(real_box[0])
        re_c = outer_grid(complex_box[0])
        im_c = outer_grid(complex_box[1])
        lower = re_real[0]
        norm_sq_upper = max_abs(re_c) ** 2 + max_abs(im_c) ** 2
        assert lower > 0, ('real-axis positivity failed', r)
        assert norm_sq_upper < lower ** 2, ('strict modulus failure not certified', r)
        assert err_r < Q(1, 10 ** 60) and err_c < Q(1, 10 ** 60)
        print(f'WITNESS alpha=3/2 beta=1/2 z=-{r}+{Q(r,2)}*i')
        print(f'  E(-r) in [{re_real[0]}, {re_real[1]}]')
        print(f'  Re E(z) in [{re_c[0]}, {re_c[1]}]')
        print(f'  Im E(z) in [{im_c[0]}, {im_c[1]}]')
        print(f'  |E(z)|^2 <= {norm_sq_upper} < {lower ** 2} <= E(-r)^2')
        print('  each parity-tail bound < 1/10^60; strict violation certified')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

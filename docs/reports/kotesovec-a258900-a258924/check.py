#!/usr/bin/env python3
"""Exact coefficient and normalization checks for OEIS A258900 and A258924.

Python standard library only. No floating-point calculation is used. The checks
are finite identities supporting the proposed analytic argument, not a proof of
an asymptotic limit. Source: the OEIS internal records, revisions 18 and 11.
"""
from fractions import Fraction as Q
from math import factorial


def mul(a, b, n):
    c = [Q(0) for _ in range(n + 1)]
    for i, ai in enumerate(a[:n + 1]):
        for j, bj in enumerate(b[:n + 1 - i]):
            c[i + j] += ai * bj
    return c


def power(a, k, n):
    assert isinstance(k, int) and k >= 0
    r = [Q(1)] + [Q(0)] * n
    while k:
        if k % 2:
            r = mul(r, a, n)
        a = mul(a, a, n)
        k //= 2
    return r


def unit_power(a, exponent, n):
    """Truncated binomial expansion of a(z)^exponent when a(0)=1."""
    assert a[0] == 1
    h = (a + [Q(0)] * (n + 1))[:n + 1]
    h[0] = 0
    term = [Q(1)] + [Q(0)] * n
    result = term.copy()
    binom = Q(1)
    for k in range(1, n + 1):
        term = mul(term, h, n)
        binom *= Q(exponent - k + 1, k)
        result = [u + binom * v for u, v in zip(result, term)]
    return result


def source_primitive(p, n):
    """F(y)=y H(y^p), H_k=(1/p)_k/(k! (pk+1))."""
    h, rising = [], Q(1)
    for k in range(n + 1):
        if k:
            rising *= Q(1, p) + k - 1
        h.append(rising / (factorial(k) * (p * k + 1)))
    return h


def from_reversion(p, n):
    h = source_primitive(p, n)
    # Lagrange: [x^m]S = [u^(m-1)](u/F(u))^m / m.
    return [factorial(p * k) * unit_power(h, -(p * k + 1), k)[k]
            for k in range(n + 1)]


def from_tangent_ode(p, n):
    # T=x q(x^p), T'=1+T^p, q(0)=1.
    q = [Q(1)] + [Q(0)] * n
    for k in range(1, n + 1):
        q[k] = power(q, p, k - 1)[k - 1] / (p * k + 1)
    qp = power(q, p, n)
    c = unit_power([Q(1)] + qp[:n], -Q(1, p), n)
    return [factorial(p * k) * c[k] for k in range(n + 1)], q, c


def check_series(p, n):
    a = from_reversion(p, n)
    b, tangent, c = from_tangent_ode(p, n)
    assert a == b, (p, 'independent constructions disagree')
    assert all(v.denominator == 1 for v in a)
    s = [a[k] / factorial(p * k + 1) for k in range(n + 1)]
    # C^p + S^p = C(z)^p + z (S(x)/x)^p = 1, z=x^p.
    cp, sp = power(c, p, n), power(s, p, n)
    assert [cp[k] + (sp[k - 1] if k else 0) for k in range(n + 1)] == [1] + [0] * n
    assert mul(c, tangent, n) == s  # T=S/C.
    # Check F(S(x))=x directly, independently of Lagrange's formula.
    composition = [Q(0)] * (n + 1)
    for k, hk in enumerate(source_primitive(p, n)):
        sk = power(s, p * k + 1, n - k)
        for j, sj in enumerate(sk):
            composition[k + j] += hk * sj
    assert composition == [1] + [0] * n
    for k in range(n + 1):
        assert (p * k + 1) * s[k] == c[k]  # S'=C.
    assert all(t > 0 for t in tangent)
    return a, tangent


PUBLISHED = {
    4: [1, -6, -1764, -7700616, -147910405104, -8310698364852576,
        -1085420895640591777344, -284168646775526186095019136,
        -134459287943928269154814258953984,
        -106506405136317713669903020280294647296],
    5: [1, -24, -169344, -25255286784, -23089632627769344,
        -79051067969864491597824, -766667475511149432871084621824,
        -17578325209217134578862801556544159744,
        -839197248407269659950832532302025663168118784],
}
PUBLISHED_TANGENT = {
    4: [1, 24, 32256, 285272064, 8967114326016],
    5: [1, 120, 3024000, 858574080000, 1226178516326400000],
}


def exponent_sum(*terms):
    """Exact bookkeeping for positive monomials after named reflection identities."""
    result = {}
    for d in terms:
        for k, v in d.items():
            result[k] = result.get(k, Q(0)) + Q(v)
    return {k: v for k, v in result.items() if v}


def check_normalization():
    # K_p=A sqrt(p)/((p-1) Gamma(1-alpha) (2*pi)^((p-1)/2)).
    # For p=4, sin(pi/4)=2^(-1/2), and
    # Gamma(2/3)^(-1)=Gamma(1/3)*sqrt(3)/(2*pi).
    k4 = exponent_sum(
        {'2': -Q(7, 6), '3': Q(1, 3), 'pi': Q(1, 3)},  # A
        {'2': 1}, {'3': -1}, {'2': -Q(3, 2), 'pi': -Q(3, 2)},
        {'Gamma(1/3)': 1, '3': Q(1, 2), '2': -1, 'pi': -1})
    assert k4 == {'Gamma(1/3)': 1, '2': -Q(8, 3), '3': -Q(1, 6), 'pi': -Q(13, 6)}
    # For p=5, Gamma(3/4)^(-1)=Gamma(1/4)/(pi*sqrt(2)), and
    # sin(pi/5)=2^(-3/2)*(5-sqrt(5))^(1/2).
    k5 = exponent_sum(
        {'2': Q(1, 2), '5': -Q(1, 2), 'pi': Q(1, 4), 'sin(pi/5)': -Q(1, 4)},
        {'5': Q(1, 2)}, {'2': -2}, {'2': -2, 'pi': -2},
        {'Gamma(1/4)': 1, 'pi': -1, '2': -Q(1, 2)},
        {'sin(pi/5)': Q(1, 4), '2': Q(3, 8), '5-sqrt(5)': -Q(1, 8)})
    assert k5 == {'Gamma(1/4)': 1, '2': -Q(29, 8), '5-sqrt(5)': -Q(1, 8), 'pi': -Q(11, 4)}
    assert Q(4 + 1, 2) + Q(1, 4 - 1) == Q(17, 6)
    assert Q(5 + 1, 2) + Q(1, 5 - 1) == Q(13, 4)
    assert 2 * 4 * 2 - Q(4, 2) == 14
    # The positive sine root is independently checked algebraically in Q(sqrt(5)):
    # cos(2*pi/5)=(sqrt(5)-1)/4; (1-cos(2*pi/5))/2=(5-sqrt(5))/8.
    assert (Q(1) + Q(1, 4)) / 2 == Q(5, 8)
    assert -Q(1, 4) / 2 == -Q(1, 8)
    print('NORMALIZATION: exact rational exponents agree with both conjectured constants')
    print('NAMED INPUTS: gamma reflection and the positive elementary sine values')
    print('ASYMPTOTIC POWERS: p=4 -> 17/6, p=5 -> 13/4; p=4 growth -> 2^(14*n)')


def main():
    for p in range(3, 9):
        n = 12 if p in PUBLISHED else 8
        a, t = check_series(p, n)
        if p in PUBLISHED:
            assert a[:len(PUBLISHED[p])] == PUBLISHED[p]
            assert [t[k] * factorial(p * k + 1) for k in range(5)] == PUBLISHED_TANGENT[p]
            print(f'FIDELITY p={p}: {len(PUBLISHED[p])} published a(n) terms and 5 published S/C terms agree')
            print(f'PREFIX p={p}: ' + ', '.join(str(int(v)) for v in a[:6]))
        print(f'IDENTITIES p={p}, 0<=n<={n}: independent constructions, F(S)=x, S\'=C, C^p+S^p=1, T=S/C, T\'=1+T^p OK')
    check_normalization()
    print('SCOPE: finite exact checks; analytic continuation, transfer, and the limit are not machine-proved')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

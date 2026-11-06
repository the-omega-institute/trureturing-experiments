#!/usr/bin/env python3
"""Exact arithmetic for Yang Remark 5.2; this does not prove the global upper bound.

Source: arXiv:2610.13103v1, Theorem 5.1 (fidelity) and Remark 5.2.
The T4/T5 parametrizations agree with Jiang--Yang (2024), Lemma 5.2.
Only integers and exact symbolic algebra occur.
"""
from math import gcd, isqrt
import sympy as sp


def ceil_two_sqrt(m):
    r = isqrt(4 * m)
    return r if r * r == 4 * m else r + 1


def feasible_a(m, s):
    """All positive a with a(s-a) >= m, obtained using an integer square root."""
    disc = s * s - 4 * m
    if disc < 0:
        return range(0)
    r = isqrt(disc)
    return range(max(1, (s - r + 1) // 2), min(s - 1, (s + r) // 2) + 1)


def check_B_gap(n):
    A = ceil_two_sqrt(n + 1)
    assert not any(gcd(a, n) == 1 for a in feasible_a(n + 1, A))
    assert any(gcd(a, n) == 1 for a in feasible_a(n + 1, A + 1))
    return A, A + 1


def check_sequence(n, x, a, b, required_e, required_xe):
    assert 1 < x < n and gcd(x, n) == 1
    assert 0 < a < n and 0 < b < n  # both support elements have order n
    # Every subsequence is e^i (xe)^j, with i <= a and j <= b.
    # Test the only two boundary faces on which support can be exhausted.
    assert (-b * x) % n == required_e > a
    assert (-a * pow(x, -1, n)) % n == required_xe > b


def fidelity():
    # Literal first instances of source Theorem 5.1, its two displayed families.
    for label, n, d, x, req_e, req_xe in [
        ('T2', 2730 * (2730 - 2835 + 736), 1313, None, None, None),
        ('T3', 2730 * (2730 - 2625 + 631), 1418, None, None, None),
    ]:
        assert n == d * (d - 1) - 26
        if label == 'T2':
            assert (12 * n - d) % 13 == 0
            x, req_e, req_xe = (12 * n - d) // 13, n - 2, n - 13
        else:
            assert (n - d) % 2 == 0
            x, req_e, req_xe = (n - d) // 2, n - 13, n - 2
        A, B = check_B_gap(n)
        assert A == 2 * d - 1 and B == 2 * d
        check_sequence(n, x, n - d, n - d + 1, req_e, req_xe)
        print(f'FIDELITY {label}: n={n}, d={d}, x={x}, A={A}, B={B}, '
              f'length={2*n-A}, source_N={2*n-A+1}')


def symbolic():
    t = sp.symbols('t', integer=True)
    for name, d in [('T4', 30030 * t + 2192), ('T5', 30030 * t - 2192)]:
        n = sp.expand(d*d - 64)
        sign = 1 if name == 'T4' else -1
        assert sp.expand(n - 30030*(30030*t*t + sign*4384*t + 160)) == 0
        r = d - 1 if name == 'T4' else d + 1
        multiplier = d + 1 if name == 'T4' else d - 1
        c = 4 if name == 'T4' else 3
        x = sp.expand((c * n - r) / 7)
        k = sp.expand((c * multiplier - 1) / 7)
        assert sp.Poly(x, t).domain == sp.ZZ
        assert sp.Poly(k, t).domain == sp.ZZ
        assert sp.expand(7*x + r - c*n) == 0
        assert sp.expand(multiplier*x + 9 - k*n) == 0
        # Any common divisor of x,n divides r and 9; r is nonzero modulo 3.
        assert sp.Poly(r, t).all_coeffs()[-1] % 3 != 0
        print(f'SYMBOLIC {name}: d={d}; n={n}; x={x}; k={k}')
        print(f'IDENTITIES {name}: 7*x + ({r}) = {c}*n; '
              f'({multiplier})*x + 9 = ({k})*n')


def candidate():
    first = {
        'T4': (0, 4804800, 2192, 2745287, 4802609, 4802607, 9605217),
        'T5': (1, 774954180, 27838, 332119243, 774926341, 774926343, 1549852685),
    }
    for name, first_t in [('T4', 0), ('T5', 1)]:
        for t in range(first_t, 1001):
            d = 30030*t + (2192 if name == 'T4' else -2192)
            n = d*d - 64
            c = 4 if name == 'T4' else 3
            r = d - 1 if name == 'T4' else d + 1
            multiplier = d + 1 if name == 'T4' else d - 1
            assert (c*n - r) % 7 == 0
            x = (c*n - r)//7
            assert (c*multiplier - 1) % 7 == 0
            k = (c*multiplier - 1)//7
            assert 7*x + r == c*n and multiplier*x + 9 == k*n
            assert gcd(r, 9) == 1 and gcd(x, n) == 1
            assert (2*d-1)**2 < 4*(n+1) <= (2*d)**2
            assert ceil_two_sqrt(n+1) == 2*d
            a = n-d+1 if name == 'T4' else n-d-1
            b = n-d-1 if name == 'T4' else n-d+1
            check_sequence(n, x, a, b, n-9, n-7)
            assert a+b == 2*n-2*d
            if t == first_t:
                row = (t, n, d, x, a, b, a+b+1)
                assert row == first[name]
                assert check_B_gap(n) == (2*d, 2*d+1)
                print(f'FIRST {name}: t={t}, n={n}, d={d}, x={x}, '
                      f'multiplicities=({a},{b}), lower_bound_N={a+b+1}, '
                      f'published_upper_bound={2*n-2*d+1}')
        print(f'FAMILY {name}: exact checks for every integer t in [{first_t},1000]')


if __name__ == '__main__':
    fidelity()
    symbolic()
    candidate()
    print('ALL_CHECKS_OK')

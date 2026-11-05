#!/usr/bin/env python3
"""Exact finite checks for the literal, normalized Girard-Waring coefficients.

Source: Wojnar, Wojnar and Brin, arXiv:1706.08381v2, equation (2),
and OEIS A301424 revision 26, with the A301417 extended definition.
These computations do not establish the all-index generating function.
Python standard library only; every numerical calculation is exact.
"""
from fractions import Fraction as Q
from math import comb, factorial, prod

PUBLISHED = [
    1, 7, 64, 609, 5846, 56161, 539540, 5183417, 49797685,
    478412117, 4596160548, 44155846113, 424210322004, 4075437640457,
    39153200900024, 376149330687809, 3613710136705565,
    34717331354145139, 333533418773956668, 3204294140706218329,
    30784024515164777522,
]


def vectors(d, k):
    """Enumerate each m in N_0^d of weighted degree k exactly once."""
    def visit(i, left, prefix):
        if i == d:
            if left % i == 0:
                yield tuple(prefix + [left // i])
            return
        for v in range(left // i + 1):
            yield from visit(i + 1, left - i * v, prefix + [v])
    yield from visit(1, k, [])


def literal(d, k):
    assert d >= 1 and k >= 1
    result = {}
    for m in vectors(d, k):
        r = sum(m)
        assert 1 <= r <= k
        c = (Q(k, d * r) * (-1) ** (k + r)
             * Q(factorial(r), prod(factorial(v) for v in m))
             * prod(comb(d, i + 1) ** v for i, v in enumerate(m)))
        assert c != 0
        result[m] = c
    return result


def newton(d, bound):
    """Independent polynomial recursion, with e_i=binomial(d,i)*t_i.

    p_k = sum_{i=1}^{k-1} (-1)^(i-1)e_i p_(k-i)
          + (-1)^(k-1)k e_k, with e_i=0 for i>d; T_k=p_k/d.
    No factorial or partition enumeration occurs here.
    """
    powers = [{}]
    for k in range(1, bound + 1):
        p = {}
        for i in range(1, min(d, k - 1) + 1):
            for m, c in powers[k - i].items():
                v = list(m)
                v[i - 1] += 1
                v = tuple(v)
                p[v] = p.get(v, Q(0)) + (-1) ** (i - 1) * comb(d, i) * c
        if k <= d:
            v = tuple(int(i == k) for i in range(1, d + 1))
            p[v] = p.get(v, Q(0)) + (-1) ** (k - 1) * k * comb(d, k)
        powers.append({m: c for m, c in p.items() if c})
    return [{m: c / d for m, c in p.items()} for p in powers]


def convolve(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def gf_coefficients(bound):
    # Coefficients in increasing degree order; denominator has constant 1.
    numerator = [1] + [-comb(6, i) for i in range(7)]
    denominator = convolve([1, -1], [1] + [-comb(7, i) for i in range(1, 8)])
    printed = [1, -8, -14, -14, 0, 14, 14, 6, 1]
    assert denominator == printed
    a = []
    for n in range(bound):
        a.append(Q(numerator[n] if n < len(numerator) else 0)
                 - sum(denominator[j] * a[n - j]
                       for j in range(1, min(n, 8) + 1)))
    return a, denominator


def masses(poly):
    pos = sum((c for c in poly.values() if c > 0), Q(0))
    neg_magnitude = -sum((c for c in poly.values() if c < 0), Q(0))
    return pos, neg_magnitude, sum(poly.values(), Q(0)), sum(map(abs, poly.values()), Q(0))


def elementary(data):
    e = [Q(1)] + [Q(0)] * len(data)
    for x in data:
        for i in range(len(data), 0, -1):
            e[i] += x * e[i - 1]
    return e


def main():
    # Paper fidelity: equation (1), plus the n=3 cubic row of Table GW.n=3.
    assert literal(3, 2) == {(2, 0, 0): Q(3), (0, 1, 0): Q(-2)}
    assert literal(3, 3) == {(3, 0, 0): Q(9), (1, 1, 0): Q(-9), (0, 0, 1): Q(1)}
    assert [masses(literal(6, k))[0] for k in range(1, 8)] == [1, 6, 46, 371, 3026, 24707, 201748]
    print('PAPER_FIDELITY_OK: n=3 quadratic/cubic dictionaries; n=6 table masses k=1..7')
    bound = 27
    newton_polys = newton(7, bound)
    gf, denominator = gf_coefficients(bound)
    print('GF_DENOMINATOR_INCREASING_DEGREE:', denominator)
    exact_polys = []
    print('k monomials positive_mass negative_magnitude signed_mass')
    for k in range(1, bound + 1):
        p = literal(7, k)
        assert p == newton_polys[k], ('Newton dictionary', k)
        pos, neg, signed, unsigned = masses(p)
        assert pos == gf[k - 1], ('GF', k)
        assert signed == 1 and neg == pos - 1 and 2 * pos == unsigned + 1
        assert pos.denominator == 1
        if k <= len(PUBLISHED):
            assert pos == PUBLISHED[k - 1], ('OEIS DATA', k)
        print(k, len(p), pos, neg, signed)
        exact_polys.append(p)
    print('OEIS_FIDELITY_OK: all 21 published A301424 DATA terms')
    data_sets = [
        [Q(0)] * 7,
        [Q(1)] * 7,
        list(map(Q, [-3, -1, 0, 0, 2, 2, 5])),
        [Q(1, 2), Q(-2, 3), Q(3, 4), Q(0), Q(5, 7), Q(-1), Q(2)],
    ]
    for data in data_sets:
        e = elementary(data)
        t = [e[i] / comb(7, i) for i in range(1, 8)]
        for k, p in enumerate(exact_polys, 1):
            evaluated = sum((c * prod(t[i] ** m[i] for i in range(7))
                             for m, c in p.items()), Q(0))
            assert evaluated == sum((x ** k for x in data), Q(0)) / 7
    assert len(exact_polys[-1]) == 1175
    assert gf[26] == 24203797992661339086776056
    print('DATA_EVALUATION_OK: 4 exact seven-data families, k=1..27, including repeated/zero/rational data')
    print('INDEPENDENT_METHODS_OK: factorial formula, Newton dictionaries, exact rational GF division')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

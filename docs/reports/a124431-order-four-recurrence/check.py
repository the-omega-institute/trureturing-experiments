#!/usr/bin/env python3
"""Exact finite checks for the literal OEIS A124431 sum and Mathar's recurrence.

The checks supply finite evidence, not an all-index proof. Only Python integers
and standard-library combinatorial arithmetic are used.
"""
from collections import defaultdict
from math import comb

LIMIT = 250
CT_LIMIT = 60
SOURCE_PREFIX = [
    1, 3, 9, 29, 97, 331, 1145, 4001, 14089, 49915, 177713,
    635293, 2278841, 8198227, 29567729, 106872961, 387038993,
    1404052659, 5101219929, 18559193245, 67605310097,
    246541193883, 899999057385, 3288522934433, 12026324883865,
]


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def literal(n):
    return sum(2**k * choose((n + k) // 2, k)
               * choose((n + k + 1) // 2, k) for k in range(n + 1))


def peakless_colored(n):
    # State (height, was_up). A down-step immediately after an up-step is forbidden.
    states = {(0, False): 1}
    for _ in range(n):
        nxt = defaultdict(int)
        for (height, was_up), count in states.items():
            nxt[height + 1, True] += count
            nxt[height, False] += (3 if height == 0 else 2) * count
            if height > 0 and not was_up:
                nxt[height - 1, False] += count
        states = nxt
    return sum(count for (height, _), count in states.items() if height == 0)


def mul(f, g, degree):
    out = [0] * (degree + 1)
    for i, a in enumerate(f[:degree + 1]):
        for j, b in enumerate(g[:degree + 1 - i]):
            out[i + j] += a * b
    return out


def add(f, g, degree, scale=1):
    return [(f[i] if i < len(f) else 0)
            + scale * (g[i] if i < len(g) else 0) for i in range(degree + 1)]


def shift(f, k, degree):
    return ([0] * k + f)[:degree + 1]


def derivative(f):
    return [i * f[i] for i in range(1, len(f))]


def constant_term_series(degree):
    # K = 1 + (2*x + x^2 + x^2*t + 2*x/t)*K, in Q[t,t^-1][[x]].
    terms = [{0: 1}]
    for n in range(1, degree + 1):
        row = defaultdict(int)
        for exponent, coefficient in terms[n - 1].items():
            row[exponent] += 2 * coefficient
            row[exponent - 1] += 2 * coefficient
        if n >= 2:
            for exponent, coefficient in terms[n - 2].items():
                row[exponent] += coefficient
                row[exponent + 1] += coefficient
        terms.append(dict(row))
    return ([row.get(0, 0) for row in terms],
            [row.get(1, 0) for row in terms],
            [row.get(-1, 0) for row in terms])


def parity_split(n):
    even = odd = 0
    for m in range(n + 1):
        k = 2 * m - n
        if 0 <= k <= m:
            even += 2**k * choose(m, k)**2
        k = 2 * m + 1 - n
        if 0 <= k <= m:
            odd += 2**k * choose(m, k) * choose(m + 1, k)
    return even + odd


def main():
    a = [literal(n) for n in range(LIMIT + 1)]
    assert a[:len(SOURCE_PREFIX)] == SOURCE_PREFIX
    print('SOURCE_PREFIX_OK: OEIS revision 41, all 25 listed terms, n=0..24')
    assert [peakless_colored(n) for n in range(13)] == a[:13]
    assert peakless_colored(3) == 3**3 + 2 == 29
    print('SOURCE_FIDELITY_OK: Deutsch peakless colored path counts, n=0..12; a(3)=27+2=29')
    print('INITIAL_VALUES:', ','.join(map(str, a[:9])))

    for n in range(4, LIMIT + 1):
        comment = ((n + 1) * a[n] - 2 * (2 * n + 1) * a[n - 1]
                   + 2 * (n - 1) * a[n - 2] + 2 * (5 - 2 * n) * a[n - 3]
                   + (n - 3) * a[n - 4])
        formula = ((n + 1) * a[n] + 2 * (-2 * n - 1) * a[n - 1]
                   + 2 * (n - 1) * a[n - 2] + 2 * (-2 * n + 5) * a[n - 3]
                   + (n - 3) * a[n - 4])
        assert comment == formula == 0, (n, comment, formula)
    print('RECURRENCE_OK: both COMMENT and FORMULA versions, n=4..250 (247 cases)')
    assert all(parity_split(n) == a[n] for n in range(101))
    print('PARITY_SPLIT_OK: n=0..100')

    degree = CT_LIMIT
    G, H, L = constant_term_series(degree)
    U, B, C, P, Q = [1, -2, -1], [0, 0, 1], [0, 2], [1, 0, 1], [1, -4, 1]
    D = mul(P, Q, 4)
    assert D == add(mul(U, U, 4), mul(B, C, 4), 4, -4) == [1, -4, 2, -4, 1]
    ct_F = add(mul([1, 1], G, degree), H, degree, 2)
    assert ct_F == a[:degree + 1]
    assert mul(B, L, degree) == mul(C, H, degree)
    assert add(mul(U, G, degree), shift(H, 1, degree), degree, -4) == [1] + [0] * degree
    assert mul(D, mul(G, G, degree), degree) == [1] + [0] * degree
    Y = add([1], shift(a, 1, degree), degree, 2)
    assert Y == mul(P, G, degree)
    assert mul(Q, mul(Y, Y, degree), degree) == P + [0] * (degree - 2)
    print('CONSTANT_TERM_AND_ALGEBRAIC_GF_OK: all coefficients x^0..x^60')

    central = [comb(2 * j, j) for j in range(101)]
    assert mul(central, central, 100) == [4**j for j in range(101)]
    assert add(mul(derivative(P), Q, 4), mul(P, derivative(Q), 4), 4, -1) == [4, 0, -4, 0, 0]
    ode = add(shift(mul(D, derivative(a), LIMIT - 1), 1, LIMIT),
              mul([1, -6, 2, -2, 1], a, LIMIT), LIMIT)
    assert ode == [1, 0, -1] + [0] * (LIMIT - 2)
    print('CENTRAL_BINOMIAL_CONVOLUTION_OK: degree=0..100')
    print('DIFFERENTIAL_OPERATOR_OK: polynomial identity exact; ODE coefficients x^0..x^250')
    print('EVIDENCE_SCOPE: finite exact checks; no unbounded theorem or Lean proof is asserted')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

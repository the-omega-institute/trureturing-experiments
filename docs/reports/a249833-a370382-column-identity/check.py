#!/usr/bin/env python3
"""Exact, independent truncations of the source ODE and Kurkov's array.

Sources: OEIS A249833 revision 30, A370382 revision 15; Paul D. Hanna
and Mikhail Kurkov. No floating-point arithmetic or external packages.
This finite consistency check is not a proof of the infinite identity.
"""
from fractions import Fraction as Q
from math import comb, factorial

MAX_DEGREE = 40
DISPLAYED_A249833 = [
    1, 1, 2, 7, 34, 210, 1574, 13866, 140340, 1604284, 20439484,
    287152488, 4409695952, 73482586464, 1320533540808, 25456195929232,
    523975944225280, 11469534961767408, 266038450202037728,
    6518167274358688512, 168209881653024622944, 4560447490191133853536,
    129593625015740116555072,
]
DISPLAYED_A370382 = [
    1, 1, 2, 1, 5, 7, 1, 13, 27, 34, 1, 33, 108, 176, 210, 1, 81,
    437, 944, 1364, 1574, 1, 193, 1758, 5154, 9144, 12292, 13866, 1,
    449, 6959, 28304, 62496, 98742, 126474, 140340, 1, 1025, 26988,
    155034, 431644, 808564, 1183264, 1463944, 1604284, 1, 2305,
    102473, 841960, 2992648, 6703438, 11276398, 15626632, 18835200,
    20439484,
]


def mul(a, b, degree):
    """Ordinary truncated polynomial multiplication over Q."""
    out = [Q(0)] * (degree + 1)
    for i, ai in enumerate(a[:degree + 1]):
        if ai:
            for j, bj in enumerate(b[:degree + 1 - i]):
                if bj:
                    out[i + j] += ai * bj
    return out


def formal_log_one(a, degree):
    """log(a) = sum_{r>=1} (-1)^(r+1)(a-1)^r/r, a[0]=1."""
    if a[0] != 1:
        raise AssertionError('formal log requires constant term one')
    q = list(a[:degree + 1])
    q[0] -= 1
    power = [Q(1)] + [Q(0)] * degree
    out = [Q(0)] * (degree + 1)
    for r in range(1, degree + 1):
        power = mul(power, q, degree)
        factor = Q((-1) ** (r + 1), r)
        out = [s + factor * t for s, t in zip(out, power)]
    return out


def source_ode(degree):
    """Solve U'=U+U^2 log U, U(0)=1, coefficient by coefficient.

    At step m only coefficients through degree m are read. No array
    entries or conjectured identities are used in this computation.
    """
    u = [Q(1)]
    for m in range(degree):
        logarithm = formal_log_one(u, m)
        rhs_nonlinear = mul(mul(u, u, m), logarithm, m)
        u.append((u[m] + rhs_nonlinear[m]) / Q(m + 1))
    return u


def source_array(total_degree):
    """B(n,k), on the triangular dependency closure n+k<=total_degree."""
    rows = [[1] * (total_degree + 1)]
    for n in range(total_degree):
        prev = rows[-1]
        rows.append([
            prev[k + 1] + (k + 1) * sum(comb(k, j) * prev[j]
                                        for j in range(k + 1))
            for k in range(total_degree - n)
        ])
    return rows


def require_equal(actual, expected, label):
    if actual != expected:
        raise AssertionError(f'{label}: {actual!r} != {expected!r}')


def main():
    u = source_ode(MAX_DEGREE)
    a = [q * factorial(n) for n, q in enumerate(u)]
    require_equal(a[:len(DISPLAYED_A249833)], DISPLAYED_A249833,
                  'A249833 displayed coefficients')
    print(f'FIDELITY A249833: all {len(DISPLAYED_A249833)} displayed coefficients agree')

    log_u = formal_log_one(u, 10)
    require_equal([q * factorial(n) for n, q in enumerate(log_u)],
                  [0, 1, 1, 3, 12, 64, 424, 3358, 30952, 325488, 3845724],
                  'A249833 displayed log expansion')
    print('FIDELITY A249833: displayed log(U) expansion through degree 10 agrees')
    square = mul(u, u, 9)
    require_equal([q * factorial(n) for n, q in enumerate(square)],
                  [1, 2, 6, 26, 148, 1040, 8688, 84068, 924384, 11381696],
                  'A249833 displayed square expansion')
    nonlinear = mul(mul(u, u, 5), formal_log_one(u, 5), 5)
    require_equal([q * factorial(n) for n, q in enumerate(nonlinear)],
                  [0, 1, 5, 27, 176, 1364], 'A249833 displayed U^2 log(U)')
    print('FIDELITY A249833: displayed U^2 and U^2 log(U) expansions agree')

    rows = source_array(MAX_DEGREE)
    antidiagonals = [rows[n][s - n] for s in range(10) for n in range(s + 1)]
    require_equal(antidiagonals, DISPLAYED_A370382,
                  'A370382 displayed downward antidiagonals')
    print(f'FIDELITY A370382: all {len(antidiagonals)} displayed entries agree')

    require_equal(a[0], 1, 'a(0)')
    for n in range(MAX_DEGREE):
        require_equal(a[n + 1], rows[n][0], f'column identity n={n}')
    print(f'IDENTITY: a(n+1)=B(n,0) for every 0<=n<{MAX_DEGREE}; a(0)=1')
    if not all(q.denominator == 1 for q in a):
        raise AssertionError('nonintegral EGF coefficient in checked window')
    print(f'INTEGRALITY: factorial-normalized source coefficients 0..{MAX_DEGREE} are integers')
    print('INITIAL:', ','.join(str(int(q)) for q in a[:8]))
    print('a(22)=', int(a[22]), sep='')
    print('a(40)=', int(a[40]), sep='')
    print('Arithmetic: fractions.Fraction, integer binomial coefficients and factorials only')
    print('Scope: finite consistency checks; no infinite-identity proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

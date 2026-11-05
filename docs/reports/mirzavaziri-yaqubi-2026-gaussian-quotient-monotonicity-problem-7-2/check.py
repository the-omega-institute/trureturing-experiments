#!/usr/bin/env python3
"""Exact coefficient certificates for arXiv:2610.19249v1, Problem 7.2.

Python standard library only. Coefficient lists are in ascending degree.
Finite family checks support, but do not replace, the unbounded argument.
"""
from itertools import combinations
from math import comb, gcd


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def coeff(a, t):
    return a[t] if 0 <= t < len(a) else 0


def multiply(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for j, x in enumerate(b):
        if x:
            for i, y in enumerate(a):
                c[i + j] += x * y
    return trim(c)


def times_binomial(a, r):
    assert r > 0
    c = list(a) + [0] * r
    for i, x in enumerate(a):
        c[i + r] -= x
    return trim(c)


def divide_binomial(a, r):
    """Exact division by 1-q^r, with a full multiplication check."""
    a = trim(a)
    assert r > 0 and len(a) > r
    b = [0] * (len(a) - r)
    for t in range(len(b)):
        b[t] = a[t] + (b[t - r] if t >= r else 0)
    assert times_binomial(b, r) == a
    return trim(b)


def divide_qinteger(a, r):
    b = divide_binomial(times_binomial(a, 1), r)
    assert multiply(b, [1] * r) == trim(a)
    return b


def gaussian_product(n, k):
    """The literal product in Section 1; every division is checked exactly."""
    assert 0 <= k <= n
    g = [1]
    for i in range(1, k + 1):
        g = divide_binomial(times_binomial(g, n - k + i), i)
    assert len(g) - 1 == k * (n - k)
    assert g == g[::-1] and sum(g) == comb(n, k)
    return g


def gaussian_pascal(n, k):
    """Independent q-Pascal recurrence, as in the D5 Gaussian owner."""
    row = [[1]] + [[0] for _ in range(k)]
    for a in range(n):
        new = [[1]]
        for b in range(1, k + 1):
            if b > a + 1:
                new.append([0])
            else:
                shift = a - b + 1
                lhs, rhs = row[b], [0] * shift + row[b - 1]
                new.append(trim([coeff(lhs, t) + coeff(rhs, t)
                                 for t in range(max(len(lhs), len(rhs)))]))
        row = new
    return row[k]


def negatives(a):
    return [(i, x) for i, x in enumerate(a) if x < 0]


def restricted_counts(weights, bound):
    c = [1] + [0] * bound
    for w in weights:
        for t in range(w, bound + 1):
            c[t] += c[t - w]
    return c


def fidelity():
    # Theorem 4.3(iii): exactly these four negative coefficients for m=4,6.
    for m, expected in [(4, [1, 5, 59, 63]), (6, [1, 7, 147, 153])]:
        g = gaussian_product(6 * m, m)
        assert g == gaussian_pascal(6 * m, m)
        p = divide_qinteger(divide_qinteger(g, 2 * m - 1), 3 * m - 1)
        assert len(p) - 1 == 5 * m * m - 5 * m + 4
        assert negatives(p) == [(t, -1) for t in expected]
        print(f'FIDELITY Theorem 4.3: m={m}, negatives={negatives(p)}')
    # The explicit m=1 identity at the end of the theorem's proof.
    p = divide_qinteger(gaussian_product(6, 1), 2)
    assert p == [1, 0, 1, 0, 1]
    print('FIDELITY P^(1)=1+q^2+q^4')


A = [1,-1,1,0,0,0,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,1,
     0,1,1,0,1,1,1,0,2,0,1,1,1,1,1,1,1,1,1,1]
B = [2,2,2,2,2,1,2,2,1,2,2,1,1,2,1,1,2,1,1,1,1,1,
     1,1,1,1,0,1,1,0,1,1,0,0,1,0,0,1,0,0,0,0]
SMALL_U = [1,0,1,1,1,1,2,2,2,3,3,3,4,4,5,5,6,6,7,7,
           8,9,9,10,11,11,12,13,13,13,14,14,14,15,15,15,15,15]
SMALL_P = [1,-1,1,0,0,0,1,0,0,1,1,-1,2,0,1,0,2,0,1,1,
           2,0,2,1,2,0,3,1,1,1,3,0,2,2,2,0,3,1,1,1,3,0,2,
           1,2,0,2,1,1,0,2,0,1,0,2,-1,1,1,0,0,1,0,0,0,1,-1,1]


def residue_certificates():
    assert len(A) == len(B) == 42
    c = restricted_counts([2, 3, 7], 42 * 50 + 41)
    d = [c[t] - coeff(c, t - 1) for t in range(len(c))]
    assert d[:42] == A
    assert all(d[t + 42] == d[t] + 1 for t in range(len(d) - 42))
    assert all(d[t] == t // 42 + A[t % 42] for t in range(len(d)))
    # E=[42]/([2][3][7]): an exact finite certificate of the period identity.
    e = [1] * 42
    for r in [2, 3, 7]:
        e = divide_qinteger(e, r)
    assert len(e) - 1 == 32 and sum(e) == 1
    assert multiply(multiply(multiply(e, [1] * 2), [1] * 3), [1] * 7) == [1] * 42
    # Multiplying d's rational generating function by (1-q^42) gives E/(1-q).
    # This finite identity and degree(E)<42 establish the affine period law.
    def dt(t):
        return t // 42 + A[t % 42]  # Euclidean residues, also for negative t.
    observed = [dt(r) - dt(r - 28) - dt(r - 29) - dt(r - 30)
                for r in range(-4, 38)]
    assert observed == B and min(B) == 0
    assert min(A[2:]) >= 0 and A[1] == -1
    print(f'RESIDUES A={A}')
    print(f'RESIDUES B={B}')
    print(f'PERIOD certificate: deg(E)={len(e)-1}, E(1)={sum(e)}, checked 50 periods')
    return dt


def small_witness():
    g = gaussian_product(30, 3)
    assert g == gaussian_pascal(30, 3)
    by_subsets = [0] * 82
    for s in combinations(range(30), 3):
        by_subsets[sum(s) - 3] += 1
    assert g == by_subsets and sum(g) == 4060
    u = divide_qinteger(g, 7)
    p = divide_qinteger(u, 10)
    assert u[:38] == SMALL_U and u == u[::-1]
    assert u[37] == u[38] and p == SMALL_P
    assert negatives(p) == [(1, -1), (11, -1), (55, -1), (65, -1)]
    print('SMALL independent subset count=4060, deg(G)=81, deg(U)=75, deg(P)=66')
    print(f'SMALL U coefficients 0..37={u[:38]}')
    print(f'SMALL P all coefficients={p}')
    print(f'SMALL negative coefficients={negatives(p)}')


def family(dt):
    # Each full Gaussian quotient is built from the source product, not from A or B.
    cases = list(range(21)) + [100]
    for m in cases:
        n = 30 + 420 * m
        assert 2 <= 3 <= n // 2 and 2 <= 7 <= 10 <= n and gcd(7, 10) == 1
        # Denominator cyclotomic multiplicities; the only repeated factor is Phi_2.
        den_mult = {h: sum(r % h == 0 for r in [2, 3, 7, 10])
                    for h in range(1, 11)}
        num_mult = {h: sum(r % h == 0 for r in [1, n - 2, n - 1, n])
                    for h in range(1, 11)}
        assert all(num_mult[h] >= den_mult[h] for h in den_mult)
        g = gaussian_product(n, 3)
        u = divide_qinteger(g, 7)
        p = divide_qinteger(u, 10)
        assert len(u) - 1 == 3 * n - 15 and len(p) - 1 == 3 * n - 24
        assert u == u[::-1] and p == p[::-1]
        mid = (len(u) - 1) // 2
        assert mid == 37 + 630 * m
        delta = [u[t] - u[t-1] for t in range(1, mid + 1)]
        assert min(delta[1:]) >= 0  # monotonicity starts at q^1, not q^0.
        assert [u[0], u[1], u[10], u[11]] == [1, 0, 3, 3]
        assert p[1] == p[11] == p[-2] == p[-12] == -1
        assert 11 not in [1, len(p) - 2]
        # A second construction of every first-half difference, from the series.
        def d(t):
            return dt(t) if t >= 0 else 0
        for t in range(2, mid + 1):
            rhs = d(t) - d(t - n + 2) - d(t - n + 1) - d(t - n)
            assert u[t] - u[t - 1] == rhs
        assert delta[n - 3] == delta[n - 2] == 10 * m
        if m:
            assert mid - 41 >= n
            assert [delta[t - 1] for t in range(mid - 41, mid + 1)] == B
            assert all(delta[t + 41] == delta[t - 1] - 2
                       for t in range(n, mid - 41))
        assert times_binomial(p, 10) == times_binomial(u, 1)
        print(f'FAMILY m={m}, N={n}, deg(U)={len(u)-1}, deg(P)={len(p)-1}, '
              f'min(first-half differences from 2)={min(delta[1:])}, '
              f'u10=u11={u[10]}, p11={p[11]}')
    # Residue cancellation in the final block is affine in m: 15m-3(5m)=0.
    assert 630 == 15 * 42 and 210 == 5 * 42
    assert all(x >= 0 for x in B)
    print(f'FAMILY exact full coefficient checks: {len(cases)} cases; m=0..20 and m=100')
    print('SCOPE unbounded family requires the stated cyclotomic/period argument; no Lean proof is claimed')


def main():
    fidelity()
    dt = residue_certificates()
    small_witness()
    family(dt)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

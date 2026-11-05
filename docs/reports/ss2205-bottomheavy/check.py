#!/usr/bin/env python3
"""Exact finite checks for Sagan--Swanson, arXiv:2205.14078v1, Conjecture 7.6.

Python 3 standard library only. Coefficients are integers, in increasing exponent
order; () denotes zero. Finite checks are not a proof of the universal claim.
"""
from itertools import combinations, combinations_with_replacement, product


def trim(p):
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, q):
    return trim((p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0)
                for i in range(max(len(p), len(q))))


def scale(w, p):
    return trim(w * a for a in p)


def mul(p, q):
    if not p or not q:
        return ()
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def mul_u(p, a):
    """Multiply by 1 + q + ... + q^(a-1), using exact sliding sums."""
    if not p or a <= 0:
        return ()
    out = []
    window = 0
    for i in range(len(p) + a - 1):
        if i < len(p):
            window += p[i]
        if 0 <= i - a < len(p):
            window -= p[i - a]
        out.append(window)
    return trim(out)


def bottom_heavy(p):
    if not p:
        return True
    d = len(p) - 1
    return all(p[r] >= p[d - r] for r in range((d + 1) // 2))


def unimodal(p):
    descending = False
    for a, b in zip(p, p[1:]):
        if b < a:
            descending = True
        elif b > a and descending:
            return False
    return True


def tables(limit):
    out = {name: [{0: (1,)}] for name in ('S', 'c', 'S_B', 'c_B')}
    for n in range(1, limit + 1):
        for name, rows in out.items():
            prev = rows[-1]
            row = {}
            for k in range(n + 1):
                a = {'S': k, 'c': n - 1, 'S_B': 2*k + 1, 'c_B': 2*n - 1}[name]
                row[k] = add(prev.get(k - 1, ()), mul_u(prev.get(k, ()), a))
            rows.append(row)
    return out


def symmetric_eval(lengths, m, elementary):
    """Independently enumerate the monomials of e_m or h_m, then convolve."""
    if m < 0:
        return ()
    select = combinations if elementary else combinations_with_replacement
    out = ()
    for indices in select(range(len(lengths)), m):
        term = (1,)
        for i in indices:
            term = mul(term, (1,) * lengths[i])
        out = add(out, term)
    return out


def main():
    data = tables(50)
    sb = (15, 24, 34, 38, 43, 42, 43, 38, 35, 26, 20, 14, 10, 6, 4, 2, 1)
    cb = (21, 36, 51, 60, 70, 74, 79, 78, 79, 74, 71, 62, 56, 44,
          35, 26, 20, 14, 10, 6, 4, 2, 1)
    assert data['S_B'][6][4] == sb
    assert data['c_B'][7][5] == cb
    assert data['S'][4][3] == (3, 2, 1)
    assert not unimodal(sb) and not unimodal(cb)
    assert bottom_heavy(sb) and bottom_heavy(cb)
    print('FIDELITY S_B[6,4] coefficients =', sb)
    print('FIDELITY c_B[7,5] coefficients =', cb)
    print('FIDELITY S[4,3] coefficients = (3, 2, 1)')
    print('FIDELITY both printed type-B examples are nonunimodal and bottom heavy')

    tested = 0
    comparisons = 0
    for name, rows in data.items():
        for n, row in enumerate(rows):
            for k in range(-3, 54):
                p = row.get(k, ())
                assert all(isinstance(a, int) and a >= 0 for a in p)
                assert bottom_heavy(p), (name, n, k, p)
                if k < 0 or k > n:
                    assert p == ()
                tested += 1
                comparisons += len(p) // 2
    print('SETTLEMENT all four recurrences: 0 <= n <= 50, -3 <= k <= 53')
    print('SETTLEMENT polynomial cases =', tested, '; reflected comparisons =', comparisons)

    identities = 0
    for n in range(9):
        for k in range(n + 2):
            for name in data:
                if name == 'S':
                    lengths, elementary = tuple(range(1, k + 1)), False
                elif name == 'c':
                    lengths, elementary = tuple(range(1, n)), True
                elif name == 'S_B':
                    lengths, elementary = tuple(range(1, 2*k + 2, 2)), False
                else:
                    lengths, elementary = tuple(range(1, 2*n, 2)), True
                expected = symmetric_eval(lengths, n - k, elementary)
                assert expected == data[name][n].get(k, ()), (name, n, k)
                identities += 1
    print('BRIDGE h/e expansions: 0 <= n <= 8, 0 <= k <= n+1; identities =', identities)

    products = comparison_cases = 0
    for count in range(5):
        for lengths in product(range(1, 7), repeat=count):
            p = (1,)
            for a in lengths:
                assert mul_u(p, a) == mul(p, (1,) * a)
                p = mul(p, (1,) * a)
            assert p == p[::-1] and unimodal(p) and p[0] == 1
            d = len(p) - 1
            for D in range(d, d + 13):
                for r in range((D + 1) // 2):
                    left = p[r] if r <= d else 0
                    right = p[D-r] if D-r <= d else 0
                    assert left >= right, (lengths, D, r)
                    comparison_cases += 1
            products += 1
    print('STRUCTURE products: 0..4 factors, each length 1..6; products =', products)
    print('STRUCTURE common-left-endpoint comparisons, d <= D <= d+12 =', comparison_cases)

    pool = [(1,)]
    for count in (1, 2):
        for lengths in product(range(1, 5), repeat=count):
            p = (1,)
            for a in lengths:
                p = mul(p, (1,) * a)
            pool.append(p)
    sums = 0
    for p, q in product(pool, repeat=2):
        for a, b in product(range(4), repeat=2):
            assert bottom_heavy(add(scale(a, p), scale(b, q)))
            sums += 1
    print('STRUCTURE weighted sums of two products, weights 0..3; cases =', sums)
    print('LIMIT finite integer checks only; no universal proof or Lean proof is asserted')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

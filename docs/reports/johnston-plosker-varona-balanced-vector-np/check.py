#!/usr/bin/env python3
"""Exact checks for arXiv:2607.06336v1; Python 3 standard library only.
Finite tests do not prove complexity bounds or universal NP-completeness.
"""
from itertools import product, combinations_with_replacement
from math import gcd


def rank(rows, width=None):
    """Integer Gaussian elimination, normalizing each row by its gcd."""
    rows = [list(r) for r in rows if any(r)]
    n = width if width is not None else (len(rows[0]) if rows else 0)
    p = 0
    for j in range(n):
        pivot = next((i for i in range(p, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[p], rows[pivot] = rows[pivot], rows[p]
        a = rows[p][j]
        for i in range(p+1, len(rows)):
            b = rows[i][j]
            if b:
                r = [a*x-b*y for x, y in zip(rows[i], rows[p])]
                g = 0
                for x in r:
                    g = gcd(g, abs(x))
                rows[i] = [x//g for x in r] if g else r
        p += 1
        if p == len(rows):
            break
    return p


def annihilators(v):
    return [e for e in product((-1, 0, 1), repeat=len(v))
            if sum(a*b for a, b in zip(e, v)) == 0]


def balanced_brute(v):
    return bool(v) and rank(annihilators(v), len(v)) == len(v)-1


def coefficients(v):
    """Coefficients of product (1 + X^a + X^(2a))."""
    c = [1]
    for a in v:
        b = [0]*(len(c)+2*a)
        for j, x in enumerate(c):
            for t in range(3):
                b[j+t*a] += x
        c = b
    return c


def remove_factor(c, a):
    """Exact polynomial division by 1+X^a+X^(2a), including remainder."""
    q = []
    for j in range(len(c)-2*a):
        q.append(c[j]-(q[j-a] if j >= a else 0)-(q[j-2*a] if j >= 2*a else 0))
    assert all(x >= 0 for x in q)
    for j in range(len(q), len(c)):
        assert c[j] == sum(q[j-t*a] for t in (1, 2) if 0 <= j-t*a < len(q))
    return q


def gram(v):
    """Exact H=sum_{e.v=0} e e^T by scalar counting and factor deletion.
This recognition oracle does not use completion, extension or saturation.
"""
    c, total = coefficients(v), sum(v)
    h = [[0]*len(v) for _ in v]
    def at(q, s):
        return q[s] if 0 <= s < len(q) else 0
    for i, a in enumerate(v):
        qi = remove_factor(c, a)
        h[i][i] = at(qi, total)+at(qi, total-2*a)
        for j in range(i):
            b, qij = v[j], remove_factor(qi, v[j])
            h[i][j] = h[j][i] = sum(x*y*at(qij, total-a-b-x*a-y*b)
                                    for x in (-1, 1) for y in (-1, 1))
    return h


def balanced_dp(v):
    return bool(v) and rank(gram(v), len(v)) == len(v)-1


def partitions(n, minimum=1):
    if n == 0:
        yield ()
    for a in range(minimum, n+1):
        for rest in partitions(n-a, a):
            yield (a,)+rest


def complete(a, target):
    """Balanced completion, its ternary matrix, prescribed-entry positions."""
    assert a and all(x > 0 for x in a) and target >= 0
    k = max(1, max([target]+list(a)).bit_length())
    m, values, sparse = 1 << k, [], []
    def pos(x):
        values.append(x)
        return len(values)-1
    def equation(terms):
        assert len(terms) == len({j for j, _ in terms})
        assert all(c in (-1, 1) for _, c in terms)
        assert sum(values[j]*c for j, c in terms) == 0
        sparse.append(dict(terms))
    u, up = pos(m), pos(m)
    w = [pos(m+(1 << j)) for j in range(k+1)]
    wp = [pos(m+(1 << j)) for j in range(k)]
    equation([(up, 1), (u, -1)])
    for j in range(k):
        equation([(wp[j], 1), (w[j], -1)])
        equation([(w[j+1], 1), (u, 1), (w[j], -1), (wp[j], -1)])
    equation([(w[k], 1), (u, -1), (up, -1)])
    copies = []
    for _ in range(k):
        q = pos(m)
        copies.append(q)
        equation([(q, 1), (u, -1)])
    original = []
    for x in a:
        shifted = pos(m+x)
        bits = [j for j in range(k) if (x >> j) & 1]
        equation([(shifted, 1), (u, -1)]+[(w[j], -1) for j in bits]
                 +[(copies[j], 1) for j in bits])
        p = pos(x)
        original.append(p)
        equation([(shifted, 1), (u, -1), (p, -1)])
    rows = [[r.get(i, 0) for i in range(len(values))] for r in sparse]
    assert len(rows) == len(values)-1 and rank(rows) == len(values)-1
    assert all(x > target for j, x in enumerate(values) if j not in original)
    assert all(x <= 2*m for x in values)
    assert len(values) == 3*k+3+2*len(a)
    return tuple(values), rows, original


def bounded_two(a, target):
    return any(sum(x*y for x, y in zip(a, c)) == target
               for c in product(range(3), repeat=len(a)))


def subset_sum(a, target):
    return any(sum(x*y for x, y in zip(a, c)) == target
               for c in product(range(2), repeat=len(a)))


def marker_reduce(w, target):
    base = max(5, 2*sum(w)+1, target+1)
    weights = tuple(x for i, a in enumerate(w, 1) for x in (base**i+a, base**i))
    return weights, sum(base**i for i in range(1, len(w)+1))+target


def main():
    expected = [1, 2, 2, 4, 3, 7, 6, 13, 12]
    actual = [sum(balanced_brute(v) for v in partitions(n)) for n in range(1, 10)]
    assert actual == expected, actual
    print('FIDELITY balanced-partition counts, totals 1..9:', actual, flush=True)
    cases = 0
    for d in range(1, 6):
        for v in combinations_with_replacement(range(1, 5), d):
            es = annihilators(v)
            direct = [[sum(e[i]*e[j] for e in es) for j in range(d)] for i in range(d)]
            assert gram(v) == direct, v
            assert balanced_dp(v) == balanced_brute(v), v
            cases += 1
    print('ORACLE explicit ternary rows versus generating-function Gram:', cases, 'cases', flush=True)
    cases = 0
    for d in range(1, 5):
        for v in combinations_with_replacement(range(1, 5), d):
            if balanced_brute(v):
                signed = {sum(x*y for x, y in zip(v, e)) for e in product((-1, 0, 1), repeat=d)}
                for c in range(1, sum(v)+3):
                    assert balanced_brute(v+(c,)) == (c in signed), (v, c)
                    cases += 1
    print('EXTENSION appended entry iff signed sum:', cases, 'cases', flush=True)
    completions = recognitions = 0
    for d in range(1, 4):
        for a in product(range(1, 4), repeat=d):
            for target in range(6):
                v, _, _ = complete(a, target)
                completions += 1
                if d <= 2:
                    c = sum(v)-target
                    assert c > 0
                    assert balanced_dp(v+(c,)) == bounded_two(a, target), (a, target)
                    assert balanced_dp(tuple(sorted(v+(c,), reverse=True))) == bounded_two(a, target)
                    recognitions += 1
    print('COMPLETION ternary rows, exact rank, auxiliary barrier:', completions, 'cases', flush=True)
    print('REDUCTION full final-vector recognition (also sorted):', recognitions, 'cases', flush=True)
    cases = 0
    for d in range(1, 4):
        for w in product(range(1, 4), repeat=d):
            for target in range(sum(w)+2):
                a, t = marker_reduce(w, target)
                assert bounded_two(a, t) == subset_sum(w, target), (w, target)
                cases += 1
    print('NO-CARRY subset sum iff bounded-two subset sum:', cases, 'cases', flush=True)
    for w, t in [((1,), 0), ((1,), 1), ((1,), 2)]:
        a, target = marker_reduce(w, t)
        v, _, _ = complete(a, target)
        final = v+(sum(v)-target,)
        actual = balanced_dp(final)
        assert actual == subset_sum(w, t)
        print('COMPOSITE', w, 'target', t, 'dimension', len(final), 'balanced', actual, flush=True)
    huge = (2**80+1, 2**127-1, 2**64+2**19)
    v, _, _ = complete(huge, 2**100+3)
    print('BINARY SCALE max input bits', max(x.bit_length() for x in huge),
          'completion dimension', len(v), 'max output bits', max(x.bit_length() for x in v), flush=True)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

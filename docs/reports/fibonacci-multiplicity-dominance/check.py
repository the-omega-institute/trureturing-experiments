#!/usr/bin/env python3
"""Exact source-fidelity and finite support checks for Petrov--Scott 8.22.

Only Python's standard library is used. These finite checks are not a proof
of the universal conjecture. B_literal preserves the printed index bound;
B_weight_preserving is explicitly a separate interpretation.
"""
from collections import Counter
from itertools import accumulate
from math import comb


def compositions(n):
    if n == 0:
        yield ()
    else:
        for a in range(1, n + 1):
            for tail in compositions(n - a):
                yield (a,) + tail


def A(k):
    m, p = len(k), len(k) // 2
    if m <= 2:
        return k
    if m == 3:
        return (k[0], k[1] + k[2])
    return (k[0],) + tuple(k[2*j-1] + k[2*j] for j in range(1, p)) + (sum(k[2*p-1:]),)


def B_literal(k):
    m, p = len(k), len(k) // 2
    if m <= 2:
        return k
    if m == 3:
        return (k[0] + k[2], k[1])
    return (1 - p + sum(k[2*j] for j in range(p)),) + tuple(1 + k[2*j-1] for j in range(1, p)) + (k[2*p-1],)


def B_weight_preserving(k):
    b = B_literal(k)
    if len(k) >= 5 and len(k) % 2 == 1:
        return (b[0] + k[-1],) + b[1:]
    return b


def fib(c):
    assert c and all(a > 0 for a in c)
    assert all(a >= 2 for a in c[1:-1])
    if len(c) == 1:
        return (1,) * c[0]
    w = (1,) * (c[0] - 1)
    for j, a in enumerate(c[1:], start=1):
        w += (2,) + (1,) * (a - (1 if j == len(c) - 1 else 2))
    assert sum(w) == sum(c)
    return w


def ribbon_composition(w):
    runs, r = [], 0
    for a in w:
        if a == 1:
            r += 1
        else:
            runs.append(r)
            r = 0
    runs.append(r)
    if len(runs) == 1:
        return (len(w),)
    return tuple(r + (1 if j in (0, len(runs)-1) else 2) for j, r in enumerate(runs))


def partitions(n):
    """Canonical block lists, via adjoining the largest element exactly once."""
    if n == 0:
        yield ()
    else:
        for pi in partitions(n - 1):
            yield pi + ((n,),)
            for j in range(len(pi)):
                yield pi[:j] + (pi[j] + (n,),) + pi[j+1:]


def arcs(pi):
    return tuple((a, b) for block in pi for a, b in zip(block, block[1:]))


def noncrossing(pi):
    es = arcs(pi)
    return not any(a < c < b < d for a, b in es for c, d in es)


def statistics(pi):
    n, es = sum(map(len, pi)), arcs(pi)
    ell, g = Counter(), Counter()
    predecessor = {b: a for a, b in es}
    for i in range(1, n + 1):
        active = sorted(a for a, b in es if a < i <= b)
        ell[len(active)] += 1
        if i in predecessor:
            gamma = active.index(predecessor[i]) + 1
            g[gamma] += 1
    return ell, g


def z(pi):
    ell, g = statistics(pi)
    p = max(ell)
    result = (ell[0],) + tuple(a for k in range(1, p+1) for a in (g[k], ell[k]-g[k]))
    if result[-1] == 0:
        result = result[:-1]
    assert all(a > 0 for a in result)
    assert sum(result) == sum(map(len, pi))
    return result


def check_support(k, b, same_weight):
    u, v = fib(A(k)), fib(b(k))
    if same_weight and sum(v) != sum(k):
        return False
    assert sum(u) == sum(v) == sum(k)
    assert u.count(2) == v.count(2)
    assert len(u) == len(v)
    assert all(x >= y for x, y in zip(accumulate(u), accumulate(v)))
    assert u >= v  # descending lexicographic order, as displayed in Figure 8
    ru, rv = ribbon_composition(u), ribbon_composition(v)
    assert ru == A(k) and rv == b(k)
    return True


def main():
    w = tuple(map(int, '121122112'))
    assert ribbon_composition(w) == (2, 4, 2, 4, 1)
    assert fib((2, 4, 2, 4, 1)) == w
    print('FIDELITY Definition 8.10 ribbon: 121122112 -> (2,4,2,4,1)')
    ell, g = statistics(((1,3,5), (2,9), (4,), (6,7,8)))
    assert ell == Counter({0:1, 1:3, 2:5}) and g == Counter({1:2, 2:3})
    print('FIDELITY Section 8.1 example: ell=(1,3,5), g=(2,3)')
    k = (1,1,1,1,1)
    assert A(k) == (1,2,2) and B_literal(k) == (1,2,1)
    assert B_weight_preserving(k) == (2,2,1)
    print('SOURCE_BOUNDARY k=(1,1,1,1,1): A=(1,2,2); B_literal=(1,2,1), weight=4; B_weight_preserving=(2,2,1), weight=5')
    matrix5 = None
    for n in range(1, 10):
        counts = Counter(z(pi) for pi in partitions(n) if noncrossing(pi))
        cat = comb(2*n, n) // (n+1)
        assert sum(counts.values()) == cat
        assert set(counts) == set(compositions(n))
        intended, literal = Counter(), Counter()
        for k, count in counts.items():
            u, v = fib(A(k)), fib(B_weight_preserving(k))
            intended[u,v] += count
            v0 = fib(B_literal(k))
            if sum(v0) == n:
                literal[u,v0] += count
        assert sum(intended.values()) == cat
        assert len(intended) == 2**(n-1)
        for mat in (intended, literal):
            for (u,v), count in mat.items():
                assert count > 0 and u >= v and u.count(2) == v.count(2)
                assert all(a >= b for a,b in zip(accumulate(u),accumulate(v)))
        print(f'NC n={n}: Catalan={cat}, compositions={len(counts)}, intended_support={len(intended)}, literal_same_weight_support={len(literal)}, literal_sum={sum(literal.values())}')
        if n == 5:
            matrix5 = intended
            literal5 = literal
    order = [tuple(map(int,s)) for s in ['221','212','2111','122','1211','1121','1112','11111']]
    expected = [
        [1,0,0,1,0,0,0,0], [0,2,0,2,0,0,0,0],
        [0,0,1,0,3,3,1,0], [0,0,0,2,0,0,0,0],
        [0,0,0,0,4,6,2,0], [0,0,0,0,0,6,3,0],
        [0,0,0,0,0,0,4,0], [0,0,0,0,0,0,0,1]]
    assert [[matrix5[u,v] for v in order] for u in order] == expected
    differences = [(u,v,matrix5[u,v],literal5[u,v]) for u in order for v in order if matrix5[u,v] != literal5[u,v]]
    assert differences == [((2,2,1),(1,2,2),1,0)]
    print('FIDELITY Figure 8: all 64 entries match the explicitly weight-preserving interpretation, sum=42')
    print('SOURCE_BOUNDARY Figure 8 versus literal block: only (221,122) differs: 1 versus 0; literal sum=41')
    total, lit = 0, 0
    for n in range(1, 15):
        images = set()
        for k in compositions(n):
            assert check_support(k, B_weight_preserving, False)
            assert (fib(A(k)),fib(B_weight_preserving(k))) not in images
            images.add((fib(A(k)),fib(B_weight_preserving(k))))
            lit += int(check_support(k, B_literal, True))
            total += 1
    assert total == 2**14 - 1
    print(f'SUPPORT all positive compositions of 1<=n<=14: intended={total}, literal_same_weight={lit}; no dominance, hike, or order violation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for ABDFR arXiv:2607.26284v1, Question 3.10.

No numerical library or floating point is used. Finite checks do not prove
universal regularity, automaticity or failure of a uniform proximity constant.
"""
from collections import deque
from itertools import product, zip_longest


def unit(d):
    return ((0,) * d, 0)


def shift(v, h):
    d = len(v)
    return tuple(v[(i - h) % d] for i in range(d))


def mul(x, y):
    v, h = x
    w, k = y
    z = shift(w, h)
    return (tuple(a + b for a, b in zip(v, z)), h + k)


def inv(x):
    v, h = x
    return (shift(tuple(-a for a in v), -h), -h)


def letter(d, g, sign):
    if g == d:
        return ((0,) * d, sign)
    return (tuple(sign if i == g else 0 for i in range(d)), 0)


def prefixes(d, word):
    out = [unit(d)]
    for g, sign in word:
        out.append(mul(out[-1], letter(d, g, sign)))
    return out


def evaluate(d, word):
    return prefixes(d, word)[-1]


def l1(x, y):
    return sum(abs(a - b) for a, b in zip(x[0], y[0])) + abs(x[1] - y[1])


def alphabet(d):
    return tuple(product(range(d + 1), (1, -1)))


def block_order(d):
    return (1, 0) + tuple(range(2, d)) + (d,)


def normal(d, x):
    v, h = x
    exps = v + (h,)
    return tuple((g, 1 if exps[g] >= 0 else -1)
                 for g in block_order(d) for _ in range(abs(exps[g])))


def block_predicate(d, word):
    rank = {g: i for i, g in enumerate(block_order(d))}
    for a, b in zip(word, word[1:]):
        if rank[a[0]] > rank[b[0]] or (a[0] == b[0] and a[1] != b[1]):
            return False
    return True


def in_m(d, x):
    v, h = x
    return h >= 0 and all(a >= 0 for a in v) and all(
        a == 0 for i, a in enumerate(v) if h < d - 1 and i > h)


def dfa(d, word, positive=False):
    """Finite states: dead, or (block rank, sign, max spatial index, capped t count).

    For L, the last two fields stay zero. For J, reject negative signs and
    accept precisely when max spatial index <= capped count of t.
    """
    ranks = {g: i for i, g in enumerate(block_order(d))}
    state = (-1, 0, 0, 0)
    for g, sign in word:
        last, old_sign, high, height = state
        r = ranks[g]
        if r < last or (r == last and sign != old_sign):
            return False
        if positive and sign != 1:
            return False
        if positive:
            if g == d:
                height = min(d - 1, height + 1)
            else:
                high = max(high, g)
        state = (r, sign, high, height)
    return not positive or state[2] <= state[3]


def fidelity():
    # Source Remark 3.3: Z^2, L all words, Q1=L, Q2={a^i b^j}.
    for n in range(1, 33):
        w = ((1, 1),) * n + ((0, 1),) * n
        q = ((0, 1),) * n + ((1, 1),) * n
        # These paths stay at height zero, so they are ordinary Z^2 paths.
        assert evaluate(2, w) == evaluate(2, q) == ((n, n), 0)
        p = prefixes(2, w)[n]
        assert min(l1(p, z) for z in prefixes(2, q)) == n
        # Q1 takes w'=w and the identity reparametrisation: k=0,c=1.
        assert all(l1(x, y) == 0 for x, y in zip(prefixes(2, w), prefixes(2, w)))
    print('FIDELITY Remark 3.3: 32 exact Z^2 prefix obstructions; Q1 identity k=0,c=1')


def metric_and_group(d):
    points = [(v, h) for v in product(range(-1, 2), repeat=d) for h in range(-1, 2)]
    letters = alphabet(d)
    edge_checks = 0
    for p in points:
        assert mul(p, inv(p)) == mul(inv(p), p) == unit(d)
        neighbours = {mul(p, letter(d, g, s)) for g, s in letters}
        grid = set()
        for i in range(d):
            for s in (1, -1):
                v = list(p[0]); v[i] += s
                grid.add((tuple(v), p[1]))
        grid.update([(p[0], p[1] - 1), (p[0], p[1] + 1)])
        assert neighbours == grid
        assert all(l1(p, q) == 1 for q in neighbours)
        edge_checks += len(neighbours)
    # Independent BFS via group multiplication, compared with the l1 formula.
    seen = {unit(d): 0}; queue = deque([unit(d)])
    while queue:
        p = queue.popleft()
        if seen[p] == 4:
            continue
        for g, s in letters:
            q = mul(p, letter(d, g, s))
            if q not in seen:
                seen[q] = seen[p] + 1; queue.append(q)
    assert all(dist == l1(unit(d), p) for p, dist in seen.items())
    # Construct a path of l1 length for every pair in the stated box.
    for p in points:
        for q in points:
            word = []
            for i in range(d):
                delta = q[0][i] - p[0][i]
                g = (i - p[1]) % d
                word.extend([(g, 1 if delta >= 0 else -1)] * abs(delta))
            delta = q[1] - p[1]
            word.extend([(d, 1 if delta >= 0 else -1)] * abs(delta))
            assert mul(p, evaluate(d, word)) == q
            assert len(word) == l1(p, q)
    print(f'METRIC d={d}: {edge_checks} grid edges; BFS radius 4 has {len(seen)} vertices; {len(points)**2} constructed pair paths')


def monoid_and_regular(d):
    # Positive steps cannot leave and later re-enter these coordinate bounds.
    origin = unit(d); reachable = {origin}; queue = deque([origin])
    while queue:
        p = queue.popleft()
        for g in (0, d):
            q = mul(p, letter(d, g, 1))
            if q[1] <= 5 and all(a <= 2 for a in q[0]) and q not in reachable:
                reachable.add(q); queue.append(q)
    checked = 0
    for v in product(range(-2, 3), repeat=d):
        for h in range(-2, 6):
            p = (v, h)
            assert (p in reachable) == in_m(d, p)
            checked += 1
    word_checks = 0
    for length in range(6):
        for word in product(alphabet(d), repeat=length):
            expected_l = block_predicate(d, word)
            assert dfa(d, word) == expected_l
            assert dfa(d, word, True) == (expected_l and in_m(d, evaluate(d, word)))
            word_checks += 1
    # Roundtrip all canonical forms in a box, including negative heights.
    forms = 0
    for v in product(range(-2, 3), repeat=d):
        for h in range(-3, 4):
            p = (v, h); w = normal(d, p)
            assert evaluate(d, w) == p and dfa(d, w)
            assert len(w) == l1(origin, p)
            assert dfa(d, w, True) == in_m(d, p)
            if in_m(d, p):
                # Proposition 3.8's representing language J uses identity f.
                assert all(l1(x, y) == 0 for x, y in zip(prefixes(d, w), prefixes(d, w)))
            forms += 1
    print(f'MONOID/REGULAR d={d}: {checked} box endpoints ({len(reachable)} positive reachable states); {word_checks} words length <=5; {forms} normal forms')


def fellow_traveller(d):
    comparisons = 0; maximum = 0
    for v in product(range(-2, 3), repeat=d):
        for h in range(-3, 4):
            p = (v, h); wp = prefixes(d, normal(d, p))
            for g, s in alphabet(d):
                q = mul(p, letter(d, g, s)); wq = prefixes(d, normal(d, q))
                # After a word ends its path is padded by its endpoint.
                distances = [l1(wp[min(i, len(wp)-1)], wq[min(i, len(wq)-1)])
                             for i in range(max(len(wp), len(wq)))]
                assert max(distances) <= 2
                maximum = max(maximum, max(distances)); comparisons += 1
    print(f'AUTOMATIC d={d}: {comparisons} signed right-adjacent normal-form pairs; max synchronous distance={maximum}')


def witnesses():
    cases = 0
    for d in range(2, 11):
        for n in range(1, 101):
            w = ((1, 1),) * n + ((0, 1),) * n + ((d, 1),)
            u = ((0, 1),) * n + ((d, 1),) + ((0, 1),) * n
            endpoint = (tuple(n if i in (0, 1) else 0 for i in range(d)), 1)
            assert evaluate(d, w) == evaluate(d, u) == endpoint
            assert dfa(d, w, True)
            assert len(w) == len(u) == 2*n + 1
            p = prefixes(d, w)[n]
            pu = prefixes(d, u)
            assert [l1(p, z) for z in pu[:n+1]] == [r+n for r in range(n+1)]
            assert [l1(p, z) for z in pu[n+1:]] == [n+abs(n-s)+1 for s in range(n+1)]
            assert min(l1(p, z) for z in pu) == n
            cases += 1
    # Exhaust all T-words up to length 11: positivity makes them geodesic.
    enum = 0; unique = 0
    for d in (2, 3):
        buckets = {}
        for length in range(12):
            for word in product(((0, 1), (d, 1)), repeat=length):
                p = evaluate(d, word)
                assert l1(unit(d), p) == length
                buckets.setdefault(p, []).append(word); enum += 1
        for n in range(1, 6):
            endpoint = (tuple(n if i in (0, 1) else 0 for i in range(d)), 1)
            assert buckets[endpoint] == [((0, 1),)*n + ((d, 1),) + ((0, 1),)*n]
            unique += 1
    print(f'REFUTATION: {cases} family instances d=2..10,N=1..100; exact minimum distance=N; both lengths=2N+1')
    print(f'POSITIVE WORDS: {enum} T-words length <=11 are geodesic; {unique} unique-witness endpoint buckets')
    print('INDEXING: witness lengths agree, so source printed |w\u2032| and map-domain |w| test identical indices')


def main():
    fidelity()
    for d in (2, 3):
        metric_and_group(d)
        monoid_and_regular(d)
        fellow_traveller(d)
    witnesses()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Exact finite checks for Pastorek Question 2; no external graph library."""
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, isqrt


def pair(a, b):
    return (a, b) if a < b else (b, a)


def constraint_rank(n, domain, image):
    edges = list(combinations(range(n), 2))
    parent = {e: e for e in edges}

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    f = dict(zip(domain, image))
    rank = changed = 0
    for a, b in combinations(domain, 2):
        e, h = pair(a, b), pair(f[a], f[b])
        changed += e != h
        x, y = root(e), root(h)
        if x != y:
            parent[x] = y
            rank += 1
    return rank, changed


def connected(n, edges):
    seen, todo = {0}, [0]
    while todo:
        a = todo.pop()
        for b in range(n):
            if a != b and pair(a, b) in edges and b not in seen:
                seen.add(b)
                todo.append(b)
    return len(seen) == n


def depth(n, edges):
    for r in range(n, 0, -1):
        for domain in combinations(range(n), r):
            for image in permutations(range(n), r):
                if all(a == b for a, b in zip(domain, image)):
                    continue
                f = dict(zip(domain, image))
                if all((pair(a, b) in edges) == (pair(f[a], f[b]) in edges)
                       for a, b in combinations(domain, 2)):
                    return n - r
    raise AssertionError('Nontrivial singleton map missing')


def source_fidelity():
    # Appendix A: the source's 20-vertex, 30-edge witness and rank-17 map.
    names = ['x', 'u', 'w', 'v', 'y', 'rho', 's1', 's2', 's3'] + [f'e{i}' for i in range(11)]
    idx = {v: i for i, v in enumerate(names)}
    text = 'x:u u:w w:v v:y rho:x rho:y s1:u s1:y s2:v s2:x s3:w e0:e3 e0:e4 e0:s3 e1:e4 e1:e7 e1:s2 e2:e3 e2:e7 e2:e9 e3:rho e4:e5 e5:e6 e5:e9 e6:e10 e6:s1 e7:e10 e8:e9 e8:e10 e8:s3'
    edges = {pair(idx[a], idx[b]) for a, b in (s.split(':') for s in text.split())}
    assert len(edges) == 30 and connected(20, edges)
    assert all(sum(v in e for e in edges) == 3 for v in range(20))
    assert not any(all(pair(a, b) in edges for a, b in zip(c, c[1:] + c[:1]))
                   for r in (3, 4) for c in permutations(range(20), r))
    assert all(pair(idx[a], idx[b]) in edges for a, b in
               zip(['x', 'u', 's1', 'y', 'rho'], ['u', 's1', 'y', 'rho', 'x']))
    f = {v: v for v in range(20) if names[v] not in ('s1', 's2', 's3')}
    for a, b in [('x', 'y'), ('u', 'v')]:
        f[idx[a]], f[idx[b]] = idx[b], idx[a]
    assert len(f) == len(set(f.values())) == 17
    assert all((pair(a, b) in edges) == (pair(f[a], f[b]) in edges)
               for a, b in combinations(f, 2))
    print('SOURCE_FIDELITY: Appendix A witness: n=20, m=30, connected, cubic, girth=5; stated map rank=17')


def small_graphs():
    count = conn = 0
    for n in range(2, 6):
        pairs = list(combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            edges = {e for i, e in enumerate(pairs) if bits & (1 << i)}
            d = depth(n, edges)
            count += 1
            assert 2 * d < n
            for u, v in pairs:
                delta = {w for w in range(n) if w not in (u, v)
                         and ((pair(u, w) in edges) != (pair(v, w) in edges))}
                f = {w: w for w in range(n) if w not in delta}
                f[u], f[v] = v, u
                assert all((pair(a, b) in edges) == (pair(f[a], f[b]) in edges)
                           for a, b in combinations(f, 2))
                assert d <= len(delta)
            if connected(n, edges):
                conn += 1
                assert n * d < 4 * len(edges)
    print(f'SMALL_GRAPHS: all {count} labelled graphs on 2..5 vertices; {conn} connected; depth and transposition/averaging bounds verified')


def injections():
    count = equality_checks = 0
    for n in range(2, 8):
        pairs = list(combinations(range(n), 2))
        for r in range(1, n + 1):
            for domain in combinations(range(n), r):
                for image in permutations(range(n), r):
                    s = sum(a != b for a, b in zip(domain, image))
                    if not s:
                        continue
                    rank, changed = constraint_rank(n, domain, image)
                    assert 2 * changed >= s * (2 * r - s - 2)
                    assert 2 * rank >= changed
                    assert 4 * rank >= s * (2 * r - s - 2)
                    count += 1
                    if n <= 4:
                        f = dict(zip(domain, image))
                        good = 0
                        for bits in range(1 << len(pairs)):
                            edges = {e for i, e in enumerate(pairs) if bits & (1 << i)}
                            good += all((pair(a, b) in edges) == (pair(f[a], f[b]) in edges)
                                        for a, b in combinations(domain, 2))
                        assert good == 2 ** (len(pairs) - rank)
                        equality_checks += 1
    print(f'PARTIAL_INJECTIONS: {count} nonidentity maps on 2..7 vertices; independent-constraint bound verified')
    print(f'EDGE_BIT_PROBABILITIES: {equality_checks} maps on 2..4 vertices; exact satisfying-graph count = 2^(binom(n,2)-rank)')


def scalar_checks():
    for g in range(1, 100001):
        h = (7 + isqrt(1 + 48 * g)) // 2
        b = 2 * h - 2
        # 6 + sqrt(12g+24) <= b, by sign and squaring; equality allowed at g=1.
        assert b >= 6 and 12 * g + 24 <= (b - 6) ** 2
    assert 6 + isqrt(36) == 12 == 2 * ((7 + isqrt(49)) // 2) - 2
    # For g>=3, sqrt(48g+1) >= 3 + sqrt(12g+24):
    # after two safe squarings this reduces to (36g-32)^2 >= 36(12g+24).
    assert 36 * 3 - 32 > 0
    assert (36 * 3 - 32) ** 2 > 36 * (12 * 3 + 24)
    # Difference polynomial has positive coefficients after g=t+3.
    assert (1296, 5040, 3616) == (1296, 2 * 1296 * 3 - 2736, 1296 * 9 - 2736 * 3 + 160)
    print('GENUS_SCALARS: exact radical comparison for g=1..100000; g=1 equality recorded; general g>=3 polynomial certificate')

    n0 = 65536
    # Sum of omission-domain counts <= 2^(n/10), from the weighted binomial theorem.
    assert 33 ** 100 * 32 < 2 ** 10 * 32 ** 100
    # n^2 <= 2^floor(n/64): 64 bases and the n->n+64 induction inequality.
    for n in range(n0, n0 + 64):
        assert n * n <= 1 << (n // 64)
    assert n0 * n0 - 128 * n0 - 4096 > 0
    # The exponent inequalities repeat modulo 960, improving on each period.
    for n in range(n0, n0 + 960):
        assert n // 5 - n // 64 >= n // 6
        assert (n + 9) // 10 + 3 <= n // 6
        assert 5 * (99 * n - 200) >= 400 * n
    assert 960 // 5 - 960 // 64 - 960 // 6 > 0
    for n in [65536, 65537, 65600, 100000, 1000000]:
        omission_count = sum(comb(n, k) for k in range(n // 100 + 1)) if n <= 100000 else None
        if omission_count is not None:
            assert omission_count <= 2 ** ((n + 9) // 10)
        bad_partial = Fraction(2 ** ((n + 9) // 10), 2 ** (n // 6) - 1)
        disconnected = Fraction(n, 2 ** (n // 2) - n)
        assert bad_partial < Fraction(1, 4)
        assert disconnected < Fraction(1, 4)
        assert bad_partial + disconnected < 1
    print('EXISTENCE_BUDGET: exact entropy, 64-step induction, 960-residue exponent certificates; union bound <1 for n=65536,65537,65600,100000,1000000')
    print('SCOPE: finite checks support the general proof outline; they do not constitute a Lean proof or enumerate the large existence graphs')


def main():
    source_fidelity()
    small_graphs()
    injections()
    scalar_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()

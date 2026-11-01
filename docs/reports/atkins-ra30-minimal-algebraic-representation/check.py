#!/usr/bin/env python3
"""Exact independent checks for Atkins arXiv:2605.14190v1, Question 7.5.
Python standard library only. No graph data, eigenvalues, or automorphisms
are downloaded. Finite checks support the proposed general argument;
they do not prove its universal structural or spectral implications.
"""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations, product
from math import isqrt


def hoffman_singleton():
    """Five pentagons P(i,j), five pentagrams Q(i,j), arithmetic modulo 5.
    Cross edges P(i,j)--Q(k,i*k+j); an explicit 50-vertex model.
    """
    vertices = list(product(range(2), range(5), range(5)))
    index = {v: i for i, v in enumerate(vertices)}
    adj = [set() for _ in vertices]

    def edge(u, v):
        x, y = index[u], index[v]
        adj[x].add(y)
        adj[y].add(x)

    for i, j in product(range(5), repeat=2):
        edge((0, i, j), (0, i, (j + 1) % 5))
        edge((1, i, j), (1, i, (j + 2) % 5))
        for k in range(5):
            edge((0, i, j), (1, k, (i * k + j) % 5))
    assert len(adj) == 50 and all(len(a) == 7 for a in adj)
    for x, y in combinations(range(50), 2):
        assert len(adj[x] & adj[y]) == (0 if y in adj[x] else 1)
    return adj


def distance_matrix(adj):
    distances = []
    for x in range(len(adj)):
        d = {x: 0}
        q = deque([x])
        while q:
            y = q.popleft()
            for z in sorted(adj[y]):
                if z not in d:
                    d[z] = d[y] + 1
                    q.append(z)
        assert len(d) == len(adj)
        distances.append([d[y] for y in range(len(adj))])
    return distances


def automorphism(colour, initial):
    """Find one exact colour-preserving permutation extending initial.
    Backtracking uses all-different and ordered-pair colour constraints.
    A completed answer is independently checked on every ordered pair.
    """
    n = len(colour)
    neighbours = [[{v for v in range(n) if colour[u][v] == h}
                   for h in range(4)] for u in range(n)]
    domains = [set(range(n)) for _ in range(n)]
    for x, y in initial.items():
        domains[x] = {y}

    def search(ds):
        done = set()
        while True:
            todo = next((i for i in range(n) if len(ds[i]) == 1 and i not in done), None)
            if todo is None:
                break
            done.add(todo)
            target = next(iter(ds[todo]))
            for i in range(n):
                if i == todo:
                    continue
                ds[i] &= neighbours[target][colour[todo][i]]
                if not ds[i]:
                    return None
        if len(done) == n:
            return tuple(next(iter(d)) for d in ds)
        i = min((i for i in range(n) if len(ds[i]) > 1), key=lambda i: (len(ds[i]), i))
        for y in sorted(ds[i]):
            child = [set(d) for d in ds]
            child[i] = {y}
            answer = search(child)
            if answer is not None:
                return answer
        return None

    p = search(domains)
    assert p is not None and len(set(p)) == n
    assert all(p[x] == y for x, y in initial.items())
    assert all(colour[p[x]][p[y]] == colour[x][y]
               for x, y in product(range(n), repeat=2))
    return p


def pair_orbit(pair, generators):
    reached = {pair}
    q = deque([pair])
    while q:
        x, y = q.popleft()
        for g in generators:
            z = (g[x], g[y])
            if z not in reached:
                reached.add(z)
                q.append(z)
    return reached


def fidelity():
    hs = hoffman_singleton()
    base = sorted(set(range(50)) - {0} - hs[0])
    assert len(base) == 42
    reindex = {v: i for i, v in enumerate(base)}
    adj = [{reindex[w] for w in hs[v] if w in reindex} for v in base]
    dist = distance_matrix(adj)
    assert max(map(max, dist)) == 3
    assert all(Counter(row) == Counter({0: 1, 1: 6, 2: 30, 3: 5}) for row in dist)
    # Intersection array, checked at every ordered pair, not just one root.
    expected_b = [6, 5, 1, 0]
    expected_c = [0, 1, 1, 6]
    for x, y in product(range(42), repeat=2):
        h = dist[x][y]
        assert sum(dist[x][z] == h + 1 for z in adj[y]) == expected_b[h]
        assert sum(dist[x][z] == h - 1 for z in adj[y]) == expected_c[h]
    print('FIDELITY: second subconstituent has 42 vertices, intersection array {6,5,1;1,1,6}')
    # 0=identity, 1=a, 2=b, 3=c, exactly as the Lean conventions.
    label = [0, 2, 1, 3]
    colour = [[label[d] for d in row] for row in dist]
    table = {
        (0, 0): {0}, (0, 1): {1}, (0, 2): {2}, (0, 3): {3},
        (1, 0): {1}, (2, 0): {2}, (3, 0): {3},
        (1, 1): {0, 1, 2, 3}, (1, 2): {1, 2, 3}, (2, 1): {1, 2, 3},
        (1, 3): {1, 2}, (3, 1): {1, 2}, (2, 2): {0, 1},
        (2, 3): {1}, (3, 2): {1}, (3, 3): {0, 3},
    }
    for x, y in product(range(42), repeat=2):
        realised = {(colour[x][z], colour[z][y]) for z in range(42)}
        required = {(i, j) for i, j in product(range(4), repeat=2)
                    if colour[x][y] in table[i, j]}
        assert realised == required
    print('FIDELITY: all 28224 ordered-base/product conditions match the full 30_65 composition table')
    classes = {tuple(y for y in range(42) if colour[x][y] in {0, 3}) for x in range(42)}
    assert len(classes) == 7 and all(len(c) == 6 for c in classes)
    for c, d in combinations(sorted(classes), 2):
        assert all(sum(y in adj[x] for y in d) == 1 for x in c)
    print('FIDELITY: seven c-cliques of size six, b-perfect-matchings between each pair')
    # Exact automorphisms witness all four ordered-pair orbitals.
    pairs = [set((x, y) for x, y in product(range(42), repeat=2)
                 if colour[x][y] == h) for h in range(4)]
    representatives = [min(p) for p in pairs]
    generators = []
    while True:
        missing = None
        for h in range(4):
            orbit = pair_orbit(representatives[h], generators)
            assert orbit <= pairs[h]
            if orbit != pairs[h]:
                missing = (representatives[h], min(pairs[h] - orbit))
                break
        if missing is None:
            break
        (x, y), (u, v) = missing
        generators.append(automorphism(colour, {x: u, y: v}))
    sizes = [len(pair_orbit(p, generators)) for p in representatives]
    assert sizes == [42, 1260, 252, 210]
    print(f'FIDELITY: {len(generators)} verified colour automorphisms generate ordered-pair orbits {sizes}')
    # T^2 = 6I + D and fibre-zero-sum quadratic T^2 + T - 6I = 0.
    for x, y in product(range(42), repeat=2):
        common = len(adj[x] & adj[y])
        assert common == (6 if x == y else int(colour[x][y] == 1))
    # Each e_u-e_v within a fibre lies in W and is annihilated by T^2+T-6I.
    for c in sorted(classes):
        u = c[0]
        for v in c[1:]:
            for x in range(42):
                coeff = lambda y: len(adj[x] & adj[y]) + int(y in adj[x]) - 6 * int(x == y)
                assert coeff(u) == coeff(v)
    print('FIDELITY: T^2 = 6I + D; T|W has annihilating polynomial X^2+X-6')


def feasibility_checks():
    # Entire small-r list arising from the general argument's t interval.
    actual = []
    for r in range(3, 6):
        for t in range(2, r - 1):
            mu = Fraction(t * t - 1, r - 1 - t)
            k = t * (t + mu)
            d = (r - 1) * (k + 1)
            mplus = d * (t + mu) / (2 * t + mu)
            actual.append((r, t, mu, k, d, mplus))
            assert mu.denominator != 1 or mplus.denominator != 1
            print(f'EXCLUSION: r={r}, t={t}, mu={mu}, k={k}, d={d}, m_plus={mplus}')
    assert actual == [
        (4, 2, Fraction(3), Fraction(10), Fraction(33), Fraction(165, 7)),
        (5, 2, Fraction(3, 2), Fraction(7), Fraction(32), Fraction(224, 11)),
        (5, 3, Fraction(8), Fraction(33), Fraction(136), Fraction(748, 7)),
    ]
    print('EXCLUSION: r=3 has no integer t in 2 <= t <= r-2')
    # Independent exact cross-check of every possible parameter below 42.
    # k>=r and N=r(k+1) follow from the unbounded structural argument.
    survivors = []
    checked = 0
    for r in range(3, 42):
        for k in range(r, 42):
            if r * (k + 1) >= 42 or (k - 1) % (r - 1):
                continue
            mu = (k - 1) // (r - 1)
            checked += 1
            disc = mu * mu + 4 * k
            s = isqrt(disc)
            if s * s != disc or (s - mu) % 2:
                continue
            t = (s - mu) // 2
            d = (r - 1) * (k + 1)
            mp = Fraction(d * (t + mu), 2 * t + mu)
            if mp.denominator == 1:
                survivors.append((r, k, mu, t))
    assert survivors == []
    print(f'CHECK: {checked} integral count-parameter sets below 42, no spectral survivor')
    # The attained bound satisfies all conditions, with roots 2,-3.
    r, k, mu, t = 6, 6, 1, 2
    d = (r - 1) * (k + 1)
    mp = Fraction(d * (t + mu), 2 * t + mu)
    mm = d - mp
    assert (r * (k + 1), d, mp, mm) == (42, 35, 21, 14)
    assert t * (t + mu) == k and mp * t == mm * (t + mu)
    print('ATTAINMENT PARAMETERS: N=42, r=k=6, mu=1, W multiplicities 21 and 14')


if __name__ == '__main__':
    fidelity()
    feasibility_checks()
    print('ALL_CHECKS_OK')
